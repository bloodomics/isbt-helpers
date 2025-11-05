"""
Script to annotate variants with rsIDs from dbSNP using GRCh38 coordinates.

This script:
1. Fetches variants from the blood group database
2. Queries NCBI Variation API for rsIDs using chromosome, position, ref, and alt alleles
3. Updates the database with PATCH requests (only for variants without existing rsIDs)

Rate Limit: NCBI allows 3 requests per second without API key (script uses 2 req/sec to be safe)

Author: Nick Gleadall
Date: November 2025
"""

import requests
import json
import pandas as pd
import time
import logging
import argparse
from functools import wraps
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# Retry decorator
def retry_with_backoff(retries=3, backoff_factor=0.3):
    """
    A decorator that implements a retry mechanism with exponential backoff
    while preserving the authenticated session.

    :param retries: Number of retries before giving up
    :param backoff_factor: Factor to increase delay between retries
    :return: Decorated function
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            session = kwargs.get("session")
            if session is None:
                raise ValueError("Session must be provided")

            # Create a new session with the same cookies and headers
            retry_session = requests.Session()
            retry_session.cookies.update(session.cookies)
            retry_session.headers.update(session.headers)

            # Set up retry strategy
            retry_strategy = Retry(
                total=retries,
                backoff_factor=backoff_factor,
                status_forcelist=[429, 500, 502, 503, 504],
                allowed_methods=["HEAD", "GET", "OPTIONS", "POST", "PATCH"],
            )
            retry_session.mount("https://", HTTPAdapter(max_retries=retry_strategy))

            # Replace the original session with the retry-enabled session
            kwargs["session"] = retry_session

            for attempt in range(retries + 1):
                try:
                    return func(*args, **kwargs)
                except requests.exceptions.RequestException as e:
                    if attempt == retries:
                        logger.error(f"All retry attempts failed: {str(e)}")
                        raise
                    wait_time = backoff_factor * (2**attempt)
                    logger.warning(
                        f"Attempt {attempt + 1} failed. Retrying in {wait_time:.2f} seconds..."
                    )
                    time.sleep(wait_time)

        return wrapper

    return decorator


def login(lead_url, email, password):
    """
    Authenticate with the blood group database.

    :param lead_url: Base URL of the API
    :param email: User's email for authentication
    :param password: User's password for authentication
    :return: Authenticated session object
    """
    login_url = f"{lead_url}/auth/login"
    login_data = {"email": email, "password": password}

    try:
        session = requests.Session()
        login_response = session.post(login_url, data=login_data)
        login_response.raise_for_status()

        logger.info("Successfully logged in")
        return session
    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to login: {str(e)}")
        if login_response is not None:
            logger.error(f"Response content: {login_response.text}")
        raise


# GRCh38/hg38 chromosome to RefSeq accession mapping
GRCH38_CHROMOSOMES = {
    '1': 'NC_000001.11', '2': 'NC_000002.12', '3': 'NC_000003.12',
    '4': 'NC_000004.12', '5': 'NC_000005.10', '6': 'NC_000006.12',
    '7': 'NC_000007.14', '8': 'NC_000008.11', '9': 'NC_000009.12',
    '10': 'NC_000010.11', '11': 'NC_000011.10', '12': 'NC_000012.12',
    '13': 'NC_000013.11', '14': 'NC_000014.9', '15': 'NC_000015.10',
    '16': 'NC_000016.10', '17': 'NC_000017.11', '18': 'NC_000018.10',
    '19': 'NC_000019.10', '20': 'NC_000020.11', '21': 'NC_000021.9',
    '22': 'NC_000022.11', 'X': 'NC_000023.11', 'Y': 'NC_000024.10',
    'MT': 'NC_012920.1', 'M': 'NC_012920.1'
}


@retry_with_backoff(retries=5, backoff_factor=0.5)
def get_rsid(chromosome, position, ref, alt, session=None):
    """
    Fetch the rsID for a given genomic variant using GRCh38 coordinates.

    :param chromosome: Chromosome number (1-22, X, Y, MT)
    :param position: Genomic position (1-based)
    :param ref: Reference allele
    :param alt: Alternate allele
    :param session: Requests session object
    :return: rsID if found, None otherwise
    """
    try:
        # Convert chromosome to RefSeq accession
        chrom_accession = GRCH38_CHROMOSOMES.get(str(chromosome))
        if not chrom_accession:
            logger.warning(f"Unknown chromosome: {chromosome}")
            return None
        
        # NCBI SPDI uses 0-based coordinates, so subtract 1 from position
        spdi_position = int(position) - 1
        
        # Query NCBI Variation API using SPDI format
        rsid_url = f"https://api.ncbi.nlm.nih.gov/variation/v0/spdi/{chrom_accession}:{spdi_position}:{ref}:{alt}/rsids"
        
        rsid_response = session.get(rsid_url, timeout=30)
        rsid_response.raise_for_status()
        rsid_json = rsid_response.json()

        if not rsid_json.get("data", {}).get("rsids"):
            return None

        rsid = rsid_json["data"]["rsids"][0]
        return rsid

    except requests.exceptions.RequestException as e:
        return None


@retry_with_backoff(retries=5, backoff_factor=0.5)
def get_variants(lead_url, session=None):
    """
    Fetch all variants from the blood group database.

    :param lead_url: Base URL of the API
    :param session: Requests session object
    :return: Pandas DataFrame containing all variants
    """
    url = f"{lead_url}/variant"
    logger.info(f"Fetching variants from: {url}")

    try:
        response = session.get(url, timeout=30)
        response.raise_for_status()
        variants = pd.DataFrame(response.json())
        logger.info(f"Successfully fetched {len(variants)} variants")
        return variants
    except requests.exceptions.RequestException as e:
        logger.error(f"Error fetching variants: {str(e)}")
        if response is not None:
            logger.error(f"Response content: {response.text}")
        raise


def annotate_rsids(variants, session, lead_url, test_mode=True, overwrite_all=False, clear_not_found=False):
    """
    Annotate variants with rsIDs from dbSNP using GRCh38 coordinates.

    :param variants: Pandas DataFrame containing variants to annotate
    :param session: Authenticated requests session
    :param lead_url: Base URL of the API
    :param test_mode: If True, logs what would be updated without making PATCH requests
    :param overwrite_all: If True, update all variants; if False, only update variants without rsIDs
    :param clear_not_found: If True, clear existing rsID when not found in dbSNP (only with overwrite_all)
    :return: Updated variants DataFrame
    """
    updated_count = 0
    skipped_count = 0
    not_found_count = 0
    cleared_count = 0
    
    for idx, row in variants.iterrows():
        # Get variant ID and coordinates
        db_variant_id = row.get('id')
        chrom = row.get('grch38_chr')
        pos = row.get('grch38_pos')
        ref = row.get('grch38_ref')
        alt = row.get('grch38_alt')
        current_rsid = row.get('rsid')
        
        # Progress indicator
        logger.info(f"[{idx + 1}/{len(variants)}] Processing variant {db_variant_id}: {chrom}:{pos} {ref}>{alt}")
        
        # Skip if missing required coordinates
        if pd.isna(chrom) or pd.isna(pos) or pd.isna(ref) or pd.isna(alt):
            logger.info(f"  ⊘ Skipping - missing GRCh38 coordinates")
            skipped_count += 1
            continue
        
        # Skip if already has rsID (unless overwrite_all is True)
        if not overwrite_all and not pd.isna(current_rsid) and current_rsid:
            logger.info(f"  ⊘ Skipping - already has rsID: {current_rsid}")
            skipped_count += 1
            continue
        
        # Skip if alleles are too long (likely indels that won't work well with the API)
        if len(str(ref)) > 50 or len(str(alt)) > 50:
            logger.info(f"  ⊘ Skipping - alleles too long (ref: {len(str(ref))}, alt: {len(str(alt))})")
            skipped_count += 1
            continue
        
        # Skip if ref and alt are the same (invalid variant)
        if str(ref) == str(alt):
            logger.info(f"  ⊘ Skipping - ref and alt are identical ({ref})")
            skipped_count += 1
            continue
        
        # Normalize chromosome format (remove 'chr' prefix if present)
        chrom_normalized = str(chrom).replace('chr', '').replace('Chr', '').replace('CHR', '')
        
        # Query dbSNP for rsID
        try:
            rsid = get_rsid(chrom_normalized, pos, ref, alt, session=session)
            
            if rsid is None:
                logger.info(f"  → Not found in dbSNP")
                not_found_count += 1
                
                # Clear existing rsID if clear_not_found flag is set and overwrite_all is enabled
                if clear_not_found and overwrite_all and current_rsid:
                    logger.info(f"  → Clearing existing rsID: {current_rsid}")
                    if not test_mode:
                        update_url = f"{lead_url}/variant/{db_variant_id}"
                        update_data = {"rsid": None}
                        update_response = session.patch(update_url, json=update_data, timeout=10)
                        update_response.raise_for_status()
                        logger.info(f"  ✓ Cleared rsID in database")
                        cleared_count += 1
                    else:
                        logger.info(f"  ✓ TEST MODE: Would clear rsID")
                        cleared_count += 1
                continue
            
            # Format rsID with 'rs' prefix
            rsid_formatted = f"rs{rsid}"
            logger.info(f"  ✓ Found rsID: {rsid_formatted}")
            
            # Update database with PATCH request
            if not test_mode:
                update_url = f"{lead_url}/variant/{db_variant_id}"
                update_data = {"rsid": rsid_formatted}
                update_response = session.patch(update_url, json=update_data, timeout=10)
                update_response.raise_for_status()
                logger.info(f"  ✓ Updated in database")
                updated_count += 1
            else:
                logger.info(f"  ✓ TEST MODE: Would update with rsID: {rsid_formatted}")
                updated_count += 1
                
        except requests.exceptions.Timeout:
            logger.error(f"  ✗ Timeout (10s)")
        except requests.exceptions.RequestException as e:
            logger.error(f"  ✗ API Error: {str(e)[:100]}")
        
        # Rate limiting: NCBI allows 3 requests/second without API key, so wait 0.5s between requests
        time.sleep(0.5)
    
    logger.info("=" * 80)
    summary = f"Summary: {updated_count} variants updated, {skipped_count} skipped, {not_found_count} not found in dbSNP"
    if clear_not_found and cleared_count > 0:
        summary += f", {cleared_count} rsIDs cleared"
    logger.info(summary)
    logger.info("=" * 80)
    
    return variants


if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Annotate variants with rsIDs from dbSNP using GRCh38 coordinates',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Test mode - shows what would be updated without making changes
  python annotate_rsid.py --test-mode
  
  # Production mode - update only variants without rsIDs (fill in gaps)
  python annotate_rsid.py
  
  # Update all variants, even those with existing rsIDs
  python annotate_rsid.py --overwrite-all
  
  # Clear rsIDs that can't be found in dbSNP (requires --overwrite-all)
  python annotate_rsid.py --overwrite-all --clear-not-found
  
  # Test with first 10 variants
  python annotate_rsid.py --test-mode --limit 10
  
  # Use custom config file
  python annotate_rsid.py --config my_config.json
  
  # Specify credentials directly
  python annotate_rsid.py --url https://api.blooddatabase.org --email user@example.com --password mypass
        """
    )
    parser.add_argument(
        '--test-mode',
        action='store_true',
        help='Run in test mode (no database updates, only show what would be changed)'
    )
    parser.add_argument(
        '--overwrite-all',
        action='store_true',
        help='Update all variants, even those with existing rsIDs (default: only update variants without rsIDs)'
    )
    parser.add_argument(
        '--clear-not-found',
        action='store_true',
        help='Clear existing rsIDs when not found in dbSNP (only works with --overwrite-all)'
    )
    parser.add_argument(
        '--limit',
        type=int,
        help='Limit number of variants to process (useful for testing)'
    )
    parser.add_argument(
        '--config',
        default='config.json',
        help='Path to configuration JSON file (default: config.json)'
    )
    parser.add_argument(
        '--url',
        help='Base URL of the API (overrides config file)'
    )
    parser.add_argument(
        '--email',
        help='Email for authentication (overrides config file)'
    )
    parser.add_argument(
        '--password',
        help='Password for authentication (overrides config file)'
    )
    
    args = parser.parse_args()
    
    # Load configuration
    if args.url and args.email and args.password:
        # Use command-line arguments
        lead_url = args.url
        email = args.email
        password = args.password
        logger.info("Using credentials from command-line arguments")
    else:
        # Load from config file
        try:
            with open(args.config) as f:
                config = json.load(f)
            
            lead_url = args.url or config['lead_url']
            email = args.email or config['email']
            password = args.password or config['password']
            logger.info(f"Loaded configuration from {args.config}")
        except FileNotFoundError:
            logger.error(f"Config file '{args.config}' not found.")
            logger.error("Either provide --config with a valid file, or use --url, --email, and --password")
            raise
        except KeyError as e:
            logger.error(f"Missing required configuration key: {e}")
            raise

    # Validate flag combination
    if args.clear_not_found and not args.overwrite_all:
        logger.error("--clear-not-found requires --overwrite-all to be set")
        logger.error("Use: python annotate_rsid.py --overwrite-all --clear-not-found")
        raise ValueError("--clear-not-found requires --overwrite-all")

    try:
        # Authenticate and fetch variants
        logger.info("Starting rsID annotation process")
        logger.info(f"Mode: {'TEST MODE (no updates)' if args.test_mode else 'PRODUCTION MODE (will update database)'}")
        logger.info(f"Overwrite existing: {'YES' if args.overwrite_all else 'NO (only update variants without rsIDs)'}")
        if args.clear_not_found:
            logger.info(f"Clear not found: YES (will clear rsIDs not found in dbSNP)")
        
        logger.info(f"Connecting to: {lead_url}")
        session = login(lead_url, email, password)
        
        logger.info("Fetching variants from database...")
        variants = get_variants(lead_url, session=session)
        
        logger.info(f"Fetched {len(variants)} variants from database")
        
        # Run rsID annotation
        logger.info("=" * 80)
        logger.info(f"BEGINNING RSID ANNOTATION")
        logger.info("=" * 80)
        
        # Limit variants if specified
        variants_to_process = variants.head(args.limit) if args.limit else variants
        if args.limit:
            logger.info(f"Processing only first {args.limit} variants")
        
        annotate_rsids(
            variants_to_process,
            session,
            lead_url,
            test_mode=args.test_mode,
            overwrite_all=args.overwrite_all,
            clear_not_found=args.clear_not_found
        )
        
        logger.info("=" * 80)
        logger.info(f"ANNOTATION COMPLETE")
        logger.info("=" * 80)
        
        if args.test_mode:
            logger.info("TEST MODE: No database updates were made")
            logger.info("Run without --test-mode to update the database")
            
    except Exception as e:
        logger.error(f"Fatal error in main process: {str(e)}")
        raise
