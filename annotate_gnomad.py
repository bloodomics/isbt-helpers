"""
Script to annotate variants with gnomAD v4 population frequencies.

This script:
1. Fetches variants from the blood group database
2. Queries gnomAD v4 GraphQL API using GRCh38 coordinates
3. Calculates Minor Allele Frequency (MAF) for each population
4. Updates the database with PATCH requests (only for variants without existing gnomAD data)

Author: Nick Gleadall
Date: October 2025
"""

import requests
import json
import pandas as pd
import time
import logging
import argparse

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)



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

def get_variants(lead_url, session):
    """
    Fetch all variants from the blood group database.

    :param lead_url: Base URL of the API
    :param session: Requests session object
    :return: Pandas DataFrame containing all variants
    """
    url = f'{lead_url}/variant'
    logger.info(f"Fetching variants from: {url}")

    try:
        response = session.get(url)
        response.raise_for_status()
        variants = pd.DataFrame(response.json())
        logger.info(f"Successfully fetched {len(variants)} variants")
        return variants
    except requests.exceptions.RequestException as e:
        logger.error(f"Error fetching variants: {str(e)}")
        if response is not None:
            logger.error(f"Response content: {response.text}")
        raise

def calculate_maf(af):
    """
    Calculate Minor Allele Frequency (MAF) from Allele Frequency (AF).
    
    MAF is the frequency of the less common allele in a population.
    If AF > 0.5, we use 1 - AF to get the minor allele frequency.
    
    :param af: Allele frequency (0 to 1)
    :return: Minor allele frequency (0 to 0.5) or None
    """
    if af is None:
        return None
    return af if af <= 0.5 else (1 - af)


def annotate_gnomad_frequencies(variants, session, lead_url, test_mode=True, overwrite_all=False, clear_not_found=False):
    """
    Fetch gnomAD v4 frequencies for variants using GRCh38 coordinates and update the database.
    
    Only processes variants that:
    - Have valid GRCh38 coordinates
    - Are found in gnomAD v4
    - Don't have gnomAD data (unless overwrite_all=True)
    
    Updates are made via PATCH requests containing only gnomAD frequency fields.

    :param variants: Pandas DataFrame containing variants with GRCh38 coordinates
    :param session: Authenticated session for database updates
    :param lead_url: Base URL of the API
    :param test_mode: If True, logs what would be updated without making PATCH requests
    :param overwrite_all: If True, updates all variants even if they have existing gnomAD data
    :param clear_not_found: If True, clear existing gnomAD data when not found in gnomAD (only with overwrite_all)
    :return: Pandas DataFrame (unchanged)
    """
    # gnomAD v4 superpopulation codes
    gnomad_populations = [
        'gnomad_all',   # All populations combined
        'gnomad_afr',   # African/African American
        'gnomad_amr',   # Admixed American
        'gnomad_asj',   # Ashkenazi Jewish
        'gnomad_eas',   # East Asian
        'gnomad_fin',   # Finnish
        'gnomad_nfe',   # Non-Finnish European
        'gnomad_oth',   # Other (population not assigned)
        'gnomad_sas'    # South Asian
    ]
    
    updated_count = 0
    skipped_count = 0
    not_found_count = 0
    cleared_count = 0
    processed_count = 0

    total_variants = len(variants)
    
    for index, row in variants.iterrows():
        processed_count += 1
        db_variant_id = row.get('id')
        
        logger.info(f"[{processed_count}/{total_variants}] Processing variant {db_variant_id}")
        
        # Check if variant already has gnomAD frequencies (unless overwrite_all is True)
        if not overwrite_all and pd.notna(row.get('gnomad_all')):
            logger.info(f"  → Already has gnomAD frequencies, skipping")
            skipped_count += 1
            continue
        
        # Use GRCh38/hg38 coordinates
        chrom = row.get('grch38_chr')
        pos = row.get('grch38_pos')
        ref = row.get('grch38_ref')
        alt = row.get('grch38_alt')
        
        # Skip if any required field is missing
        if pd.isna(chrom) or pd.isna(pos) or pd.isna(ref) or pd.isna(alt):
            logger.warning(f"  → Missing GRCh38 coordinates, skipping")
            skipped_count += 1
            continue
        
        # Skip if REF or ALT alleles are too long (gnomAD API has size limits)
        # Large variants (>1000 bp) often cause 413 Request Entity Too Large errors
        max_allele_length = 1000
        if len(str(ref)) > max_allele_length or len(str(alt)) > max_allele_length:
            logger.warning(f"  → REF or ALT allele too long (REF:{len(str(ref))} bp, ALT:{len(str(alt))} bp), skipping")
            skipped_count += 1
            continue

        # gnomAD v4 GraphQL API
        gnomad_api_url = "https://gnomad.broadinstitute.org/api"
        gnomad_variant_id = f"{chrom}-{pos}-{ref}-{alt}"  # gnomAD variant ID format
        
        # GraphQL query for gnomAD v4
        # Request only main superpopulations to keep response size manageable
        query = """
        query GnomadVariant($variantId: String!, $datasetId: DatasetId!) {
          variant(variantId: $variantId, dataset: $datasetId) {
            variant_id
            genome {
              af
              populations {
                id
                ac
                an
              }
            }
          }
        }
        """
        
        variables = {
            "variantId": gnomad_variant_id,
            "datasetId": "gnomad_r4"
        }
        
        logger.info(f"  → Querying gnomAD API: {gnomad_variant_id}")

        try:
            response = requests.post(
                gnomad_api_url,
                json={"query": query, "variables": variables},
                headers={"Content-Type": "application/json"},
                timeout=30  # Add 30 second timeout
            )
            response.raise_for_status()
            data = response.json()
            
            # Check for GraphQL errors
            if 'errors' in data:
                logger.info(f"  → Not found in gnomAD (GraphQL error)")
                not_found_count += 1
                
                # Clear existing gnomAD data if clear_not_found flag is set
                if clear_not_found and overwrite_all and pd.notna(row.get('gnomad_all')):
                    logger.info(f"  → Clearing existing gnomAD frequencies")
                    if not test_mode:
                        update_url = f"{lead_url}/variant/{db_variant_id}"
                        clear_data = {pop: None for pop in gnomad_populations}
                        update_response = session.patch(update_url, json=clear_data)
                        update_response.raise_for_status()
                        logger.info(f"  ✓ Cleared gnomAD data in database")
                        cleared_count += 1
                    else:
                        logger.info(f"  ✓ TEST MODE: Would clear gnomAD data")
                        cleared_count += 1
                continue
            
            # Check if variant exists in gnomAD
            variant_data = data.get('data', {}).get('variant')
            if not variant_data:
                logger.info(f"  → Not found in gnomAD")
                not_found_count += 1
                
                # Clear existing gnomAD data if clear_not_found flag is set
                if clear_not_found and overwrite_all and pd.notna(row.get('gnomad_all')):
                    logger.info(f"  → Clearing existing gnomAD frequencies")
                    if not test_mode:
                        update_url = f"{lead_url}/variant/{db_variant_id}"
                        clear_data = {pop: None for pop in gnomad_populations}
                        update_response = session.patch(update_url, json=clear_data)
                        update_response.raise_for_status()
                        logger.info(f"  ✓ Cleared gnomAD data in database")
                        cleared_count += 1
                    else:
                        logger.info(f"  ✓ TEST MODE: Would clear gnomAD data")
                        cleared_count += 1
                continue
            
            # Extract gnomAD frequencies from GraphQL response
            genome_data = variant_data.get('genome')
            if not genome_data:
                logger.info(f"  → Found but no genome data available")
                not_found_count += 1
                
                # Clear existing gnomAD data if clear_not_found flag is set
                if clear_not_found and overwrite_all and pd.notna(row.get('gnomad_all')):
                    logger.info(f"  → Clearing existing gnomAD frequencies")
                    if not test_mode:
                        update_url = f"{lead_url}/variant/{db_variant_id}"
                        clear_data = {pop: None for pop in gnomad_populations}
                        update_response = session.patch(update_url, json=clear_data)
                        update_response.raise_for_status()
                        logger.info(f"  ✓ Cleared gnomAD data in database")
                        cleared_count += 1
                    else:
                        logger.info(f"  ✓ TEST MODE: Would clear gnomAD data")
                        cleared_count += 1
                continue
            populations = genome_data.get('populations', [])
            
            # Calculate MAF for each population from AC/AN (allele count / allele number)
            pop_mafs = {}
            
            for pop in populations:
                ac = pop.get('ac')
                an = pop.get('an')
                pop_id = pop['id']
                
                if ac is not None and an is not None and an > 0:
                    af = ac / an
                    pop_mafs[pop_id] = calculate_maf(af)
                else:
                    pop_mafs[pop_id] = None
            
            # Calculate MAF for overall frequency
            overall_af = genome_data.get('af')
            overall_maf = calculate_maf(overall_af)
            
            # Map gnomAD population IDs to database fields
            # Note: gnomAD v4 uses 'remaining' for "Other" population
            # Use explicit None check to preserve 0 values (0 is a valid MAF)
            remaining_maf = pop_mafs.get('remaining')
            gnomad_freqs = {
                'gnomad_all': overall_maf,
                'gnomad_afr': pop_mafs.get('afr'),
                'gnomad_amr': pop_mafs.get('amr'),
                'gnomad_asj': pop_mafs.get('asj'),
                'gnomad_eas': pop_mafs.get('eas'),
                'gnomad_fin': pop_mafs.get('fin'),
                'gnomad_nfe': pop_mafs.get('nfe'),
                'gnomad_oth': remaining_maf if remaining_maf is not None else pop_mafs.get('oth'),  # gnomAD v4 uses 'remaining'
                'gnomad_sas': pop_mafs.get('sas')
            }
            
            logger.info(f"  ✓ Found in gnomAD (MAF: {overall_maf:.6f})")
            
            # Update database with PATCH request (only gnomAD fields)
            if not test_mode:
                update_url = f"{lead_url}/variant/{db_variant_id}"
                update_response = session.patch(update_url, json=gnomad_freqs)
                update_response.raise_for_status()
                logger.info(f"  ✓ Updated in database")
                updated_count += 1
            else:
                logger.info(f"  ✓ TEST MODE: Would update with MAF data")
                updated_count += 1
            
        except requests.exceptions.Timeout:
            logger.error(f"  ✗ Timeout (30s)")
        except requests.exceptions.RequestException as e:
            logger.error(f"  ✗ API Error: {str(e)[:100]}")
        
        # Rate limiting: gnomAD allows 10 requests per 60 seconds, so wait 6.5 seconds between requests
        # This runs for ALL requests (success or failure) since failed requests also count against the limit
        time.sleep(6.5)

    logger.info("=" * 80)
    summary = f"Summary: {updated_count} variants updated, {skipped_count} skipped (already have data), {not_found_count} not found in gnomAD"
    if clear_not_found and cleared_count > 0:
        summary += f", {cleared_count} gnomAD records cleared"
    logger.info(summary)
    logger.info("=" * 80)
    
    return variants

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Annotate variants with gnomAD v4 population frequencies',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Test mode - shows what would be updated without making changes
  python annotate_gnomad.py --test-mode
  
  # Production mode - updates database (use config.json)
  python annotate_gnomad.py
  
  # Production mode with command-line credentials
  python annotate_gnomad.py --url https://api.blooddatabase.org --email user@example.com --password mypass
  
  # Overwrite all existing gnomAD values
  python annotate_gnomad.py --overwrite-all
  
  # Clear gnomAD data that can't be found (requires --overwrite-all)
  python annotate_gnomad.py --overwrite-all --clear-not-found
  
  # Test overwrite mode
  python annotate_gnomad.py --test-mode --overwrite-all
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
        help='Update all variants, even those with existing gnomAD frequencies'
    )
    parser.add_argument(
        '--clear-not-found',
        action='store_true',
        help='Clear existing gnomAD data when not found in gnomAD (only works with --overwrite-all)'
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

    try:
        # Authenticate and fetch variants
        logger.info("Starting gnomAD annotation process")
        logger.info(f"Mode: {'TEST MODE (no updates)' if args.test_mode else 'PRODUCTION MODE (will update database)'}")
        logger.info(f"Overwrite existing: {'YES' if args.overwrite_all else 'NO (only update variants without gnomAD data)'}")
        
        logger.info(f"Connecting to: {lead_url}")
        session = login(lead_url, email, password)
        
        logger.info("Fetching variants from database...")
        variants = get_variants(lead_url, session)
        
        logger.info(f"Fetched {len(variants)} variants from database")
        
        # Run gnomAD annotation
        logger.info("=" * 80)
        logger.info(f"BEGINNING GNOMAD ANNOTATION")
        logger.info("=" * 80)
        
        # Limit variants if specified
        variants_to_process = variants.head(args.limit) if args.limit else variants
        if args.limit:
            logger.info(f"Processing only first {args.limit} variants")
        
        annotate_gnomad_frequencies(
            variants_to_process, 
            session, 
            lead_url, 
            test_mode=args.test_mode,
            overwrite_all=args.overwrite_all
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