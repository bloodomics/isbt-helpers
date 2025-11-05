"""
Script to annotate variants with exon and intron numbers using HGVS transcript coordinates.

This script:
1. Fetches variants from the blood group database
2. Queries VariantValidator API for exon/intron positions using HGVS transcript
3. Updates the database with PATCH requests (only for variants without existing exon/intron data)

Rate Limit: VariantValidator recommends max 4 requests per second

Author: Nick Gleadall
Date: November 2025
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
        login_response = session.post(login_url, json=login_data, timeout=30)
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
    :param session: Authenticated requests session
    :return: Pandas DataFrame containing all variants
    """
    url = f'{lead_url}/variant'
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


def annotate_exons_introns(variants, session, lead_url, test_mode=True, overwrite_all=False, clear_not_found=False):
    """
    Annotate variants with exon and intron numbers from VariantValidator.

    :param variants: Pandas DataFrame containing variants to annotate
    :param session: Authenticated requests session
    :param lead_url: Base URL of the API
    :param test_mode: If True, logs what would be updated without making PATCH requests
    :param overwrite_all: If True, update all variants; if False, only update variants without exon/intron data
    :param clear_not_found: If True, clear existing exon/intron when not found in VariantValidator (only with overwrite_all)
    :return: Updated variants DataFrame
    """
    updated_count = 0
    skipped_count = 0
    not_found_count = 0
    cleared_count = 0
    
    for idx, row in variants.iterrows():
        # Get variant ID and HGVS transcript
        db_variant_id = row.get('id')
        hgvs_transcript = row.get('hgvs_transcript')
        current_exon = row.get('exon')
        current_intron = row.get('intron')
        
        # Progress indicator
        logger.info(f"[{idx + 1}/{len(variants)}] Processing variant {db_variant_id}")
        
        # Skip if no HGVS transcript
        if pd.isna(hgvs_transcript) or not hgvs_transcript:
            logger.info(f"  ⊘ Skipping - no HGVS transcript")
            skipped_count += 1
            continue
        
        # Skip if already has exon/intron data (unless overwrite_all is True)
        if not overwrite_all and ((not pd.isna(current_exon) and current_exon) or (not pd.isna(current_intron) and current_intron)):
            logger.info(f"  ⊘ Skipping - already has exon/intron: {current_exon}/{current_intron}")
            skipped_count += 1
            continue
        
        # Extract transcript ID from HGVS
        try:
            transcript = hgvs_transcript.split(":")[0]
        except Exception as e:
            logger.error(f"  ✗ Could not parse HGVS transcript: {hgvs_transcript}")
            not_found_count += 1
            continue
        
        # Query VariantValidator API
        vv_url = f"https://rest.variantvalidator.org/VariantValidator/variantvalidator/GRCh38/{hgvs_transcript}/{transcript}"
        
        try:
            response = session.get(vv_url, timeout=30)
            
            if response.status_code != 200:
                logger.info(f"  → Not found in VariantValidator (status {response.status_code})")
                not_found_count += 1
                
                # Clear existing exon/intron if clear_not_found flag is set
                if clear_not_found and overwrite_all and (current_exon or current_intron):
                    logger.info(f"  → Clearing existing exon/intron: {current_exon}/{current_intron}")
                    if not test_mode:
                        update_url = f"{lead_url}/variant/{db_variant_id}"
                        clear_data = {"exon": None, "intron": None}
                        update_response = session.patch(update_url, json=clear_data, timeout=10)
                        update_response.raise_for_status()
                        logger.info(f"  ✓ Cleared exon/intron in database")
                        cleared_count += 1
                    else:
                        logger.info(f"  ✓ TEST MODE: Would clear exon/intron")
                        cleared_count += 1
                continue
            
            response_data = response.json()
            
            if hgvs_transcript not in response_data:
                logger.info(f"  → HGVS transcript not in response")
                not_found_count += 1
                
                # Clear existing exon/intron if clear_not_found flag is set
                if clear_not_found and overwrite_all and (current_exon or current_intron):
                    logger.info(f"  → Clearing existing exon/intron: {current_exon}/{current_intron}")
                    if not test_mode:
                        update_url = f"{lead_url}/variant/{db_variant_id}"
                        clear_data = {"exon": None, "intron": None}
                        update_response = session.patch(update_url, json=clear_data, timeout=10)
                        update_response.raise_for_status()
                        logger.info(f"  ✓ Cleared exon/intron in database")
                        cleared_count += 1
                    else:
                        logger.info(f"  ✓ TEST MODE: Would clear exon/intron")
                        cleared_count += 1
                continue
            
            variant_data = response_data[hgvs_transcript]
            variant_exonic_positions = variant_data.get("variant_exonic_positions", {})
            
            if not variant_exonic_positions:
                logger.info(f"  → No exonic positions found")
                not_found_count += 1
                
                # Clear existing exon/intron if clear_not_found flag is set
                if clear_not_found and overwrite_all and (current_exon or current_intron):
                    logger.info(f"  → Clearing existing exon/intron: {current_exon}/{current_intron}")
                    if not test_mode:
                        update_url = f"{lead_url}/variant/{db_variant_id}"
                        clear_data = {"exon": None, "intron": None}
                        update_response = session.patch(update_url, json=clear_data, timeout=10)
                        update_response.raise_for_status()
                        logger.info(f"  ✓ Cleared exon/intron in database")
                        cleared_count += 1
                    else:
                        logger.info(f"  ✓ TEST MODE: Would clear exon/intron")
                        cleared_count += 1
                continue
            
            # Get the latest GRCh38 version (highest NC_*.{version})
            start_exon = None
            end_exon = None
            start_intron = None
            end_intron = None
            latest_version = 0
            
            for key in variant_exonic_positions.keys():
                if key.startswith("NC"):
                    try:
                        version = int(key.split(".")[1])
                        if version > latest_version:
                            latest_version = version
                            pos_data = variant_exonic_positions[key]
                            start_exon = pos_data.get("start_exon")
                            end_exon = pos_data.get("end_exon")
                            start_intron = pos_data.get("start_intron")
                            end_intron = pos_data.get("end_intron")
                    except (IndexError, ValueError):
                        continue
            
            # Format exon and intron values
            new_exon_value = None
            new_intron_value = None
            
            if start_exon and end_exon:
                if start_exon == end_exon:
                    new_exon_value = str(start_exon)
                else:
                    new_exon_value = f"{start_exon}-{end_exon}"
            
            if start_intron and end_intron:
                if start_intron == end_intron:
                    new_intron_value = str(start_intron)
                else:
                    new_intron_value = f"{start_intron}-{end_intron}"
            
            # Check if we found anything
            if not new_exon_value and not new_intron_value:
                logger.info(f"  → No exon or intron data found")
                not_found_count += 1
                
                # Clear existing exon/intron if clear_not_found flag is set
                if clear_not_found and overwrite_all and (current_exon or current_intron):
                    logger.info(f"  → Clearing existing exon/intron: {current_exon}/{current_intron}")
                    if not test_mode:
                        update_url = f"{lead_url}/variant/{db_variant_id}"
                        clear_data = {"exon": None, "intron": None}
                        update_response = session.patch(update_url, json=clear_data, timeout=10)
                        update_response.raise_for_status()
                        logger.info(f"  ✓ Cleared exon/intron in database")
                        cleared_count += 1
                    else:
                        logger.info(f"  ✓ TEST MODE: Would clear exon/intron")
                        cleared_count += 1
                continue
            
            logger.info(f"  ✓ Found - Exon: {new_exon_value}, Intron: {new_intron_value}")
            
            # Update database with PATCH request
            if not test_mode:
                update_url = f"{lead_url}/variant/{db_variant_id}"
                update_data = {}
                if new_exon_value:
                    update_data["exon"] = new_exon_value
                if new_intron_value:
                    update_data["intron"] = new_intron_value
                
                update_response = session.patch(update_url, json=update_data, timeout=10)
                update_response.raise_for_status()
                logger.info(f"  ✓ Updated in database")
                updated_count += 1
            else:
                logger.info(f"  ✓ TEST MODE: Would update with exon={new_exon_value}, intron={new_intron_value}")
                updated_count += 1
                
        except requests.exceptions.Timeout:
            logger.error(f"  ✗ Timeout (30s)")
        except requests.exceptions.RequestException as e:
            logger.error(f"  ✗ API Error: {str(e)[:100]}")
        except Exception as e:
            logger.error(f"  ✗ Unexpected error: {str(e)[:100]}")
        
        # Rate limiting: VariantValidator recommends max 4 requests per second
        time.sleep(0.25)
    
    logger.info("=" * 80)
    summary = f"Summary: {updated_count} variants updated, {skipped_count} skipped, {not_found_count} not found"
    if clear_not_found and cleared_count > 0:
        summary += f", {cleared_count} exon/intron records cleared"
    logger.info(summary)
    logger.info("=" * 80)
    
    return variants


if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Annotate variants with exon and intron numbers using HGVS transcript',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Test mode - shows what would be updated without making changes
  python annotate_exons.py --test-mode
  
  # Production mode - update only variants without exon/intron data (fill in gaps)
  python annotate_exons.py
  
  # Update all variants, even those with existing exon/intron data
  python annotate_exons.py --overwrite-all
  
  # Clear exon/intron that can't be found in VariantValidator (requires --overwrite-all)
  python annotate_exons.py --overwrite-all --clear-not-found
  
  # Test with first 10 variants
  python annotate_exons.py --test-mode --limit 10
  
  # Use custom config file
  python annotate_exons.py --config my_config.json
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
        help='Update all variants, even those with existing exon/intron data (default: only update variants without data)'
    )
    parser.add_argument(
        '--clear-not-found',
        action='store_true',
        help='Clear existing exon/intron when not found in VariantValidator (only works with --overwrite-all)'
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
        logger.error("Use: python annotate_exons.py --overwrite-all --clear-not-found")
        raise ValueError("--clear-not-found requires --overwrite-all")

    try:
        # Authenticate and fetch variants
        logger.info("Starting exon/intron annotation process")
        logger.info(f"Mode: {'TEST MODE (no updates)' if args.test_mode else 'PRODUCTION MODE (will update database)'}")
        logger.info(f"Overwrite existing: {'YES' if args.overwrite_all else 'NO (only update variants without exon/intron data)'}")
        if args.clear_not_found:
            logger.info(f"Clear not found: YES (will clear exon/intron not found in VariantValidator)")
        
        logger.info(f"Connecting to: {lead_url}")
        session = login(lead_url, email, password)
        
        logger.info("Fetching variants from database...")
        variants = get_variants(lead_url, session)
        
        logger.info(f"Fetched {len(variants)} variants from database")
        
        # Run exon/intron annotation
        logger.info("=" * 80)
        logger.info(f"BEGINNING EXON/INTRON ANNOTATION")
        logger.info("=" * 80)
        
        # Limit variants if specified
        variants_to_process = variants.head(args.limit) if args.limit else variants
        if args.limit:
            logger.info(f"Processing only first {args.limit} variants")
        
        annotate_exons_introns(
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
