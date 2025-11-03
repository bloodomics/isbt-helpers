"""
Script to annotate variants with rsIDs from dbSNP.

This script:
1. Fetches variants from the blood group database
2. Queries dbSNP/Ensembl API for rsIDs
3. Updates the database with PATCH requests

Author: Nick Gleadall
Date: October 2025
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


@retry_with_backoff(retries=5, backoff_factor=0.5)
def get_rsid(chromosome, position, ref, alt, session=None):
    """
    Fetch the rsID for a given genomic variant.

    :param chromosome: Chromosome number
    :param position: Genomic position
    :param ref: Reference allele
    :param alt: Alternate allele
    :param session: Requests session object
    :return: rsID if found, None otherwise
    """
    try:
        # Fetch contextual information
        contextuals_url = f"https://api.ncbi.nlm.nih.gov/variation/v0/vcf/{chromosome}/{position}/{ref}/{alt}/contextuals"
        logger.info(f"Fetching contextuals from: {contextuals_url}")

        response = session.get(contextuals_url, timeout=10)
        response.raise_for_status()
        response_json = response.json()

        if not response_json.get("data", {}).get("spdis"):
            logger.warning(f"No SPDI found for {chromosome}:{position} {ref}>{alt}")
            return None

        # Extract SPDI and fetch rsID
        spdi = response_json["data"]["spdis"][0]
        rsid_url = f"https://api.ncbi.nlm.nih.gov/variation/v0/spdi/{spdi['seq_id']}:{spdi['position']}:{spdi['deleted_sequence']}:{spdi['inserted_sequence']}/rsids"
        logger.info(f"Fetching rsID from: {rsid_url}")

        rsid_response = session.get(rsid_url, timeout=10)
        rsid_response.raise_for_status()
        rsid_json = rsid_response.json()

        if not rsid_json.get("data", {}).get("rsids"):
            logger.warning(
                f"No rsID found for SPDI: {spdi['seq_id']}:{spdi['position']}:{spdi['deleted_sequence']}:{spdi['inserted_sequence']}"
            )
            return None

        rsid = rsid_json["data"]["rsids"][0]
        logger.info(f"Found rsID: {rsid}")
        return rsid

    except requests.exceptions.RequestException as e:
        logger.error(f"Error fetching rsID: {str(e)}")
        return None


@retry_with_backoff(retries=5, backoff_factor=0.5)
def get_variants(session=None):
    """
    Fetch all variants from the blood group database.

    :param session: Requests session object
    :return: Pandas DataFrame containing all variants
    """
    url = "https://api.blooddatabase.org/variant"
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


def parse_chromosomal(chromosomal):
    """
    Parse the chromosomal representation of a variant.

    :param chromosomal: String representation of the variant (e.g., "NC_000001.11:g.12345A>G")
    :return: Tuple of (chromosome, position, ref_allele, alt_allele)
    """
    try:
        parts = chromosomal.split(":g.")
        chromosome = parts[0]
        position_and_alleles = parts[1].split(">")
        position = "".join(char for char in position_and_alleles[0] if char.isdigit())
        ref_allele = position_and_alleles[0][-1]
        alt_allele = position_and_alleles[1]
        return chromosome, position, ref_allele, alt_allele
    except Exception as e:
        logger.error(f"Error parsing chromosomal {chromosomal}: {str(e)}")
        return None, None, None, None


@retry_with_backoff(retries=5, backoff_factor=0.5)
def update_variant(row_id, rsid, lead_url, session=None):
    """
    Update a variant in the blood group database with its rsID.

    :param row_id: ID of the variant to update
    :param rsid: rsID to be added to the variant
    :param lead_url: Base URL of the API
    :param session: Requests session object
    """
    update_url = f"{lead_url}/variant/{row_id}"
    update_data = {"rsid": f"rs{rsid}" if rsid is not None else None}

    try:
        update_response = session.patch(update_url, json=update_data, timeout=10)
        update_response.raise_for_status()
        logger.info(f"Successfully updated variant {row_id} with rsID {rsid}")
    except requests.exceptions.RequestException as e:
        logger.error(f"Error updating variant {row_id}: {str(e)}")
        if update_response is not None:
            logger.error(f"Response content: {update_response.text}")
        raise


def annotate_rsid(row, get_rsid_func, lead_url, session):
    """
    Annotate a single variant with its rsID.

    :param row: Pandas Series representing a single variant
    :param get_rsid_func: Function to fetch rsID
    :param lead_url: Base URL of the API
    :param session: Requests session object
    """
    chromosomal = row["hgvs_genomic_grch38"]
    if pd.isna(chromosomal):
        return

    chromosome, position, ref_allele, alt_allele = parse_chromosomal(chromosomal)
    if chromosome is None:
        return

    if len(ref_allele) > 50 or len(alt_allele) > 50:
        logger.info(f"Skipping variant {row['id']} due to long alleles")
        return

    try:
        rsid = get_rsid_func(
            chromosome, position, ref_allele, alt_allele, session=session
        )

        if rsid == row["rsid"]:
            logger.info(f"rsID {rsid} for variant {row['id']} is unchanged")
            return

        update_variant(row["id"], rsid, lead_url, session=session)
    except Exception as e:
        logger.error(f"Error processing variant {row['id']}: {str(e)}")


if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Annotate variants with rsIDs from dbSNP',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use config.json file (default)
  python annotate_rsid.py
  
  # Specify config file
  python annotate_rsid.py --config my_config.json
  
  # Specify credentials directly
  python annotate_rsid.py --url https://api.blooddatabase.org --email user@example.com --password mypass
        """
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
        logger.info("Starting rsID annotation process")
        
        # Login and get authenticated session
        session = login(lead_url, email, password)

        # Fetch all variants
        variants = get_variants(session=session)
        logger.info(f"Fetched {len(variants)} variants from database")

        # Process each variant
        logger.info("=" * 80)
        logger.info("BEGINNING RSID ANNOTATION")
        logger.info("=" * 80)
        
        for _, row in variants.iterrows():
            annotate_rsid(row, get_rsid, lead_url, session)
        
        logger.info("=" * 80)
        logger.info("ANNOTATION COMPLETE")
        logger.info("=" * 80)
        
    except Exception as e:
        logger.error(f"Fatal error in main process: {str(e)}")
        raise
