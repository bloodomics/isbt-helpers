# ISBT Helpers

Helper scripts for managing and annotating blood group variant data in the ISBT database.

## Scripts

### `annotate_gnomad.py`

Annotates variants with population frequency data from gnomAD v4 using GRCh38 coordinates.

**Features:**

- Uses gnomAD v4 GraphQL API (GRCh38/hg38 reference genome)
- Calculates Minor Allele Frequency (MAF) for all populations
- Test mode for dry-run without database updates
- Optional overwrite mode to update all variants (default: only update variants without gnomAD data)
- Configurable processing limit for testing
- Rate-limited to 10 requests per 60 seconds (gnomAD API limit)

**Populations annotated:**

- `gnomad_all` - All populations combined
- `gnomad_afr` - African/African American
- `gnomad_amr` - Admixed American
- `gnomad_asj` - Ashkenazi Jewish
- `gnomad_eas` - East Asian
- `gnomad_fin` - Finnish
- `gnomad_nfe` - Non-Finnish European
- `gnomad_oth` - Other (remaining)
- `gnomad_sas` - South Asian

**Usage:**

```bash
# Test mode - shows what would be updated without making changes
python annotate_gnomad.py --test-mode

# Production mode - update only variants without gnomAD data (fill in gaps)
python annotate_gnomad.py

# Update all variants, even those with existing gnomAD data
python annotate_gnomad.py --overwrite-all

# Test with first 10 variants
python annotate_gnomad.py --test-mode --limit 10

# Use custom config file
python annotate_gnomad.py --config my_config.json
```

### `annotate_rsid.py`

Annotates variants with rsIDs from dbSNP using GRCh38 coordinates.

**Features:**

- Uses NCBI Variation API with SPDI format and GRCh38 RefSeq accessions
- Queries dbSNP by chromosome, position, reference, and alternate alleles
- Test mode for dry-run without database updates
- Optional overwrite mode to update all variants (default: only update variants without rsIDs)
- Optional clear mode to remove rsIDs that can't be found in dbSNP
- Configurable processing limit for testing
- Rate-limited to 2 requests per second (NCBI allows 3 req/sec without API key)
- Automatic chromosome format normalization (strips 'chr' prefix)
- Validates variants (skips if ref == alt or alleles too long)
- Automatic retry with exponential backoff for failed requests
- 30-second timeout with retry logic

**Usage:**

```bash
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
```

### `annotate_exons.py`

Annotates variants with exon and intron numbers using HGVS transcript coordinates via VariantValidator API.

**Features:**

- Uses VariantValidator API with GRCh38 coordinates
- Annotates both exon and intron numbers
- Handles ranges when variants span multiple exons/introns (e.g., "5-6")
- Test mode for dry-run without database updates
- Optional overwrite mode to update all variants (default: only update variants without exon/intron data)
- Configurable processing limit for testing
- Rate-limited to 4 requests per second (VariantValidator recommendation)
- Automatic retry with exponential backoff for failed requests

**Usage:**

```bash
# Test mode - shows what would be updated without making changes
python annotate_exons.py --test-mode

# Production mode - update only variants without exon/intron data (fill in gaps)
python annotate_exons.py

# Update all variants, even those with existing exon/intron data
python annotate_exons.py --overwrite-all

# Test with first 10 variants
python annotate_exons.py --test-mode --limit 10

# Use custom config file
python annotate_exons.py --config my_config.json
```

### `export_for_isbt.py`

Exports blood group allele data in Excel format for ISBT submission.

**Features:**

- Exports allele tables for blood group systems
- Multi-threaded processing for faster exports
- Can export specific systems or all systems
- Generates formatted Excel files with proper column alignment

**Usage:**

```bash
# Export all systems
python export_for_isbt.py --lead-url https://api.example.com --output ./exports

# Export specific system
python export_for_isbt.py --lead-url https://api.example.com --output ./exports --system ABO

# Use multiple threads for faster processing
python export_for_isbt.py --lead-url https://api.example.com --output ./exports --threads 4
```

## Setup

1. Install dependencies:

```bash
pip install -r requirements.txt
```

2. Create a `config.json` file (copy from `config.json.example`):

```json
{
  "lead_url": "https://api.example.com",
  "email": "your-email@example.com",
  "password": "your-password"
}
```

**Important:** The `config.json` file is gitignored and should never be committed to the repository.

3. Alternatively, you can pass credentials via command-line arguments:

```bash
python annotate_gnomad.py --url https://api.example.com --email user@example.com --password yourpass
```

## Requirements

- Python 3.8+
- pandas >= 2.0.0
- requests >= 2.31.0
- openpyxl (for export_for_isbt.py)

See `requirements.txt` for full dependency list.

## API Details

### gnomAD v4 GraphQL API

- Endpoint: `https://gnomad.broadinstitute.org/api`
- Dataset: `gnomad_r4`
- Reference: GRCh38/hg38
- Rate Limit: 10 requests per 60 seconds
- Authentication: None required

### NCBI Variation API

- Endpoint: `https://api.ncbi.nlm.nih.gov/variation/v0/`
- Format: SPDI (Sequence Position Deletion Insertion)
- Reference: GRCh38 RefSeq accessions (NC_000001.11 through NC_000024.10)
- Coordinate System: 0-based (script handles conversion from 1-based)
- Rate Limit: 3 requests per second (without API key)
- Authentication: None required

### VariantValidator API

- Endpoint: `https://rest.variantvalidator.org/`
- Reference: GRCh38
- Rate Limit: 4 requests per second (recommended)
- Authentication: None required

## Data Validation

The annotation scripts include automatic validation:

### annotate_rsid.py validations:

- Skips variants with missing GRCh38 coordinates
- Skips variants where ref and alt alleles are identical
- Skips variants with alleles longer than 50 characters
- Automatically normalizes chromosome format (removes 'chr' prefix)

### annotate_gnomad.py validations:

- Skips variants with missing GRCh38 coordinates
- Converts allele frequency (AF) to Minor Allele Frequency (MAF): if AF > 0.5, MAF = 1 - AF

### annotate_exons.py validations:

- Skips variants without HGVS transcript annotations
- Handles variants spanning multiple exons/introns

## Troubleshooting

### Rate Limiting and Timeouts

All annotation scripts implement rate limiting and automatic retry logic:

- **gnomAD**: 6.5-second delay between requests (10 req/60s limit)
- **NCBI**: 0.5-second delay between requests (2 req/s, under 3 req/s limit)
- **VariantValidator**: 0.25-second delay between requests (4 req/s)

If you encounter timeout errors:

- Scripts automatically retry with exponential backoff (up to 5 retries)
- Timeout is set to 30 seconds for NCBI queries
- Rate limiting applies even on failed requests to avoid hammering APIs

### Common Issues

**"Not found in dbSNP"**: This can occur if:

- Variant is very rare and not in dbSNP
- Coordinates don't exactly match dbSNP entry
- Variant is invalid (ref == alt)

**"Timeout (30s)"**: Occasional timeouts from NCBI are normal:

- Script will automatically retry
- If persistent, may indicate API service issues

**Invalid variants skipped**: The scripts automatically skip:

- Variants with identical ref and alt alleles
- Variants with missing coordinates
- Variants with very long alleles (>50 characters for rsid annotation)

### Clearing Invalid Data

Use the `--clear-not-found` flag with `--overwrite-all` to remove rsIDs that can't be found:

```bash
python annotate_rsid.py --overwrite-all --clear-not-found
```

This is useful for cleaning up outdated or incorrect rsIDs in the database.

## Development

**Security:**

- Never commit credentials or API keys to the repository
- Use `config.json` for local development (gitignored)
- Use environment variables or secrets management for production

**Testing:**

- Always use `--test-mode` flag when testing scripts to avoid database updates
- Use `--limit` flag to process a small number of variants for testing
- Test with `--limit 10` to verify API connectivity and authentication before full runs

**Best Practices:**

- Run scripts during off-peak hours for large datasets
- Use `--overwrite-all` sparingly - it processes all variants even those with existing data
- Monitor logs for skipped variants and validation messages
- Keep track of processed variants in case runs are interrupted
