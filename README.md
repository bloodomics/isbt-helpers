# ISBT Helpers

Helper scripts for managing and annotating blood group variant data in the ISBT database.

## Scripts

### `annotate_gnomad.py`

Annotates variants with population frequency data from gnomAD v4.

**Features:**

- Uses gnomAD v4 GraphQL API with GRCh38 coordinates
- Calculates Minor Allele Frequency (MAF) for all populations
- Only updates variants without existing gnomAD data
- Skips variants not found in gnomAD
- Rate-limited to respect gnomAD API

**Populations annotated:**

- `gnomad_all` - All populations combined
- `gnomad_afr` - African/African American
- `gnomad_amr` - Admixed American
- `gnomad_asj` - Ashkenazi Jewish
- `gnomad_eas` - East Asian
- `gnomad_fin` - Finnish
- `gnomad_nfe` - Non-Finnish European
- `gnomad_oth` - Other
- `gnomad_sas` - South Asian

**Usage:**

```bash
python annotate_gnomad.py
```

**Requirements:**

- Python 3.8+
- pandas
- requests

### `annotate_exons.py`

_Description to be added_

### `annotate_rsid.py`

_Description to be added_

### `export_for_isbt.py`

_Description to be added_

## Setup

1. Install dependencies:

```bash
pip install pandas requests
```

2. Configure credentials (do not commit to git):
   - Edit scripts to use environment variables or a config file
   - Or update credentials directly in the scripts (not recommended)

## Development

Please ensure credentials are never committed to the repository.
