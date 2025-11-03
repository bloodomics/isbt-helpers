# Pre-Commit Checklist for annotate_gnomad.py

## ⚠️ IMPORTANT: Remove Credentials Before Committing!

Before pushing to git, you MUST remove your credentials from `annotate_gnomad.py`:

### Option 1: Use Environment Variables (Recommended)

Replace the credentials section with:

```python
import os

if __name__ == "__main__":
    # Configuration from environment variables
    lead_url = os.getenv("BLOOD_DB_URL", "https://api.blooddatabase.org")
    email = os.getenv("BLOOD_DB_EMAIL")
    password = os.getenv("BLOOD_DB_PASSWORD")

    if not email or not password:
        raise ValueError("Please set BLOOD_DB_EMAIL and BLOOD_DB_PASSWORD environment variables")
```

Then run with:

```bash
export BLOOD_DB_EMAIL="your_email@example.com"
export BLOOD_DB_PASSWORD="your_password"
python annotate_gnomad.py
```

### Option 2: Use a Config File (Alternative)

Create a `config.json` file (add to .gitignore):

```json
{
  "lead_url": "https://api.blooddatabase.org",
  "email": "your_email@example.com",
  "password": "your_password"
}
```

Then update the script:

```python
import json

if __name__ == "__main__":
    # Load configuration from file
    with open('config.json') as f:
        config = json.load(f)

    lead_url = config['lead_url']
    email = config['email']
    password = config['password']
```

And add to `.gitignore`:

```
config.json
```

## Files Ready to Commit:

- ✅ `annotate_gnomad.py` (after removing credentials!)
- ✅ `README.md` (updated with documentation)
- ✅ `requirements.txt` (new file)
- ✅ `.gitignore` (already exists)

## Git Commands:

```bash
# Check what will be committed
git status

# Review changes
git diff annotate_gnomad.py

# Add files
git add annotate_gnomad.py README.md requirements.txt

# Commit
git commit -m "Add gnomAD v4 annotation script with MAF calculation

- Fetch population frequencies from gnomAD v4 GraphQL API
- Calculate Minor Allele Frequency (MAF) for 9 populations
- Update only variants without existing gnomAD data
- Use PATCH requests to update only frequency fields
- Add proper logging and error handling"

# Push to remote
git push origin main
```

## Final Checks:

- [ ] Credentials removed from code
- [ ] No passwords in commit history
- [ ] README.md updated
- [ ] requirements.txt created
- [ ] Code tested and working
- [ ] Logging output is clean and informative
