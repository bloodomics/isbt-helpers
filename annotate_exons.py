import pandas as pd
import requests
import json
import numpy as np

lead_url = "https://api.blooddatabase.org"

# Ask user for email and password at terminal
email = input("Enter your email: ")
password = input("Enter your password: ")


# Make a requests session
session = requests.session()

# Post to login
response = session.post(
    "{}/auth/login".format(lead_url),
    json={"email": f"{email}", "password": f"{password}"},
)

# Ceck the response
print(json.loads(response.text))

# Get a list of variants
variants = session.get("{}/variant".format(lead_url))
variants = json.loads(variants.text)


# Define a process variant function
def annotate_exon(variant):

    # If the exon is already set, skip
    if variant["exon"]:
        return

    # if there is no hgvs_transcript, skip
    if not variant["hgvs_transcript"]:
        return

    hgvs_transcript = variant["hgvs_transcript"]
    transcript = hgvs_transcript.split(":")[0]

    vv_url = f"https://rest.variantvalidator.org/VariantValidator/variantvalidator/GRCh38/{hgvs_transcript}/{transcript}"

    response = session.get(vv_url)

    if response.status_code == 200:
        response = json.loads(response.text)[hgvs_transcript]

        # get the variant_exonic_positions
        variant_exonic_positions = response["variant_exonic_positions"]

        # Get the hg38 start_exon (look for chromosome ending in .12)
        start_exon = None
        latest_version = 0
        # Check all keys which start with NC , get the highest value NC_xxxx.{value}
        for key in variant_exonic_positions.keys():
            if key.startswith("NC"):
                version = int(key.split(".")[1])
                if version > latest_version:
                    latest_version = version
                    start_exon = variant_exonic_positions[key]["start_exon"]
                    end_exon = variant_exonic_positions[key]["end_exon"]

        # IF the start exon and end exon are the same then we can just use the start exon
        if start_exon == end_exon:
            new_exon_value = start_exon
        else:
            new_exon_value = f"{start_exon}-{end_exon}"

        # check if sta rt exon is None
        if not new_exon_value:
            print(
                f"No exon found for variant - id{variant['id']} transcript: {variant['hgvs_transcript']}"
            )
            return
        # Patch the variant to update the exon
        patch_url = f"{lead_url}/variant/{variant['id']}"
        patch_data = {"exon": new_exon_value}
        response = session.patch(patch_url, json=patch_data)

        # print(response.text)
    else:

        print(
            f"No exon found for variant - id{variant['id']} transcript: {variant['hgvs_transcript']}"
        )


import time

to_go = len(variants)
# Apply this to all variants but limit to 4 requests per second
for variant in variants:
    to_go -= 1
    print(f"{to_go} to go")
    annotate_exon(variant)
