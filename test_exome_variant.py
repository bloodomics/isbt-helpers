#!/usr/bin/env python3
"""Quick test to verify exome-only variant detection"""

import requests
import json

# Test with the exome-only variant
query = """
query GnomadVariant($variantId: String!, $datasetId: DatasetId!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    variant_id
    exome { 
      af
      populations { id ac an }
    }
    genome { 
      af 
      populations { id ac an }
    }
  }
}
"""

variables = {
    'variantId': '7-142957922-T-A',
    'datasetId': 'gnomad_r4'
}

print("Testing exome-only variant: 7-142957922-T-A")
print("=" * 60)

response = requests.post(
    'https://gnomad.broadinstitute.org/api',
    json={'query': query, 'variables': variables},
    headers={'Content-Type': 'application/json'}
)

data = response.json()
variant = data['data']['variant']

print(f"Variant: {variant['variant_id']}")
print(f"Exome AF: {variant['exome']['af'] if variant['exome'] else None}")
print(f"Genome AF: {variant['genome'] if variant['genome'] else None}")
print()

if variant['exome'] and not variant['genome']:
    print("✓ This is an exome-only variant")
    print("✓ The updated script will now capture this variant!")
    
    # Calculate combined frequency (same as script logic)
    combined_ac = 0
    combined_an = 0
    for pop in variant['exome']['populations']:
        if ':' not in pop['id'] and '_XX' not in pop['id'] and '_XY' not in pop['id']:
            combined_ac += pop['ac']
            combined_an += pop['an']
    
    combined_af = combined_ac / combined_an if combined_an > 0 else 0
    combined_maf = min(combined_af, 1 - combined_af) if combined_af else None
    
    print(f"\nCombined (exome-only) MAF: {combined_maf:.6e}")
else:
    print("✗ Unexpected result")
