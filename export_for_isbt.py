import os
import pandas as pd
import requests
from openpyxl.utils import get_column_letter
from openpyxl.styles import Alignment
import concurrent.futures
import re


# Define function to generate args
def generate_args():
    import argparse

    parser = argparse.ArgumentParser(description="Export for ISBT")
    parser.add_argument("--lead-url", type=str, required=True, help="Lead URL")
    parser.add_argument("--output", type=str, required=True, help="Output directory")
    parser.add_argument(
        "--threads", type=int, default=1, help="Number of threads to use"
    )

    return parser.parse_args()


# Function for exporting allele tables for each system
def export_allele_tables(system, lead_url, session, output_dir):

    # First get all the alleles for a system
    response = session.get(f"{lead_url}/allele/search?system_symbol={system}")

    # Check if request was successful
    if response.status_code == 200:
        data = response.json()
    else:
        raise Exception(f"Request failed with status code {response.status_code}")

    # Check there are alleles
    if len(data) == 0:
        print(f"No alleles found for system {system}")
        return

    # Get a list of allele ids
    allele_ids = [allele["id"] for allele in data]

    # Get the allele tables
    allele_data = []

    # Print a message saying fetching X alleles
    print(f"Fetching {len(allele_ids)} alleles for system {system}")
    # Create counter
    counter = 1
    for allele_id in allele_ids:
        # Print Counter
        print(f"Fetching {system} allele {counter} of {len(allele_ids)}")

        response = session.get(f"{lead_url}/allele/{allele_id}")

        # Check if request was successful
        if response.status_code == 200:
            allele_data.append(response.json())
        else:
            raise Exception(f"Request failed with status code {response.status_code}")

        counter += 1

    # Create a dataframe
    allele_df = pd.DataFrame(allele_data)

    # Filter for specific columns - id, isbt_allele, alternate_names , isbt_phenotype, reference_allele, sv_allele, null_allele, mod_allele, partial_allele, weak_allele, el_allele, notes, comment, system, gene, variants, genbanks, publications
    allele_df = allele_df[
        [
            "id",
            "isbt_allele",
            "alternate_names",
            "isbt_phenotype",
            "reference_allele",
            "sv_allele",
            "null_allele",
            "mod_allele",
            "partial_allele",
            "weak_allele",
            "el_allele",
            "notes",
            "comment",
            "variants",
            "genbanks",
            "publications",
            "gene",
        ]
    ]

    # Replace these columns that contain booleans with 1 or 0, replace na with 0 reference_allele	sv_allele	null_allele	mod_allele	partial_allele	weak_allele	el_allele
    # Make a list of columns to replace
    bool_cols = [
        "reference_allele",
        "sv_allele",
        "null_allele",
        "mod_allele",
        "partial_allele",
        "weak_allele",
        "el_allele",
    ]

    # Replace True with 1 and False with 0
    allele_df[bool_cols] = (
        allele_df[bool_cols].replace({True: 1, False: 0}).infer_objects(copy=False)
    )

    # Replace NA with 0
    allele_df[bool_cols] = allele_df[bool_cols].fillna(0).infer_objects(copy=False)

    # Get variant from each row
    for idx, row in allele_df.iterrows():
        # Get the variants
        variants = row["variants"]

        # v_list
        v_list = []
        p_list = []
        r_list = []
        chromosomal_list = []
        #! COLLECT EXONS HERE
        exon_intron_list = []

        # If there are variants
        if len(variants) > 0:

            # First sort variants by hgvs_transcript - NM_017436.7:c.[-188]+3010G>T
            try:
                variants = sorted(
                    variants,
                    key=lambda x: float(
                        re.search(r"c\.(-?\d+)", x["hgvs_transcript"]).group(1)
                    ),
                )
            except:
                pass

            for variant in variants:

                # Get the variant transcript
                variant_name = variant["hgvs_transcript"]

                # Get the variant protein
                variant_protein = variant["hgvs_predicted_protein"]

                # Get the variant rsid
                variant_rsid = variant["rsid"]

                # Get hgvs_genomic_grch38
                variant_genomic = variant["hgvs_genomic_grch38"]

                # Get the exon intron
                exon_intron = variant["exon"]

                # if any are "" or NA then replace with "-"
                if variant_name == "":
                    variant_name = "-"
                if variant_protein == "":
                    variant_protein = "-"
                if variant_rsid == "":
                    variant_rsid = "-"
                if variant_genomic == "":
                    variant_genomic = "-"

                # Append to list
                v_list.append(variant_name)
                p_list.append(variant_protein)
                r_list.append(variant_rsid)
                chromosomal_list.append(variant_genomic)
                exon_intron_list.append(exon_intron)

        v_list = [str(x) if pd.notna(x) else "-" for x in v_list]
        p_list = [str(x) if pd.notna(x) else "-" for x in p_list]
        r_list = [str(x) if pd.notna(x) else "-" for x in r_list]
        chromosomal_list = [str(x) if pd.notna(x) else "-" for x in chromosomal_list]
        exon_intron_list = [str(x) if pd.notna(x) else "-" for x in exon_intron_list]

        # For each string in v_list and p_list split on ":" if in there and take the second element
        for i in range(len(v_list)):
            if ":" in v_list[i]:
                v_list[i] = v_list[i].split(":")[1]

        for i in range(len(p_list)):
            if ":" in p_list[i]:
                p_list[i] = p_list[i].split(":")[1]

        # Add to the dataframe
        allele_df.loc[idx, "variant"] = "\n".join(v_list)
        allele_df.loc[idx, "protein_variant"] = "\n".join(p_list)
        allele_df.loc[idx, "rsid"] = "\n".join(r_list)
        allele_df.loc[idx, "genomic_variant"] = "\n".join(chromosomal_list)
        allele_df.loc[idx, "exon_intron"] = "\n".join(exon_intron_list)

    # For publications
    pubs = []
    for publication in allele_df["publications"]:
        pub_list = []
        # Check if publication is None/NA
        if publication is not None:
            # Sort publications by type PMID before Abstract, then alphanumeric
            publication = sorted(
                publication,
                key=lambda x: (x.get("type", "-"), x.get("identifier", "-")),
            )

            for pub in publication:
                try:
                    if pub is not None:
                        # if type is pmid then get the identifier
                        if pub.get("type") == "pmid":
                            identifier = pub.get("identifier", "-")
                            pub_list.append(f"PMID:{identifier}")
                        # If type is abstract, get the identifier and the citation
                        elif pub.get("type") == "abstract":
                            identifier = pub.get("identifier", "-")
                            citation = pub.get("citation", "-")
                            pub_list.append(f"{identifier}: {citation}")
                        # For all other types
                        else:
                            citation = pub.get("citation", "-")
                            pub_list.append(citation)
                except AttributeError:
                    # Handle case where pub is not a dictionary
                    pub_list.append("-")

        # If no valid publications were found, append a dash
        pubs.append("\n".join(pub_list) if pub_list else "-")

    # Add to the dataframe
    allele_df["publications"] = pubs

    # For Genbanks
    genbanks = []
    for genbank in allele_df["genbanks"]:

        gb_list = []
        # Check if genbank is None/NA
        if genbank is not None:
            # Sort genbanks by accession number
            genbank = sorted(genbank, key=lambda x: x.get("accession", "-"))

            for gb in genbank:
                try:
                    if gb is not None:
                        accession = gb.get("accession", "-")
                        gb_list.append(accession)
                except AttributeError:
                    # Handle case where gb is not a dictionary
                    gb_list.append("-")

        # If no valid genbanks were found, append a dash
        genbanks.append("\n".join(gb_list) if gb_list else "-")

    # Add to the dataframe
    allele_df["genbanks"] = genbanks

    # For alternative names, join the lists with ;
    allele_df["alternate_names"] = allele_df["alternate_names"].apply(
        lambda x: ";".join(x) if x else "-"
    )

    # Replace gene column with gene name if available
    gene_names = []
    for gene in allele_df["gene"]:
        try:
            gene_names.append(gene["name"])
        except:
            gene_names.append("-")

    allele_df["gene"] = gene_names

    # Order the columns
    allele_df = allele_df[
        [
            "id",
            "isbt_phenotype",
            "isbt_allele",
            "alternate_names",
            "reference_allele",
            "variant",
            "exon_intron",
            "protein_variant",
            "genomic_variant",
            "rsid",
            "gene",
            "genbanks",
            "publications",
            "null_allele",
            "mod_allele",
            "partial_allele",
            "weak_allele",
            "el_allele",
            "sv_allele",
            "notes",
            "comment",
        ]
    ]

    # Rename columns
    allele_df.rename(
        columns={
            "id": "Database ID",
            "isbt_phenotype": "Phenotype",
            "isbt_allele": "ISBT Allele",
            "alternate_names": "Alternate Names",
            "reference_allele": "Reference Allele",
            "variant": "DNA change",
            "exon_intron": "Exon/Intron",
            "protein_variant": "Protein change",
            "genomic_variant": "Genomic Variant",
            "rsid": "RSID",
            "gene": "Gene",
            "genbanks": "GenBank",
            "publications": "Publications",
            "sv_allele": "Structural Variant Allele",
            "null_allele": "Null Allele",
            "mod_allele": "Modified Allele",
            "partial_allele": "Partial Allele",
            "weak_allele": "Weak Allele",
            "el_allele": "Elution Allele",
            "notes": "Annotation Notes",
            "comment": "Comments",
        },
        inplace=True,
    )

    # Save to Excel with auto-adjusted columns, rows, and text wrapping
    with pd.ExcelWriter(f"{output_dir}/{system}.xlsx", engine="openpyxl") as writer:
        allele_df.to_excel(writer, index=False)
        worksheet = writer.sheets["Sheet1"]

        # Enable text wrapping for all cells and adjust column widths
        for column in worksheet.columns:
            max_length = 0
            column_letter = get_column_letter(column[0].column)

            for cell in column:
                try:
                    # Enable text wrapping for each cell
                    cell.alignment = Alignment(wrap_text=True)

                    if cell.value:
                        # Count the maximum number of lines and length
                        length = max(
                            len(str(line)) for line in str(cell.value).split("\n")
                        )
                        if length > max_length:
                            max_length = length
                except:
                    pass

            adjusted_width = max_length + 2
            worksheet.column_dimensions[column_letter].width = adjusted_width

        # Auto-adjust row heights
        for row in worksheet.rows:
            max_height = 1
            for cell in row:
                try:
                    if cell.value:
                        # Count number of line breaks
                        num_lines = str(cell.value).count("\n") + 1
                        if num_lines > max_height:
                            max_height = num_lines
                except:
                    pass
            # Set row height based on number of lines
            worksheet.row_dimensions[row[0].row].height = max_height * 15


def get_system_list(lead_url, session):

    # Get data from lead URL
    response = session.get(f"{lead_url}/system")

    # Check if request was successful
    if response.status_code == 200:
        data = response.json()
    else:
        raise Exception(f"Request failed with status code {response.status_code}")

    # Get list of system names
    systems = [system["symbol"] for system in data]

    print(f"Found {len(systems)} systems: {systems}")

    return list(set(systems))


# Define main function
def main():

    args = generate_args()
    # Get absolute path of output directory and make if does not exist
    output_dir = os.path.abspath(args.output)
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Creat a session with the lead URL
    session = requests.Session()

    # Get list of systems
    systems = get_system_list(args.lead_url, session)
    # Print number of systems
    print(f"Found {len(systems)} systems")

    # Try just for KEL
    export_allele_tables("FY", args.lead_url, session, output_dir)

    # Do in parallel
    # with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
    #     futures = [
    #         executor.submit(
    #             export_allele_tables, system, args.lead_url, session, output_dir
    #         )
    #         for system in systems
    #     ]
    #     for future in concurrent.futures.as_completed(futures):
    #         try:
    #             future.result()
    #         except Exception as exc:
    #             print(f"Generated an exception: {exc}")


if __name__ == "__main__":
    main()
