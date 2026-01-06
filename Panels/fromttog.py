import pandas as pd
import requests
import time

print("RUNNING FILE:", __file__)

ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/id/{}"

def transcript_to_gene(ensmust_id):
    """Query Ensembl REST API to convert ENSMUST → ENSMUSG"""
    url = ENSEMBL_LOOKUP_URL.format(ensmust_id)
    headers = {"Content-Type": "application/json"}

    r = requests.get(url, headers=headers)
    if not r.ok:
        return None

    data = r.json()
    return data.get("Parent")  # this is the ENSMUSG


def convert_csv(csv_path, output_path):
    df = pd.read_csv(csv_path, index_col=False)

    print("Columns:", df.columns.tolist())
    print(df.head())
    print("ID column sample:", df["id"].head().tolist())

    if "id" not in df.columns:
        raise ValueError("Column 'id' not found in CSV")

    # Clean transcript IDs
    transcripts = (
        df["id"]
        .dropna()
        .astype(str)
        .str.replace(r"\.\d+$", "", regex=True)
        .unique()
    )

    gene_ids = set()

    for t in transcripts:
        print("Processing transcript ID:", t)
        gene = transcript_to_gene(t)
        if gene:
            gene_ids.add(gene)
        time.sleep(0.1)  # be nice to Ensembl servers

    out_df = pd.DataFrame({"Ensemble_Id": sorted(gene_ids)})
    out_df.to_csv(output_path, index=False)

    print(f"Saved {len(out_df)} gene IDs to {output_path}")


print("SCRIPT PARTITO")
convert_csv(
    csv_path="raw_panels/202312_codebook_0_mouse-brain-final_CP1489.csv",
    output_path="202312_codebook_0_mouse-brain-final_CP1489.csv"
)
print("SCRIPT PARTITO")