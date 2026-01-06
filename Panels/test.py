import requests

def entrez_species(entrez_id):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        f"?db=gene&id={entrez_id}&retmode=json"
    )
    r = requests.get(url, timeout=10)
    data = r.json()

    uid = str(entrez_id)
    if uid not in data["result"]:
        return None

    return data["result"][uid].get("organism", {}).get("scientificname")


# test
for test_id in [565089, 526146, 930, 12478]:
    print(test_id, "→", entrez_species(test_id))
