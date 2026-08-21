import urllib.request
from pathlib import Path

# Official PubChem CIDs for the 5 drugs
DRUG_CIDS = {
    "gefitinib": 3088125,
    "erlotinib": 176870,
    "afatinib": 10184653,
    "dacomitinib": 11511120,
    "osimertinib": 71496458,
}

output_dir = Path("data/raw/ligands")
output_dir.mkdir(parents=True, exist_ok=True)

for drug_name, cid in DRUG_CIDS.items():
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF/?record_type=3d"
    file_path = output_dir / f"{drug_name}.sdf"

    print(f"Downloading official 3D SDF for {drug_name} (CID: {cid})...")
    urllib.request.urlretrieve(url, file_path)

print("\nAll 5 official 3D drug datasets downloaded successfully to data/raw/ligands/")