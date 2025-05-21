import os
from pathlib import Path
import requests
import re
import time

def get_latest_results_tbl():
    """Find the most recent hmmsearch_swissprot.tbl in results/run_* subdirs."""
    results_dir = Path("results")
    subdirs = [d for d in results_dir.iterdir() if d.is_dir() and d.name.startswith("run_")]
    if not subdirs:
        raise FileNotFoundError("No run_* subdirectories in results/")
    latest_subdir = max(subdirs, key=os.path.getmtime)
    tbl_file = latest_subdir / "hmmsearch_swissprot.tbl"
    if not tbl_file.exists():
        raise FileNotFoundError(f"{tbl_file} does not exist.")
    return tbl_file

def extract_entry_name(hit_id):
    # e.g., "APLP2_HUMAN/309-361" -> "APLP2_HUMAN"
    return hit_id.split('/')[0]

def extract_uniprot_accession(hit_id):
    # Handles IDs like 'sp|F6ULY1|WFC6B_MOUSE' or just 'F6ULY1'
    if hit_id.startswith("sp|") or hit_id.startswith("tr|"):
        # Standard UniProt FASTA header
        parts = hit_id.split("|")
        if len(parts) >= 2:
            return parts[1]
    # If already an accession
    return hit_id.split("/")[0]  # Also strips any residue range

def load_hmm_hits(hmm_results_file):
    hits = set()
    with open(hmm_results_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            raw_id = line.split()[0]
            entry_name = extract_entry_name(raw_id)
            hits.add(entry_name)
    return hits

def is_kunitz_annotated(uniprot_acc):
    """Query UniProt API to check for Kunitz annotation in the entry."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_acc}.json"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            return False  # Could not fetch or not found
        data = r.json()
        # Check features for 'Kunitz'
        for feat in data.get("features", []):
            if "Kunitz" in feat.get("description", ""):
                return True
        # Check keywords for 'Kunitz'
        for keyword in data.get("keywords", []):
            if "Kunitz" in keyword.get("value", ""):
                return True
        # Check Pfam domains for PF00014
        for dbref in data.get("dbReferences", []):
            if dbref.get("type") == "Pfam" and dbref.get("id") == "PF00014":
                return True
        return False
    except Exception as e:
        print(f"Warning: failed to fetch or parse UniProt entry for {uniprot_acc}: {e}")
        return False

def main():
    tbl_file = get_latest_results_tbl()
    print(f"Using results file: {tbl_file}")
    entry_names = load_hmm_hits(tbl_file)
    print(f"Total unique entry names found: {len(entry_names)}")
    novel = []
    checked = 0
    for entry_name in entry_names:
        checked += 1
        print(f"[{checked}/{len(entry_names)}] Checking {entry_name}...", end=" ")
        accession = entry_name_to_accession(entry_name)
        if not accession:
            print("COULD NOT FIND ACCESSION -- SKIPPING")
            continue
        time.sleep(0.2)  # Be kind to UniProt servers
        if not is_kunitz_annotated(accession):
            print("NOVEL")
            novel.append(f"{entry_name} ({accession})")
        else:
            print("known")
        time.sleep(0.2)  # Be kind to UniProt servers
    print(f"\nNovel candidate hits (not annotated as Kunitz): {len(novel)}")
    for n in novel:
        print(n)
    with open("results/novel_kunitz_candidates.txt", "w") as out:
        for n in novel:
            out.write(n + "\n")
    print("\nDone. See results/novel_kunitz_candidates.txt for the list.")

if __name__ == "__main__":
    main()
