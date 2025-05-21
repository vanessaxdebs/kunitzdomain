#!/usr/bin/env python3
"""
Novelty check for Kunitz HMM hits:
Flags hits not already annotated as Kunitz (by keyword, Pfam, or domain description) in UniProt.
Works with HMMER .tblout with IDs like sp|P84875|PCPI_SABMA (extracts accession directly).
"""

import os
from pathlib import Path
import requests
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

def extract_uniprot_accession(hit_id):
    """
    For a hit like sp|P84875|PCPI_SABMA, returns P84875.
    For a hit like APLP2_HUMAN/309-361, returns APLP2_HUMAN (though this is not a real accession).
    """
    if hit_id.startswith("sp|") or hit_id.startswith("tr|"):
        parts = hit_id.split("|")
        if len(parts) >= 2:
            return parts[1]
    # Fallback: return the part before '/' (if present)
    return hit_id.split("/")[0]

def load_hmm_hits(hmm_results_file):
    """
    Loads all unique accessions from the results .tblout file.
    """
    hits = set()
    with open(hmm_results_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            raw_id = line.split()[0]
            accession = extract_uniprot_accession(raw_id)
            hits.add(accession)
    return hits

def is_kunitz_annotated(uniprot_acc):
    """
    Query UniProt API to check for Kunitz annotation in keywords, Pfam, or domain features.
    """
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
    accessions = load_hmm_hits(tbl_file)
    print(f"Total unique accessions found: {len(accessions)}")
    novel = []
    checked = 0
    for accession in sorted(accessions):
        checked += 1
        print(f"[{checked}/{len(accessions)}] Checking {accession}...", end=" ")
        time.sleep(0.2)  # Be kind to UniProt servers
        if not is_kunitz_annotated(accession):
            print("NOVEL")
            novel.append(accession)
        else:
            print("known")
        time.sleep(0.2)
    print(f"\nNovel candidate hits (not annotated as Kunitz): {len(novel)}")
    for n in novel:
        print(n)
    with open("results/novel_kunitz_candidates.txt", "w") as out:
        for n in novel:
            out.write(n + "\n")
    print("\nDone. See results/novel_kunitz_candidates.txt for the list.")

if __name__ == "__main__":
    main()
