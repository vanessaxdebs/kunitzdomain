import requests

def load_hmm_hits(hmm_results_file):
    hits = set()
    with open(hmm_results_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            hits.add(line.split()[0])
    return hits

def is_kunitz_annotated(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    r = requests.get(url)
    if r.status_code != 200:
        return False  # Could not fetch or not found
    data = r.json()
    # Look in features or keywords for 'Kunitz'
    for feat in data.get("features", []):
        if "Kunitz" in feat.get("description", ""):
            return True
    for keyword in data.get("keywords", []):
        if "Kunitz" in keyword.get("value", ""):
            return True
    # Check Pfam annotations (if available)
    for dbref in data.get("dbReferences", []):
        if dbref.get("type") == "Pfam" and dbref.get("id") == "PF00014":
            return True
    return False

def main():
    hits = load_hmm_hits("results/hmmsearch_swissprot.tbl")
    novel = []
    for uniprot_id in hits:
        if not is_kunitz_annotated(uniprot_id):
            novel.append(uniprot_id)
    print(f"Novel candidate hits (not annotated as Kunitz):")
    for n in novel:
        print(n)
    with open("results/novel_kunitz_candidates.txt", "w") as out:
        for n in novel:
            out.write(n + "\n")

if __name__ == "__main__":
    main()
