import os

REQUIRED_FILES = {
    "data/kunitz_seed.sto": "# STOCKHOLM",
    "data/validation.fasta": ">",
    "data/validation_labels.txt": "",
}

def check_file(path, must_startwith):
    if not os.path.isfile(path):
        print(f"ERROR: {path} is missing.")
        return False
    if os.path.getsize(path) == 0:
        print(f"ERROR: {path} is empty.")
        return False
    if must_startwith:
        with open(path) as f:
            first_line = f.readline()
            if not first_line.startswith(must_startwith):
                print(f"ERROR: {path} does not start with {must_startwith}. Found: {first_line.strip()}")
                return False
    print(f"OK: {path}")
    return True

if __name__ == "__main__":
    ok = True
    for file, startswith in REQUIRED_FILES.items():
        ok = check_file(file, startswith) and ok
    if not ok:
        print("Some data files are missing or invalid. Please follow the README to set them up correctly.")
        exit(1)
    else:
        print("All data files are present and look valid.")


