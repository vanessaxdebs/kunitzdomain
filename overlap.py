# Check for overlap between seed and validation sets.
seed_ids = set()
with open("data/kunitz_seed.sto") as f:
    for line in f:
        if line.startswith("#=GS") or line.startswith(">"):
            name = line.split()[1] if line.startswith("#=GS") else line[1:].strip()
            seed_ids.add(name.split()[0])

val_ids = set()
with open("data/validation.fasta") as f:
    for line in f:
        if line.startswith(">"):
            name = line[1:].strip()
            val_ids.add(name.split()[0])

overlap = seed_ids & val_ids
print("Overlap:", overlap)
