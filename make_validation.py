import gzip

# Step 1: Collect training sequence IDs
training_ids = set()
with open("data/kunitz_seed.sto") as f:
    for line in f:
        if line.startswith("#=GS"):
            training_ids.add(line.split()[1])
        elif line and not line.startswith(("#", "//", "\n")):
            training_ids.add(line.split()[0])

# Step 2: Extract and filter validation sequences
seqs = {}
with gzip.open("data/PF00014.alignment.seed.gz", "rt") as f:
    for line in f:
        if line.startswith("#") or line.startswith("//") or not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) == 2:
            seq_id, seq = parts
            if seq_id not in training_ids:
                seqs.setdefault(seq_id, "")
                seqs[seq_id] += seq

# Step 3: Write filtered FASTA
with open("data/validation.fasta", "w") as fasta:
    for seq_id, seq in seqs.items():
        fasta.write(f">{seq_id}\n{seq}\n")

# Step 4: Write labels (all 1s)
with open("data/validation_labels.txt", "w") as labels:
    for seq_id in seqs.keys():
        labels.write(f"{seq_id}\t1\n")
