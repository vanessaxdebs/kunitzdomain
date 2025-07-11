import gzip
from Bio import SeqIO

# 1. Extract training accessions from kunitz_seed.sto
training_accessions = set()
with open("data/kunitz_seed.sto") as sto:
    for line in sto:
        if line.startswith("#=GS"):
            parts = line.strip().split()
            if "AC" in parts:
                ac_index = parts.index("AC")
                acc = parts[ac_index + 1].split('.')[0]
                training_accessions.add(acc)

# 2. Read gzipped UniProt Kunitz domain FASTA and filter out training set
input_fasta = "data/uniprotkb_PF00014_2025_07_11.fasta.gz"
output_fasta = "data/validation_positives.fasta"

with gzip.open(input_fasta, "rt") as in_f, open(output_fasta, "w") as out_f:
    for record in SeqIO.parse(in_f, "fasta"):
        # UniProt FASTA headers look like: >sp|O17644|O17644_CAEEL ...
        if "|" in record.id:
            acc = record.id.split("|")[1]
        else:
            acc = record.id.split("_")[0]
        if acc not in training_accessions:
            SeqIO.write(record, out_f, "fasta")

print(f"Filtered validation set written to {output_fasta}")
