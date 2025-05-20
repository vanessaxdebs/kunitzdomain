#!/usr/bin/env python3
"""
Build a Stockholm-format seed alignment for the Kunitz domain
based on structural alignment of known 3D PDB chains.
"""

import os
import subprocess
from Bio import AlignIO, SeqIO
from Bio.PDB import PDBList

# Step 1: Download PDB files
pdb_ids = ["3TGI", "1BPI", "5PTI"]  # Add more if desired
chains = {"3TGI": "I", "1BPI": "A", "5PTI": "A"}

pdb_dir = "pdb_structures"
chain_dir = os.path.join(pdb_dir, "chains")
os.makedirs(chain_dir, exist_ok=True)
pdbl = PDBList()

for pdb_id in pdb_ids:
    pdbl.retrieve_pdb_file(pdb_id, pdir=pdb_dir, file_format="pdb")
    full_pdb_path = os.path.join(pdb_dir, f"pdb{pdb_id.lower()}.ent")
    chain_file = os.path.join(chain_dir, f"{pdb_id}_{chains[pdb_id]}.pdb")
    with open(chain_file, "w") as out:
        subprocess.run(["pdb_selchain", f"-{chains[pdb_id]}", full_pdb_path], stdout=out)

# Step 2: Align with MUSTANG
aligned_dir = os.path.join(pdb_dir, "aligned")
os.makedirs(aligned_dir, exist_ok=True)
aligned_fasta = os.path.join(aligned_dir, "kunitz_alignment.fasta")

cmd = [
    "mustang",
    "-i",
    *[os.path.join(chain_dir, f"{pdb_id}_{chains[pdb_id]}.pdb") for pdb_id in pdb_ids],
    "-o", os.path.join(aligned_dir, "kunitz_alignment")
]
print("Running MUSTANG:", " ".join(cmd))
subprocess.run(cmd, check=True)

# Step 3: Convert aligned FASTA to Stockholm
input_fasta = os.path.join(aligned_dir, "kunitz_alignment.fasta")
output_sto = "stockholm/kunitz_seed.sto"
os.makedirs("stockholm", exist_ok=True)

align = AlignIO.read(input_fasta, "fasta")
with open(output_sto, "w") as out:
    AlignIO.write(align, out, "stockholm")

print(f" Structural Stockholm alignment created at {output_sto}")
