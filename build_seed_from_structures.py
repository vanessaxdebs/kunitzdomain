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
# These are the PDB IDs and their respective chains for the Kunitz domains
pdb_ids = ["3TGI", "1BPI", "5PTI"]  # Add more if desired
chains = {"3TGI": "I", "1BPI": "A", "5PTI": "A"}

pdb_dir = "pdb_structures"
chain_dir = os.path.join(pdb_dir, "chains")
os.makedirs(chain_dir, exist_ok=True) # Ensure the directory for chain files exists
pdbl = PDBList() # Initialize PDBList for downloading

for pdb_id in pdb_ids:
    print(f"Downloading PDB structure '{pdb_id}'...")
    # Download the full PDB file
    pdbl.retrieve_pdb_file(pdb_id, pdir=pdb_dir, file_format="pdb")
    # Construct the path to the downloaded PDB file
    full_pdb_path = os.path.join(pdb_dir, f"pdb{pdb_id.lower()}.ent")
    
    # Define the output path for the selected chain file
    chain_file = os.path.join(chain_dir, f"{pdb_id}_{chains[pdb_id]}.pdb")
    
    # Use pdb_selchain to extract the specific chain
    print(f"Extracting chain {chains[pdb_id]} from {pdb_id}...")
    with open(chain_file, "w") as out:
        # subprocess.run executes an external command
        # pdb_selchain needs to be installed and in your system's PATH
        subprocess.run(["pdb_selchain", f"-{chains[pdb_id]}", full_pdb_path], stdout=out, check=True)
    print(f"Chain {chains[pdb_id]} saved to {chain_file}")

# Step 2: Align structures using MUSTANG
aligned_dir = os.path.join(pdb_dir, "aligned")
os.makedirs(aligned_dir, exist_ok=True) # Ensure the directory for aligned structures exists

# Define the common output prefix for MUSTANG's files (e.g., kunitz_alignment.fasta, kunitz_alignment.pdb)
output_prefix = os.path.join(aligned_dir, "kunitz_alignment")

# Define the path where the FASTA alignment output is expected (corrected from .fasta to .afasta)
aligned_fasta = os.path.join(aligned_dir, "kunitz_alignment.afasta") # Corrected from .fasta to .afasta based on MUSTANG output

# Construct the command for MUSTANG
cmd = [
    "mustang",
    "-i", # Input structures flag
    # List all the generated chain PDB files as inputs
    *[os.path.join(chain_dir, f"{pdb_id}_{chains[pdb_id]}.pdb") for pdb_id in pdb_ids],
    "-o", output_prefix, # Output identifier/prefix flag
    "-F", "fasta"         # Crucial: This tells MUSTANG to output in FASTA format
]

# Run the MUSTANG command
print("\nRunning MUSTANG:", " ".join(cmd))
# subprocess.run(check=True) will raise an error if MUSTANG fails
subprocess.run(cmd, check=True) 
print("MUSTANG alignment completed successfully.")

# Step 3: Convert aligned FASTA to Stockholm format
input_fasta = aligned_fasta # Use the path to the FASTA file MUSTANG just created
output_sto = "stockholm/kunitz_seed.sto"
os.makedirs("stockholm", exist_ok=True) # Ensure the directory for Stockholm file exists

print(f"\nConverting FASTA alignment from '{input_fasta}' to Stockholm format...")
# Read the FASTA alignment using Biopython's AlignIO
align = AlignIO.read(input_fasta, "fasta")

# Write the alignment to Stockholm format
with open(output_sto, "w") as out:
    AlignIO.write(align, out, "stockholm")

print(f"Structural Stockholm alignment created at {output_sto}")
# Write the alignment to Stockholm format
with open(output_sto, "w") as out:
    AlignIO.write(align, out, "stockholm")

print(f"Structural Stockholm alignment created at {output_sto}")
