#!/usr/bin/env python3
"""
Kunitz-type Protease Inhibitor Domain - Profile HMM Project

Bioinformatics Course Project
Author: Vanessa El Debs 
Inspired by: Prof. Capriotti's course

Main Aims:
- Build a profile HMM for the Kunitz domain
- Validate the model and annotate SwissProt sequences
"""

import os
import subprocess
import sys
from pathlib import Path
import yaml
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# === 1. CONFIGURATION ===
def load_config(config_file="config.yaml"):
    """Load parameters from YAML config file"""
    try:
        with open(config_file) as f:
            config = yaml.safe_load(f)
        
        # Create output directory with timestamp
        run_id = datetime.now().strftime("%Y%m%d_%H%M")
        output_dir = Path(config['output_dir']) / f"run_{run_id}"
        output_dir.mkdir(parents=True, exist_ok=True)
        config['output_dir'] = output_dir
        
        return config
    except Exception as e:
        print(f"Error loading config: {e}")
        sys.exit(1)

# Load config at startup
CONFIG = load_config()

# === 2. STRUCTURAL ANALYSIS ===
def fetch_pdb(pdb_id, out_dir):
    """Download PDB file with error handling"""
    try:
        from Bio.PDB import PDBList
        pdbl = PDBList()
        pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir=out_dir, file_format='pdb')
        print(f"PDB file saved at {pdb_file}")
        return pdb_file
    except Exception as e:
        print(f"Failed to fetch PDB: {e}")
        sys.exit(1)

# === 3. HMM BUILDING & VALIDATION === 
def build_hmm(seed_alignment, output_dir):
    """Build HMM with better error handling"""
    hmm_file = output_dir / "kunitz.hmm"
    cmd = ["hmmbuild", str(hmm_file), str(seed_alignment)]
    
    print(f"Building HMM: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"HMM built successfully at {hmm_file}")
        return hmm_file
    except subprocess.CalledProcessError as e:
        print(f"HMM build failed. Error:\n{e.stderr}")
        sys.exit(1)

# === 4. METRICS & VISUALIZATION ===
def plot_confusion(metrics, outdir):
    """Improved confusion matrix plot"""
    cm = np.array([[metrics["TP"], metrics["FP"]], 
                  [metrics["FN"], metrics["TN"]]])
    
    fig, ax = plt.subplots(figsize=(6,6))
    im = ax.imshow(cm, cmap='Blues')
    
    # Add text annotations
    for i in range(2):
        for j in range(2):
            ax.text(j, i, cm[i,j], 
                   ha='center', va='center',
                   color='black', fontsize=14)
    
    # Labels and title
    ax.set_xlabel("Predicted", fontsize=12)
    ax.set_ylabel("Actual", fontsize=12)
    ax.set_title("Confusion Matrix", fontsize=14)
    ax.set_xticks([0,1])
    ax.set_yticks([0,1])
    ax.set_xticklabels(["Positive", "Negative"])
    ax.set_yticklabels(["Positive", "Negative"])
    
    plt.colorbar(im)
    plt.tight_layout()
    plot_file = outdir / "confusion_matrix.png"
    plt.savefig(plot_file, dpi=150)
    plt.close()
    print(f"Saved confusion matrix to {plot_file}")

# === 5. MAIN PIPELINE ===
def main():
    print(f"\n{'='*50}")
    print("Kunitz-type Protease Inhibitor HMM Project")
    print(f"{'='*50}\n")
    
    try:
        # 1. Structural analysis
        print("[1/6] Structural Analysis")
        fetch_pdb(CONFIG['pdb_id'], CONFIG['output_dir'])
        
        # 2. Seed alignment check
        print("\n[2/6] Checking Seed Alignment")
        if not Path(CONFIG['seed_alignment']).exists():
            raise FileNotFoundError(f"Seed alignment not found at {CONFIG['seed_alignment']}")
        
        # 3. HMM building
        print("\n[3/6] Building Profile HMM")
        hmm_file = build_hmm(CONFIG['seed_alignment'], CONFIG['output_dir'])
        
        # 4. Validation
        print("\n[4/6] Running Validation")
        tblout = run_hmmsearch(hmm_file, CONFIG['validation_fasta'], CONFIG['output_dir'])
        predicted = parse_tblout(tblout)
        positives, negatives = load_labels(CONFIG['validation_labels'])
        metrics = evaluate_performance(predicted, positives, negatives)
        plot_confusion(metrics, CONFIG['output_dir'])
        
        # 5. SwissProt annotation
        print("\n[5/6] Annotating SwissProt")
        swiss_tbl = annotate_swissprot(hmm_file, CONFIG['swissprot_fasta'], CONFIG['output_dir'])
        analyze_swissprot(swiss_tbl)
        
        # 6. HMM logo
        print("\n[6/6] Generating HMM Logo")
        run_hmmlogo(hmm_file, CONFIG['output_dir'])
        
        print(f"\n{'='*50}")
        print(f"Pipeline completed successfully!")
        print(f"Results saved to: {CONFIG['output_dir']}")
        print(f"{'='*50}")
        
    except Exception as e:
        print(f"\nERROR: Pipeline failed - {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
