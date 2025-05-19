#!/usr/bin/env python3
"""
Kunitz-type Protease Inhibitor Domain - Profile HMM Project

Complete implementation with all functions and type hints
"""

import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, Set, Tuple
import yaml
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# === 1. CONFIGURATION ===
def find_config() -> Path:
    """Locate config.yaml whether running from bin/ or project root"""
    script_dir = Path(__file__).parent
    for path in [
        script_dir / "../config/config.yaml",  # If running from bin/
        script_dir / "config/config.yaml",     # If running from root
        Path("config/config.yaml")             # Fallback
    ]:
        if path.exists():
            return path.resolve()
    raise FileNotFoundError("config.yaml not found!")

def load_config(config_file: Path = None) -> Dict:
    """Load and validate configuration from YAML file
    
    Args:
        config_file: Path to YAML configuration file
        
    Returns:
        Dictionary containing all configuration parameters
        
    Raises:
        SystemExit: If config file is invalid or missing required fields
    """
    if config_file is None:
        config_file = find_config()
    
    required_fields = {
        'output_dir': str,
        'seed_alignment': str,
        'validation_fasta': str,
        'validation_labels': str,
        'swissprot_fasta': str,
        'e_value_cutoff': float,
        'pdb_id': str
    }
    
    try:
        with open(config_file) as f:
            config = yaml.safe_load(f)
            
        # Validate required fields
        for field, field_type in required_fields.items():
            if field not in config:
                raise ValueError(f"Missing required config field: {field}")
            if not isinstance(config[field], field_type):
                raise ValueError(f"Invalid type for {field}, expected {field_type.__name__}")
                
        # Create output directory with timestamp
        run_id = datetime.now().strftime("%Y%m%d_%H%M")
        output_dir = Path(config['output_dir']) / f"run_{run_id}"
        output_dir.mkdir(parents=True, exist_ok=True)
        config['output_dir'] = output_dir
        
        return config
        
    except Exception as e:
        print(f"Configuration error: {e}", file=sys.stderr)
        sys.exit(1)

# === 2. CORE FUNCTIONS ===
def run_hmmsearch(hmm_file: Path, fasta_file: Path, output_dir: Path, tag: str = "validation") -> Path:
    """Run hmmsearch and save table output
    
    Args:
        hmm_file: Path to HMM file
        fasta_file: Path to FASTA file to search
        output_dir: Directory to save results
        tag: Identifier for output files
        
    Returns:
        Path to results table file
        
    Raises:
        SystemExit: If hmmsearch fails
    """
    tblout = output_dir / f"hmmsearch_{tag}.tbl"
    cmd = [
        "hmmsearch",
        "--tblout", str(tblout),
        "-E", str(CONFIG['e_value_cutoff']),
        str(hmm_file), str(fasta_file)
    ]
    
    try:
        print(f"Running hmmsearch: {' '.join(cmd)}")
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return tblout
    except subprocess.CalledProcessError as e:
        print(f"hmmsearch failed. Error:\n{e.stderr}", file=sys.stderr)
        sys.exit(1)

def parse_tblout(tbl_file: Path) -> Set[str]:
    """Parse HMMER tblout file and extract hits
    
    Args:
        tbl_file: Path to HMMER output file
        
    Returns:
        Set of sequence IDs that were hits
    """
    hits = set()
    try:
        with open(tbl_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) > 0:
                    hits.add(parts[0])
        return hits
    except Exception as e:
        print(f"Error parsing tblout file: {e}", file=sys.stderr)
        sys.exit(1)

def load_labels(label_file: Path) -> Tuple[Set[str], Set[str]]:
    """Load validation labels from tab-separated file
    
    Args:
        label_file: Path to label file (format: seqid[tab]1/0)
        
    Returns:
        Tuple of (positive_ids, negative_ids)
    """
    pos, neg = set(), set()
    try:
        with open(label_file) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                seqid, label = line.strip().split()
                (pos if label == "1" else neg).add(seqid)
        return pos, neg
    except Exception as e:
        print(f"Error loading labels: {e}", file=sys.stderr)
        sys.exit(1)

def evaluate_performance(predicted: Set[str], positives: Set[str], negatives: Set[str]) -> Dict:
    """Calculate performance metrics
    
    Args:
        predicted: Set of predicted positive IDs
        positives: Set of true positive IDs
        negatives: Set of true negative IDs
        
    Returns:
        Dictionary containing performance metrics
    """
    tp = len(predicted & positives)
    fp = len(predicted & negatives)
    fn = len(positives - predicted)
    tn = len(negatives - predicted)
    
    metrics = {
        'TP': tp,
        'FP': fp,
        'FN': fn, 
        'TN': tn,
        'accuracy': (tp + tn) / (tp + tn + fp + fn) if (tp+tn+fp+fn) else 0,
        'precision': tp / (tp + fp) if (tp + fp) else 0,
        'recall': tp / (tp + fn) if (tp + fn) else 0
    }
    metrics['f1'] = 2 * (metrics['precision'] * metrics['recall']) / \
                   (metrics['precision'] + metrics['recall']) if (metrics['precision'] + metrics['recall']) else 0
    
    return metrics

def annotate_swissprot(hmm_file: Path, swissprot_fasta: Path, output_dir: Path) -> Path:
    """Run hmmsearch against SwissProt database
    
    Args:
        hmm_file: Path to HMM file
        swissprot_fasta: Path to SwissProt FASTA
        output_dir: Directory to save results
        
    Returns:
        Path to results table file
    """
    return run_hmmsearch(hmm_file, swissprot_fasta, output_dir, tag="swissprot")

def analyze_swissprot(tbl_file: Path) -> Set[str]:
    """Analyze SwissProt search results
    
    Args:
        tbl_file: Path to hmmsearch results
        
    Returns:
        Set of protein IDs found
    """
    hits = set()
    try:
        with open(tbl_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                prot = line.split()[0]
                hits.add(prot)
        print(f"Found {len(hits)} putative Kunitz domains in SwissProt")
        return hits
    except Exception as e:
        print(f"Error analyzing SwissProt results: {e}", file=sys.stderr)
        sys.exit(1)

def run_hmmlogo(hmm_file: Path, output_dir: Path) -> None:
    """Generate sequence logo from HMM
    
    Args:
        hmm_file: Path to HMM file
        output_dir: Directory to save logo
    """
    logo_file = output_dir / "hmm_logo.png"
    cmd = ["hmmlogo", "-o", str(logo_file), str(hmm_file)]
    
    try:
        print(f"Generating HMM logo: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        print(f"HMM logo saved to {logo_file}")
    except subprocess.CalledProcessError as e:
        print(f"Failed to generate HMM logo: {e}", file=sys.stderr)
        sys.exit(1)

# Load config globally
CONFIG = load_config()

# === MAIN PIPELINE === 
if __name__ == "__main__":
    # === Load config and define paths ===
    seed_alignment = Path(CONFIG["seed_alignment"])
    validation_fasta = Path(CONFIG["validation_fasta"])
    validation_labels = Path(CONFIG["validation_labels"])
    swissprot_fasta = Path(CONFIG["swissprot_fasta"])
    output_dir = CONFIG["output_dir"]

    # === Step 1: Build HMM from seed alignment ===
    hmm_file = output_dir / "kunitz.hmm"
    print(f"Building HMM from seed alignment: {seed_alignment}")
    try:
        subprocess.run(["hmmbuild", str(hmm_file), str(seed_alignment)], check=True)
    except subprocess.CalledProcessError as e:
        print(f"hmmbuild failed: {e}", file=sys.stderr)
        sys.exit(1)

    # === Step 2: Run hmmsearch on validation set ===
    val_tbl = run_hmmsearch(hmm_file, validation_fasta, output_dir, tag="validation")

    # === Step 3: Parse hits and evaluate ===
    predicted = parse_tblout(val_tbl)
    positives, negatives = load_labels(validation_labels)
    metrics = evaluate_performance(predicted, positives, negatives)

    print("\nValidation Performance:")
    for key, val in metrics.items():
        print(f"  {key}: {val:.3f}" if isinstance(val, float) else f"  {key}: {val}")

    # === Step 4: Annotate SwissProt ===
    swiss_tbl = annotate_swissprot(hmm_file, swissprot_fasta, output_dir)
    swiss_hits = analyze_swissprot(swiss_tbl)

    # === Step 5: Generate HMM logo ===
    run_hmmlogo(hmm_file, output_dir)

    print(f"\n✅ Pipeline finished. Results saved to: {output_dir}")
