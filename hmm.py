#!/usr/bin/env python3
"""
Kunitz-type Protease Inhibitor Domain - HMM Profile Pipeline
Compatible with standard bioinformatics and hmmologs project practices.
"""

import os
import sys
import subprocess
from pathlib import Path
from typing import Dict, Set, Tuple
import yaml
from datetime import datetime

def check_file(path: str, format: str = "auto") -> bool:
    """Check for file existence, non-emptiness, and optional format."""
    if not os.path.isfile(path):
        print(f"ERROR: Required file '{path}' is missing.")
        return False
    if os.path.getsize(path) == 0:
        print(f"ERROR: Required file '{path}' is empty.")
        return False

    with open(path) as f:
        first_line = f.readline().strip()
        # Check Stockholm
        if format == "stockholm":
            if not first_line.startswith("# STOCKHOLM 1.0"):
                print(f"ERROR: '{path}' is not in Stockholm format.")
                return False
        # Check FASTA
        if format == "fasta":
            if not first_line.startswith(">"):
                print(f"ERROR: '{path}' does not start with '>' (not FASTA format).")
                return False
            for i, line in enumerate(f, 2):
                if line.startswith("#=GF") or line.startswith("# STOCKHOLM"):
                    print(f"ERROR: '{path}' contains Stockholm or annotation lines at line {i}. Only plain FASTA allowed.")
                    return False
    return True

def find_config() -> Path:
    """Locate config.yaml in the project."""
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
    """Load YAML config, check required fields, create output dir."""
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
        for field, field_type in required_fields.items():
            if field not in config:
                raise ValueError(f"Missing required config field: {field}")
            if not isinstance(config[field], field_type):
                # Allow int for float fields
                if field_type == float and isinstance(config[field], int):
                    config[field] = float(config[field])
                else:
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

def run_hmmbuild(seed_alignment: Path, hmm_file: Path):
    """Build HMM from Stockholm alignment."""
    cmd = ["hmmbuild", str(hmm_file), str(seed_alignment)]
    print(f"Building HMM: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"hmmbuild failed: {e}", file=sys.stderr)
        sys.exit(1)

def run_hmmsearch(hmm_file: Path, fasta_file: Path, output_dir: Path, tag: str = "validation", e_value: float = 1e-5) -> Path:
    """Run hmmsearch and save table output."""
    tblout = output_dir / f"hmmsearch_{tag}.tbl"
    cmd = [
        "hmmsearch",
        "--tblout", str(tblout),
        "-E", str(e_value),
        str(hmm_file), str(fasta_file)
    ]
    print(f"Running hmmsearch: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        return tblout
    except subprocess.CalledProcessError as e:
        print(f"hmmsearch failed. Error:\n{e.stderr}", file=sys.stderr)
        sys.exit(1)

def parse_tblout(tbl_file: Path) -> Set[str]:
    """Parse HMMER tblout file, extract sequence IDs."""
    hits = set()
    try:
        with open(tbl_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split()
                if parts:
                    hits.add(parts[0])
        return hits
    except Exception as e:
        print(f"Error parsing tblout file: {e}", file=sys.stderr)
        sys.exit(1)

def load_labels(label_file: Path) -> Tuple[Set[str], Set[str]]:
    """Load validation labels (tab-separated: seqid [tab] 1/0)."""
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
    """Calculate performance metrics."""
    tp = len(predicted & positives)
    fp = len(predicted & negatives)
    fn = len(positives - predicted)
    tn = len(negatives - predicted)
    total = tp + tn + fp + fn
    acc = (tp + tn) / total if total else 0
    prec = tp / (tp + fp) if (tp + fp) else 0
    rec = tp / (tp + fn) if (tp + fn) else 0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0
    return dict(TP=tp, FP=fp, FN=fn, TN=tn, accuracy=acc, precision=prec, recall=rec, f1=f1)

def run_hmmlogo(hmm_file: Path, output_dir: Path):
    """Generate sequence logo from HMM."""
    logo_file = output_dir / "hmm_logo.png"
    cmd = ["hmmlogo", "-o", str(logo_file), str(hmm_file)]
    print(f"Generating HMM logo: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
        print(f"HMM logo saved to {logo_file}")
    except subprocess.CalledProcessError as e:
        print(f"Failed to generate HMM logo: {e}", file=sys.stderr)
        sys.exit(1)

# === MAIN PIPELINE ===
if __name__ == "__main__":
    # --- File checks before config load ---
    if not check_file("data/kunitz_seed.sto", format="stockholm"):
        sys.exit(1)
    if not check_file("data/validation.fasta", format="fasta"):
        sys.exit(1)
    if not check_file("data/validation_labels.txt"):
        sys.exit(1)

    CONFIG = load_config()
    seed_alignment = Path(CONFIG["seed_alignment"])
    validation_fasta = Path(CONFIG["validation_fasta"])
    validation_labels = Path(CONFIG["validation_labels"])
    swissprot_fasta = Path(CONFIG["swissprot_fasta"])
    output_dir = CONFIG["output_dir"]
    e_value = CONFIG["e_value_cutoff"]

    # Step 1: Build HMM
    hmm_file = output_dir / "kunitz.hmm"
    run_hmmbuild(seed_alignment, hmm_file)

    # Step 2: Run hmmsearch on validation set
    val_tbl = run_hmmsearch(hmm_file, validation_fasta, output_dir, tag="validation", e_value=e_value)

    # Step 3: Parse hits and evaluate
    predicted = parse_tblout(val_tbl)
    positives, negatives = load_labels(validation_labels)
    metrics = evaluate_performance(predicted, positives, negatives)

    print("\nValidation Performance:")
    for key, val in metrics.items():
        print(f"  {key}: {val:.3f}" if isinstance(val, float) else f"  {key}: {val}")

    # Step 4: Annotate SwissProt
    swiss_tbl = run_hmmsearch(hmm_file, swissprot_fasta, output_dir, tag="swissprot", e_value=e_value)
    swiss_hits = parse_tblout(swiss_tbl)
    print(f"\nSwissProt hits: {len(swiss_hits)}")

    # Step 5: Generate HMM logo
    run_hmmlogo(hmm_file, output_dir)

    print(f"\n✅ Pipeline finished. Results saved to: {output_dir}")
