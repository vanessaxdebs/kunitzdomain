#!/usr/bin/env python3
"""
Kunitz-type Protease Inhibitor Domain - HMM Profile Pipeline (Biopython version)
Compatible with hmmologs project standards.
Includes positive and negative validation.
"""

import os
import sys
import subprocess
from pathlib import Path
from typing import Dict, Set, Tuple
import yaml
from datetime import datetime
from Bio import SeqIO

# ---------- Biopython-based utilities ----------

def sto_to_ungapped_fasta(sto_path, fasta_path):
    seqs = {}
    with open(sto_path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("//") or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) == 2:
                name, seq = parts
                seqs.setdefault(name, "")
                seqs[name] += seq
    with open(fasta_path, "w") as outfa:
        for name, seq in seqs.items():
            seq = seq.replace("-", "").replace(".", "").replace(" ", "")
            outfa.write(f">{name}\n{seq}\n")
    print(f"Converted {len(seqs)} sequences from {sto_path} to {fasta_path}")

def sto_to_fasta(sto_path, fasta_path):
    count = SeqIO.write(SeqIO.parse(sto_path, "stockholm"), fasta_path, "fasta")
    print(f"Converted {count} sequences from {sto_path} to {fasta_path}")

def fasta_to_label_txt(fasta_path, txt_path, label="1"):
    count = 0
    with open(txt_path, "w") as txt:
        for record in SeqIO.parse(fasta_path, "fasta"):
            txt.write(f"{record.id}\t{label}\n")
            count += 1
    print(f"Wrote {count} sequence labels to {txt_path}")

def check_stockholm(path):
    try:
        with open(path) as f:
            next(SeqIO.parse(f, "stockholm"))
        return True
    except Exception:
        print(f"ERROR: {path} is not a valid Stockholm file.")
        return False

def check_fasta(path):
    try:
        with open(path) as f:
            next(SeqIO.parse(f, "fasta"))
        return True
    except Exception:
        print(f"ERROR: {path} is not a valid FASTA file.")
        return False

def check_label_txt(path):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        print(f"ERROR: {path} is missing or empty.")
        return False
    with open(path) as f:
        for line in f:
            if line.strip() and len(line.strip().split()) != 2:
                print(f"ERROR: {path} has a line with wrong format: '{line.strip()}'")
                return False
    return True

# ---------- Config loading & output dir creation ----------

def find_config() -> Path:
    script_dir = Path(__file__).parent
    for path in [
        script_dir / "../config/config.yaml",
        script_dir / "config/config.yaml",
        Path("config/config.yaml")
    ]:
        if path.exists():
            return path.resolve()
    raise FileNotFoundError("config.yaml not found!")

def load_config(config_file: Path = None) -> Dict:
    if config_file is None:
        config_file = find_config()
    required_fields = {
        'output_dir': str,
        'seed_alignment': str,
        'validation_fasta': str,
        'validation_labels': str,
        'negative_fasta': str,
        'negative_labels': str,
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
        run_id = datetime.now().strftime("%Y%m%d_%H%M")
        output_dir = Path(config['output_dir']) / f"run_{run_id}"
        output_dir.mkdir(parents=True, exist_ok=True)
        config['output_dir'] = output_dir
        return config
    except Exception as e:
        print(f"Configuration error: {e}", file=sys.stderr)
        sys.exit(1)

# ---------- HMM and metrics functions ----------

def run_hmmbuild(seed_alignment: Path, hmm_file: Path):
    cmd = ["hmmbuild", str(hmm_file), str(seed_alignment)]
    print(f"Building HMM: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"hmmbuild failed: {e}", file=sys.stderr)
        sys.exit(1)

def run_hmmsearch(hmm_file: Path, fasta_file: Path, output_dir: Path, tag: str = "validation", e_value: float = 1e-5) -> Path:
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
        print(f"hmmsearch failed. Error:\n{e.stderr.decode()}", file=sys.stderr)
        sys.exit(1)

def parse_tblout(tbl_file: Path) -> Set[str]:
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
    logo_file = output_dir / "hmm_logo.png"
    cmd = ["hmmlogo", "-o", str(logo_file), str(hmm_file)]
    print(f"Generating HMM logo: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
        print(f"HMM logo saved to {logo_file}")
    except subprocess.CalledProcessError as e:
        print(f"Failed to generate HMM logo: {e}", file=sys.stderr)
        sys.exit(1)

# ---------- Main pipeline ----------

if __name__ == "__main__":
    CONFIG = load_config()

    sto_file = Path(CONFIG["seed_alignment"])
    validation_fasta = Path(CONFIG["validation_fasta"])
    validation_labels = Path(CONFIG["validation_labels"])
    negative_fasta = Path(CONFIG["negative_fasta"])
    negative_labels = Path(CONFIG["negative_labels"])
    swissprot_fasta = Path(CONFIG["swissprot_fasta"])
    output_dir = CONFIG["output_dir"]
    e_value = CONFIG["e_value_cutoff"]

    # Convert seed alignment Stockholm to FASTA if needed
    fasta_seed = output_dir / "seed.fasta"
    if not fasta_seed.exists() or fasta_seed.stat().st_size == 0:
        sto_to_fasta(sto_file, fasta_seed)

    # Convert validation Stockholm to FASTA if needed
    if not validation_fasta.exists() or validation_fasta.stat().st_size == 0:
        print("Validation FASTA missing or empty. Please provide it.")
        sys.exit(1)

    # Check files
    if not check_stockholm(sto_file):
        sys.exit(1)
    if not check_fasta(validation_fasta):
        sys.exit(1)
    if not check_fasta(negative_fasta):
        sys.exit(1)
    if not check_label_txt(validation_labels):
        sys.exit(1)
    if not check_label_txt(negative_labels):
        sys.exit(1)

    # Build HMM
    hmm_file = output_dir / "kunitz.hmm"
    run_hmmbuild(sto_file, hmm_file)

    # Run hmmsearch on validation positives
    val_tbl = run_hmmsearch(hmm_file, validation_fasta, output_dir, tag="validation", e_value=e_value)
    val_predicted = parse_tblout(val_tbl)
    val_positives, val_negatives = load_labels(validation_labels)

    # Run hmmsearch on negative sequences
    neg_tbl = run_hmmsearch(hmm_file, negative_fasta, output_dir, tag="negative", e_value=e_value)
    neg_predicted = parse_tblout(neg_tbl)
    neg_positives, neg_negatives = load_labels(negative_labels)

    # Combine predictions and labels for evaluation
    all_predicted = val_predicted.union(neg_predicted)
    all_positives = val_positives  # Positives only from validation labels
    all_negatives = val_negatives.union(neg_negatives)  # Combine negative seq IDs

    # Evaluate combined performance
    metrics = evaluate_performance(all_predicted, all_positives, all_negatives)

    print("\nCombined Validation Performance (Positives + Negatives):")
    for key, val in metrics.items():
        print(f"  {key}: {val:.3f}" if isinstance(val, float) else f"  {key}: {val}")

    # Save false positives and false negatives for analysis
    fp = all_predicted & all_negatives
    fn = all_positives - all_predicted

    fp_path = output_dir / "false_positives.txt"
    fn_path = output_dir / "false_negatives.txt"

    with open(fp_path, "w") as f:
        for seq_id in sorted(fp):
            f.write(seq_id + "\n")

    with open(fn_path, "w") as f:
        for seq_id in sorted(fn):
            f.write(seq_id + "\n")

    print(f"False Positives: {len(fp)} saved to {fp_path.name}")
    print(f"False Negatives: {len(fn)} saved to {fn_path.name}")

    # Annotate SwissProt
    swiss_tbl = run_hmmsearch(hmm_file, swissprot_fasta, output_dir, tag="swissprot", e_value=e_value)
    swiss_hits = parse_tblout(swiss_tbl)
    print(f"\nSwissProt hits: {len(swiss_hits)}")

    # Generate HMM logo
    run_hmmlogo(hmm_file, output_dir)

    print(f"\n✅ Pipeline finished. Results saved to: {output_dir}")
