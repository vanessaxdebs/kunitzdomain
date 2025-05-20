# Kunitz Domain Profile HMM Pipeline

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Overview

This repository provides a **bioinformatics pipeline** for the identification of Kunitz-type protease inhibitor domains using a structure-informed profile Hidden Markov Model (HMM). The pipeline guides users from structural data collection and alignment through HMM construction, validation, and domain discovery/annotation in SwissProt.  
The project is inspired by best practices in protein domain annotation and is designed for reproducibility and clarity.

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Summary](#pipeline-summary)
- [Directory Structure](#directory-structure)
- [Setup and Installation](#setup-and-installation)
- [Input Data Requirements](#input-data-requirements)
- [Pipeline Steps](#pipeline-steps)
- [Example Usage](#example-usage)
- [Expected Outputs](#expected-outputs)
- [Troubleshooting](#troubleshooting)
- [References](#references)
- [License](#license)

---

## Pipeline Summary

1. **Structural Data Collection:** Retrieve Kunitz-domain protein structures from the PDB and UniProt.
2. **Multiple Alignment:** Align representative sequences/structures and create a seed alignment (Stockholm format).
3. **Profile HMM Construction:** Build a profile HMM using HMMER.
4. **Model Validation:** Test the HMM on a labeled validation set, compute performance metrics.
5. **SwissProt Annotation:** Annotate SwissProt proteins using the trained HMM.
6. **Visualization & Reporting:** Generate sequence logos and performance plots; export results and a final report.

---

## Directory Structure

```text
.
├── config/
│   └── config.yaml             # (Optional) Pipeline configuration file
├── data/
│   ├── kunitz_seed.sto         # Seed alignment (Stockholm format)
│   ├── swissprot.fasta         # SwissProt sequences (FASTA)
│   ├── validation.fasta        # Validation set sequences (FASTA)
│   └── validation_labels.txt   # Validation set labels (tab-separated)
├── results/
│   ├── kunitz.hmm              # Trained profile HMM
│   ├── hmmsearch_validation.tbl # Validation search results
│   ├── hmmsearch_swissprot.tbl # SwissProt search results
│   ├── hmm_logo.png            # Sequence logo (if generated)
│   └── ...                     # Other result files and plots
├── hmm.py                      # Main pipeline script
├── check_data.py               # (Optional) Data integrity checker
├── plot_confusion_matrix.py    # (Optional) Script for confusion matrix visualization
├── environment.yaml            # Conda environment specification
├── README.md                   # This file
└── LICENSE
```

---

## Setup and Installation
## SwissProt Data Download

Due to its large size, the SwissProt FASTA file (`data/swissprot.fasta`) is **not included** in this repository.

**To download the SwissProt FASTA:**

1. Visit the [UniProt Downloads page](https://www.uniprot.org/downloads).
2. Download the latest **SwissProt (reviewed) FASTA** file.  
   - The file is typically named `uniprot_sprot.fasta.gz`.
3. Uncompress the file:
   ```bash
   gunzip uniprot_sprot.fasta.gz
   ```
4. Move or copy the decompressed file to `data/swissprot.fasta` in this repository:
   ```bash
   mv uniprot_sprot.fasta data/swissprot.fasta
   ```

Alternatively, you can use the command line:
```bash
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
mv uniprot_sprot.fasta data/swissprot.fasta
```

**Note:**  
- If you only want a subset of SwissProt (e.g., specific taxa or reviewed human proteins), use UniProt's [Advanced Search](https://www.uniprot.org/uniprotkb?query=reviewed:true) and download the filtered FASTA file.
- Make sure the file is named exactly `swissprot.fasta` and placed in the `data/` directory.

---
### 1. Clone the repository

```bash
git clone https://github.com/vanessaxdebs/kunitzdomain.git
cd kunitzdomain
```

### 2. Create the environment

We recommend using **conda** for reproducibility.  
Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) if needed.

```bash
conda env create -f environment.yaml
conda activate kunitz
```

### 3. Install External Tools

- **HMMER** (required):  
  Download from [http://hmmer.org/](http://hmmer.org/) and add to your `PATH`.
- **MUSTANG** (for structural alignment, optional for re-creating seed):  
  [https://github.com/rdk/pymustang](https://github.com/rdk/pymustang)
- **JalView** (optional, for alignment editing):  
  [https://www.jalview.org/](https://www.jalview.org/)

---

## Input Data Requirements

- `data/kunitz_seed.sto`:  
  Multiple sequence alignment of Kunitz domains in **Stockholm format**.
- `data/validation.fasta` and `data/validation_labels.txt`:  
  Validation sequences and label file (tab-separated: `seqid [tab] 1/0`).
- `data/swissprot.fasta`:  
  SwissProt protein sequences in FASTA format.

If these files are not present, see the Methods section in the [report](./report.md) for instructions on how to generate them.

---

## Pipeline Steps

The main pipeline is run via the `hmm.py` script.

### **Running the Pipeline**

```bash
python hmm.py
```

This script:

1. Checks input files.
2. Builds the profile HMM using `hmmbuild`.
3. Validates the HMM with `hmmsearch`, computes metrics (accuracy, precision, recall, F1).
4. Searches SwissProt for Kunitz domains.
5. Attempts to generate a sequence logo (`hmmlogo`).
6. Outputs results and summary to the `results/` directory.

**Note:** Some steps (like logo generation) require additional tools (`hmmlogo`, etc.).  
If a step fails, the script will print a warning, but main results will still be produced.

---

## Example Usage

### **Basic run**

```bash
python hmm.py
```

### **Custom configuration**

If using a custom YAML config file:

```bash
python hmm.py --config config/config.yaml
```

---

## Expected Outputs

- `results/kunitz.hmm`: The profile HMM.
- `results/hmmsearch_validation.tbl`: Table of validation search results.
- `results/hmmsearch_swissprot.tbl`: Table of SwissProt search results.
- `results/confusion_matrix.png`: Confusion matrix plot of validation results.
- `results/hmm_logo.png`: Sequence logo (if generated).
- `results/final_report.md`: Pipeline report in markdown format (if script supports).
- Other intermediate files and visualizations.

---

## Troubleshooting

- **Some files/folders not appearing in GitHub:**  
  Make sure you have committed and pushed all files. If `results/` is in your `.gitignore`, remove or comment out that line and push again.
- **HMMER not found:**  
  Ensure `hmmbuild` and `hmmsearch` are installed and in your `PATH`.
- **hmmlogo fails:**  
  This tool is not always bundled with HMMER. If unavailable, you can visualize the domain logo using [Pfam](https://pfam.xfam.org/family/PF00014#tabview=tab5).
- **Validation set metrics are perfect (e.g., accuracy=1.0):**  
  This may indicate a small or non-diverse test set. Try using a more challenging or larger validation set if needed.

---

## License

This project is licensed under the [MIT License](./LICENSE).


