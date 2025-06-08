# Kunitz Domain Profile HMM Pipeline

**Author:** Vanessa debs  
**Affiliation:** FABIT, University of Bologna

---

## Overview

This repository provides a robust, reproducible bioinformatics pipeline for the **identification and annotation of Kunitz-type protease inhibitor domains** using a structure-guided profile Hidden Markov Model (HMM). The workflow guides users through structural data acquisition, multiple alignment, HMM construction, rigorous validation, and large-scale domain annotation in SwissProt. All steps are documented and implemented in Python, with results ready for direct interpretation or further analysis.

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Summary](#pipeline-summary)
- [Directory Structure](#directory-structure)
- [Installation and Setup](#installation-and-setup)
- [Data Acquisition](#data-acquisition)
- [Pipeline Execution](#pipeline-execution)
- [Result Visualization](#result-visualization)
- [SwissProt FASTA Download](#swissprot-fasta-download)
- [Novelty Check](#novelty-check)
- [Troubleshooting](#troubleshooting)
- [References](#references)
- [License](#license)
- [Report Download](#report-download)

---

## Pipeline Summary

The pipeline comprises the following major steps:

1. **Structural Data Retrieval:** Obtain and curate high-quality Kunitz domain structures from UniProt and the PDB.
2. **Multiple Structural Alignment:** Align domain instances using structural alignment tools to generate a representative seed alignment.
3. **Profile HMM Construction:** Build a profile HMM from the alignment using the HMMER suite.
4. **Validation:** Evaluate the model with an independent, labeled test set and compute standard performance metrics.
5. **SwissProt Annotation:** Apply the HMM to SwissProt protein sequences to identify and annotate Kunitz domains on a proteome-wide scale.
6. **Visualization and Reporting:** Generate plots (e.g., confusion matrix) and export all relevant outputs for publication or further study.

---

## Directory Structure

```text
.
├── config/
│   └── config.yaml               # (Optional) Pipeline configuration
├── data/
│   ├── kunitz_seed.sto           # Seed alignment (Stockholm format)
│   ├── swissprot.fasta           # SwissProt sequences (user-provided, see below)
│   ├── validation.fasta          # Validation set sequences (FASTA)
│   └── validation_labels.txt     # Validation set labels
├── results/
│   ├── kunitz.hmm                # Trained profile HMM
│   ├── hmmsearch_validation.tbl  # Validation set search results
│   ├── hmmsearch_swissprot.tbl   # SwissProt search results
│   ├── confusion_matrix.png      # Confusion matrix plot
│   ├── hmm_logo.png              # Sequence logo (if generated)
│   └── ...                       # Additional outputs
├── hmm.py                        # Main pipeline script
├── novelty_check.py              # Novelty candidate identification script
├── check_data.py                 # Data integrity checker
├── plot_confusion_matrix.py      # Confusion matrix plotting script
├── environment.yaml              # Conda environment specification
├── README.md                     # This file
├── LICENSE
└── report.zip                    # PDF report and supplementary materials (see below)
```

---

## Installation and Setup

### 1. Clone the Repository

```bash
git clone https://github.com/vanessaxdebs/kunitzdomain.git
cd kunitzdomain
```

### 2. Create and Activate the Environment

We recommend using [conda](https://docs.conda.io/en/latest/miniconda.html):

```bash
conda env create -f environment.yaml
conda activate kunitz
```

### 3. Install External Dependencies

- **HMMER** (required): [http://hmmer.org/](http://hmmer.org/)
- **MUSTANG** (optional, for seed alignment): [https://github.com/rdk/pymustang](https://github.com/rdk/pymustang)
- **JalView** (optional, for alignment editing): [https://www.jalview.org/](https://www.jalview.org/)

Ensure `hmmbuild` and `hmmsearch` are on your `PATH`.

---

## Data Acquisition

- **Seed Alignment:**  
  `data/kunitz_seed.sto`  
  Provided in the repository (or see report for methods to generate).

- **Validation Set:**  
  `data/validation.fasta` and `data/validation_labels.txt`  
  Provided in the repository.

- **SwissProt FASTA:**  
  **Not included** due to file size. See [SwissProt FASTA Download](#swissprot-fasta-download) below.

---

## Pipeline Execution

### 1. Run the Main Pipeline

```bash
python hmm.py
```

This script will:

- Build the profile HMM from the seed alignment.
- Validate the model using the provided validation set.
- Search for Kunitz domains in SwissProt (if `data/swissprot.fasta` is present).
- Attempt to generate a sequence logo (if `hmmlogo` is available).
- Save all outputs to the `results/` directory.

### 2. Plot the Confusion Matrix

After running the main pipeline, visualize model performance:

```bash
python plot_confusion_matrix.py
```

This will generate `results/confusion_matrix.png`, summarizing true/false positives and negatives.

### 3. (Optional) Data Integrity Check

If you wish to verify your input data:

```bash
python check_data.py
```

---

## Result Visualization

- **Confusion Matrix:**  
  Visualizes model performance on the validation set.
- **Sequence Logo:**  
  If the logo generation fails, consult [Pfam Kunitz logo](https://pfam.xfam.org/family/PF00014#tabview=tab5) for a reference image.
- **Tabular Results:**  
  `results/hmmsearch_validation.tbl` and `results/hmmsearch_swissprot.tbl` contain all matches and scores.

---

## SwissProt FASTA Download

Due to size constraints, `data/swissprot.fasta` must be downloaded by the user:

1. **From UniProt:**
   - Visit the [UniProt Downloads page](https://www.uniprot.org/downloads).
   - Download the latest **SwissProt (reviewed) FASTA** file (typically `uniprot_sprot.fasta.gz`).
   - Uncompress and move it:
     ```bash
     gunzip uniprot_sprot.fasta.gz
     mv uniprot_sprot.fasta data/swissprot.fasta
     ```

2. **Command Line Alternative:**
   ```bash
   wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
   gunzip uniprot_sprot.fasta.gz
   mv uniprot_sprot.fasta data/swissprot.fasta
   ```

3. **Subset Download:**  
   Use UniProt's [Advanced Search](https://www.uniprot.org/uniprotkb?query=reviewed:true) to download only a portion of SwissProt (e.g., for a specific organism).

**Ensure the file is named `swissprot.fasta` and is placed in the `data/` directory.**

---

## Novelty Check

After running the main pipeline and generating SwissProt annotation results, you can identify potentially novel Kunitz domain candidates by running the novelty check script:

```bash
python novelty_check.py
```

This script will:

- Parse the latest `results/run_*/hmmsearch_swissprot.tbl` output.
- Extract UniProt accessions from the HMM search results.
- Query UniProt for each hit to check if it is already annotated as Kunitz by keyword, Pfam domain, or feature.
- List all hits **not** annotated as Kunitz in `results/novel_kunitz_candidates.txt`.

Use this list to prioritize candidates for further manual inspection or experimental validation.

---

## Troubleshooting

- **File Not Found Errors:**  
  Verify that all required input files are present and correctly named.
- **SwissProt Search Fails:**  
  Ensure you have downloaded and placed the SwissProt FASTA as described above.
- **HMMER Not Found:**  
  Confirm that the HMMER executables are installed and available in your environment.
- **hmmlogo Issues:**  
  If logo generation fails, use the official Pfam logo as a substitute for reports.
- **Plotting Errors:**  
  Ensure all plotting dependencies (`matplotlib`, `seaborn`) are installed via conda.

---

## References

- Armstrong, D. R., et al. (2019). PDBe: improved findability of macromolecular structure data in the PDB. _Nucleic Acids Research_.
- Blum, M., et al. (2024). InterPro: the protein sequence classification resource in 2025. _Nucleic Acids Research_.
- Eddy, S. R. (2011). Accelerated profile HMM searches. _PLoS Computational Biology_.
- Konagurthu, A. S., et al. (2006). MUSTANG: A multiple structural alignment algorithm. _Proteins: Structure, Function, and Bioinformatics_.
- Krogh, A., et al. (1994). Hidden Markov models in computational biology. _Journal of Molecular Biology_.

---

## License

This project is licensed under the [MIT License](./LICENSE).

---

## Report Download

A detailed PDF report and supplementary materials for this project are available for download in Structure-Informed Profile HMMs for Kunitz Domain Detection_ A Reproducible Pipeline and Its Benchmarking (1).zip

---

**For questions, issues, or contributions, please use the GitHub issue tracker or contact the repository owner.** 
