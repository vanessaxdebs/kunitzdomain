# Kunitz Domain HMM Pipeline

> **Author:** vanessaxdebs  
> **Repository:** [github.com/vanessaxdebs/kunitzdomain](https://github.com/vanessaxdebs/kunitzdomain)
>
> **Description:**  
> This repository contains a reproducible pipeline to build, validate, and apply a profile Hidden Markov Model (HMM) for the Kunitz domain family using HMMER 3.x, as well as code for plotting and visualizing results such as confusion matrices.

---

## Quick Start

### 1. **Clone the Repository**

```sh
git clone https://github.com/vanessaxdebs/kunitzdomain.git
cd kunitzdomain
```

---

### 2. **Install Requirements**

#### HMMER 3.4+ (Required)

- [Download and install HMMER](http://hmmer.org/)
    - On Ubuntu/Debian:
    ```sh
    sudo apt-get install hmmer
    ```
    - Or build from source (see the [HMMER website](http://hmmer.org/)).

#### Python 3.7+ and packages

Create a virtual environment (optional but recommended):

```sh
python3 -m venv .venv
source .venv/bin/activate
```

Install required Python packages:

```sh
pip install matplotlib seaborn
```

#### (Optional) WebLogo for sequence logos

```sh
pip install weblogo
```

---

### 3. **Prepare Input Files**

Place your files in the `data/` directory:
- `data/kunitz_seed.sto` — Stockholm alignment of your Kunitz domain family (seed)
- `data/validation.fasta` — FASTA for validation (with known positive/negative examples)
- `data/validation_labels.txt` — True class labels for your validation set (for confusion matrix plotting)
- `data/swissprot.fasta` — SwissProt FASTA (can be downloaded, see below)

#### Download SwissProt FASTA (if needed)

```sh
mkdir -p data
curl -L "https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=reviewed:true" -o data/swissprot.fasta
```
**Note:**  
The SwissProt FASTA is very large! **Do NOT add it to git.** (see `.gitignore`)

---

### 4. **Run the Pipeline**

The main script is `hmm.py`. It will:
- Build an HMM profile from your seed alignment
- Validate against your validation set
- Search the SwissProt database
- Generate logo data and consensus sequences

```sh
python3 hmm.py
```

Results will appear in a timestamped subdirectory of `results/`, e.g. `results/run_20250519_1906/`.

---

### 5. **Plot a Confusion Matrix**

You can visualize your model’s performance on the validation set using `plot_confusion_matrix.py`.

**How to Use:**

1. Edit `plot_confusion_matrix.py` with your actual TP, FP, TN, FN values, or let it read them from your results.
2. Run:
    ```sh
    python3 plot_confusion_matrix.py
    ```
3. It will create `confusion_matrix_kunitz.png` like the one in 


---

### 6. **Create a Graphical Logo (Optional)**

You can create a PNG/SVG logo using the consensus FASTA output by `hmmemit`:

```sh
weblogo -F png < results/run_YYYYMMDD_HHMM/hmmemit.fasta > results/run_YYYYMMDD_HHMM/hmm_logo.png
```
Or, for SVG:
```sh
weblogo -F svg < results/run_YYYYMMDD_HHMM/hmmemit.fasta > results/run_YYYYMMDD_HHMM/hmm_logo.svg
```

Or use the [WebLogo web interface](https://weblogo.berkeley.edu/logo.cgi).

---

## Project Structure

```
.
├── assets/
│   ├── confusion_matrix_kunitz.png
│   └── pipeline.png
├── data/
│   ├── kunitz_seed.sto
│   ├── validation.fasta
│   ├── validation_labels.txt
│   └── swissprot.fasta
├── results/
│   └── run_YYYYMMDD_HHMM/
│         ├── kunitz.hmm
│         ├── hmmsearch_validation.tbl
│         ├── hmmsearch_swissprot.tbl
│         ├── hmm_logo.txt
│         ├── hmmemit.fasta
│         └── hmm_logo.png (optional)
├── hmm.py
├── plot_confusion_matrix.py
├── check_data.py
├── .gitignore
├── README.md
└── ... etc
```

---

##  Troubleshooting & FAQ

- **File too large error on push:**  
  Do not commit large FASTA files to git. Download them as needed and add to `.gitignore`.

- **hmmlogo/hmmemit not found:**  
  Ensure HMMER 3.4+ is installed and in your `$PATH`.

- **SwissProt download too large?**  
  Consider using a smaller FASTA subset for testing.

- **How do I get TP/FP/FN/TN for confusion matrix?**  
  Parse your `validation.tbl` and match IDs with `validation_labels.txt`.  
  (See `plot_confusion_matrix.py` for examples or ask for more help!)

---

##  Credits

- [HMMER](http://hmmer.org/)
- [WebLogo](https://weblogo.berkeley.edu/)
- [UniProt](https://www.uniprot.org/)

---

