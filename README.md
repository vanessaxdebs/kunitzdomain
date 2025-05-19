## Large File Notice

This repository includes a large data file: `data/swissprot.fasta` (~273 MB).

GitHub limits file sizes to 100 MB. To handle this file properly, it is managed using **Git Large File Storage (Git LFS)**.

### How to work with this large file

- Install Git LFS on your machine:  
  https://git-lfs.github.com/

- After cloning this repo, run:

```bash
git lfs install
git lfs pull
```

This will download the large file efficiently without hitting GitHub’s file size limits.

