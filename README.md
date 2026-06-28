# Searching-Responsible-Gene-for-Colour-Blindness

# Burrows–Wheeler Transform (BWT) Read Alignment and Exon Quantification

## Overview

This project implements a **Burrows–Wheeler Transform (BWT)** based DNA read aligner from first principles. The implementation performs efficient mapping of sequencing reads to the human **Chromosome X** reference genome using an FM-index inspired approach built with milestone tables, rank queries, and backward search.

Following alignment, the mapped reads are assigned to predefined exon regions of the **red** and **green opsin genes**, enabling exon-level read quantification and estimation of gene configuration probabilities.

Unlike conventional genomic analysis pipelines, this project avoids specialized alignment libraries and demonstrates the underlying algorithms behind modern short-read mappers.

---

## Objectives

* Load the BWT last-column representation of the reference genome.
* Construct milestone tables for efficient rank queries.
* Implement Rank and Select operations.
* Perform backward search using the FM-index.
* Locate sequencing reads within the reference genome.
* Support reverse-complement read alignment.
* Count reads overlapping red and green gene exons.
* Estimate gene configuration probabilities from exon read counts.

---

## Methodology

### Step 1: Data Loading

The program loads

* Chromosome X reference sequence
* Burrows–Wheeler last column
* Suffix array mapping
* Sequencing reads

---

### Step 2: FM-index Construction

A milestone table is constructed every **δ** characters.

The milestone matrix stores cumulative nucleotide counts

```text
Milestone(A)
Milestone(C)
Milestone(G)
Milestone(T)
```

allowing efficient rank computation without scanning the entire genome.

---

### Step 3: Rank Query

For a nucleotide **c**, the rank operation computes

```text
Rank(c,i)
=
Number of occurrences of c
from position 0 to i
```

using

```text
Rank(c,i)
=
Milestone(c)
+
Local Count
```

This reduces the computational cost from **O(n)** to approximately **O(δ)**.

---

### Step 4: Select Query

The select operation determines the genomic position corresponding to the **k-th** occurrence of a nucleotide.

Combined with Rank, this enables traversal of the FM-index.

---

### Step 5: Backward Search

Reads are aligned by processing nucleotides from right to left.

For every nucleotide

```text
Lower' = Select(c, Rank(c, Lower))

Upper' = Select(c, Rank(c, Upper))
```

The search interval is updated until

* a unique match is found,
* multiple candidate matches remain,
* or no alignment exists.

---

### Step 6: Reverse Complement Search

If a read cannot be aligned,

its reverse complement is generated

```text
A ↔ T

C ↔ G
```

and the alignment procedure is repeated.

---

### Step 7: Exon Quantification

Mapped read positions are compared against exon coordinates of

* Red opsin gene
* Green opsin gene

The number of reads overlapping each exon is computed to estimate exon-specific expression.

---

### Step 8: Configuration Probability

Given exon read counts

```text
Red Count = R

Green Count = G
```

the probability contribution of an exon is

```text
P =
Configuration ×
R / (R + G)
```

The final probability for a gene configuration is obtained by multiplying the exon-wise probabilities.

---

## Features

* Burrows–Wheeler Transform based read alignment
* FM-index implementation from scratch
* Rank and Select queries
* Milestone table optimization
* Reverse-complement alignment
* Partial read matching
* Exon-level read counting
* Gene configuration probability estimation

---

## Technologies Used

* Python
* NumPy
* Pandas
* Matplotlib
* tqdm
* SymPy

---

## Project Structure

```text
.
├── data/
│   ├── chrX.fa
│   ├── chrX_last_col.txt
│   ├── chrX_map.txt
│   └── reads
│
├── code.py
├── README.md
│
└── results/
```

---

## Installation

Clone the repository

```bash
git clone https://github.com/yourusername/bwt-read-aligner.git

cd bwt-read-aligner
```

Install dependencies

```bash
pip install numpy pandas matplotlib tqdm sympy
```

---

## Running the Project

Update the file paths inside `code.py` and execute

```bash
python code.py
```

The program will

* load the reference genome,
* construct milestone tables,
* perform FM-index based read alignment,
* align reverse complements when required,
* compute exon read counts,
* estimate gene configuration probabilities.

---

## Output

The implementation produces

* Read alignment positions
* Red gene exon counts
* Green gene exon counts
* Gene configuration probabilities
* Alignment statistics

---

## Computational Complexity

| Operation              | Complexity   |
| ---------------------- | ------------ |
| Milestone Construction | O(n)         |
| Rank Query             | O(δ)         |
| Select Query           | O(1)         |
| Backward Search        | O(m × δ)     |
| Read Alignment         | O(k × m × δ) |

where

* **n** = genome length
* **m** = read length
* **k** = number of reads
* **δ** = milestone interval

---

## Future Improvements

* Full FM-index implementation with wavelet trees.
* Approximate alignment supporting mismatches.
* Smith–Waterman refinement after seed matching.
* Paired-end read alignment.
* Parallel processing for large sequencing datasets.
* Support for compressed FASTA and FASTQ files.

---

## References

1. Burrows, M., & Wheeler, D. J. (1994). *A Block-Sorting Lossless Data Compression Algorithm.*
2. Ferragina, P., & Manzini, G. (2000). *Opportunistic Data Structures with Applications.*
3. Li, H., & Durbin, R. (2009). *Fast and Accurate Short Read Alignment with Burrows–Wheeler Transform.*

---

## License

This project is released under the MIT License.

---

## Author

Developed as part of a bioinformatics project to implement a Burrows–Wheeler Transform–based DNA read aligner and exon quantification pipeline from first principles using scientific computing in Python.
