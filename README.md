# Contig Builder

## Overview

Contig Builder is a Python project for constructing a larger contiguous sequence from an initial query sequence and a set of sequencing reads. The program starts from a known query sequence, identifies candidate overlapping reads using kmer lookup tables, scores candidate reads by overlap quality and sequence added, and greedily extends the query into a final contig.

The project was developed for a biomedical informatics programming assignment focused on building the largest possible contig containing an initial query sequence.

## Method Summary

The assembler uses a greedy query seeded overlap extension strategy.

The main steps are:

1. Load the query sequence and sequencing reads from FASTA files.
2. Build kmer lookup tables from the reads.
3. Use a prefix table to extend the query to the right.
4. Use a suffix table to extend the resulting contig to the left.
5. Score candidate reads based on overlap length, mismatch rate, and added sequence.
6. Track read placement coordinates during contig construction.
7. Convert internal coordinates into the final alignment format.
8. Write the final contig to `ALLELES.fasta` and read placements to `ALLELES.aln`.

The final assembler removes reads from future consideration after they are incorporated into the contig. This prevents the same read identifier from being used multiple times in the final assembly.

## Repository Structure

```text
FinalRepo/
    README.md
    requirements.txt

    data/
        QUERY.fasta
        READS.fasta

    outputs/
        ALLELES.fasta
        ALLELES.aln

    src/
        contigbuilder/
            __init__.py
            io.py
            sequenceUtils.py
            kmer_index.py
            scoring.py
            assembler.py
            alignment.py
            cli.py

    tests/
        test_sequence_utils.py
        test_alignment.py
        test_scoring.py
        test_kmer_index.py
        test_assembler.py
        test_io.py
        test_end_to_end.py

    notebooks/
        ExperimentalContigBuilder.ipynb

## Installation

Clone the repository and move into the project directory.

Create and activate a virtual environment.

```bash
python -m venv venv
source venv/bin/activate
```

Install dependencies.

```bash
pip install -r requirements.txt
```

## Input Files

The program expects two FASTA files.

### Query FASTA

`QUERY.fasta` contains the initial seed query sequence.

Example:

```text
>query1
ACGTCATGCA
```

### Reads FASTA

`READS.fasta` contains sequencing reads.

Example:

```text
>read1
AAAACGTCAT
>read2
TCATGCATTT
```

The current FASTA parser assumes each sequence is stored on a single line.

## Usage

Run the assembler from the repository root using:

```bash
PYTHONPATH=src python -m contigbuilder.cli \
  --query data/QUERY.fasta \
  --reads data/READS.fasta \
  --out_fasta outputs/ALLELES.fasta \
  --out_aln outputs/ALLELES.aln \
  --k_right 7 \
  --n_right 4 \
  --k_left 7 \
  --n_left 4 \
  --mRate 0.1
```

## Parameters

| Argument | Description | Default |
|---|---|---|
| `--query` | Path to query FASTA file | Required |
| `--reads` | Path to reads FASTA file | Required |
| `--out_fasta` | Output path for final contig FASTA | `ALLELES.fasta` |
| `--out_aln` | Output path for alignment table | `ALLELES.aln` |
| `--k_right` | Kmer length used for right extension | `7` |
| `--n_right` | Maximum read window used for right extension | `6` |
| `--k_left` | Kmer length used for left extension | `7` |
| `--n_left` | Maximum read window used for left extension | `5` |
| `--mRate` | Maximum allowed mismatch rate in the overlap | `0.1` |
| `--maxsteps` | Maximum number of extension steps per direction | `10000` |

## Output Files

### `ALLELES.fasta`

FASTA file containing the final assembled contig.

Example:

```text
>contig1
AAAACGTCATGCATTT
```

### `ALLELES.aln`

Tab separated alignment file describing which reads were incorporated into the final contig.

The output includes:

| Column | Description |
|---|---|
| `sseqid` | Read identifier |
| `qseqid` | Contig identifier |
| `sstart` | Start coordinate in the sequencing read |
| `send` | End coordinate in the sequencing read |
| `qstart` | Start coordinate in the final contig |
| `qend` | End coordinate in the final contig |
| `strand` | Forward or reverse complement orientation |
| `overlap_len` | Length of validated overlap |
| `score` | Candidate score used during extension |

Forward reads are reported with `send > sstart`. Reverse complement reads are reported with `send < sstart`.

## Scoring Function

Candidate reads are scored using:

```text
score = overlap length * (1 minus mismatch rate) + 0.5 * added sequence length
```

Candidates are skipped if:

1. The mismatch rate is greater than the allowed threshold.
2. The candidate does not add new sequence.
3. The candidate read identifier has already been used.

This score rewards strong overlap, penalizes disagreement, and gives partial credit for reads that extend the contig.

## Testing

Unit tests are included in the `tests/` directory.

Run all tests from the repository root with:

```bash
PYTHONPATH=src pytest tests/
```

Run a specific test file with:

```bash
PYTHONPATH=src pytest tests/test_alignment.py
```

The current test suite covers:

1. Prefix and suffix sequence utilities.
2. Reverse complement generation.
3. Mismatch counting.
4. Kmer lookup table creation.
5. Candidate scoring.
6. Right and left contig extension.
7. Coordinate conversion for `ALLELES.aln`.
8. FASTA and alignment output writing.
9. A small end to end assembly example.

## Design Notes and Limitations

This project uses a greedy local extension strategy rather than a graph based assembler. At each step, the highest scoring candidate read is selected and added to the contig. This makes the method simple and computationally efficient, but it also makes the final result path dependent.

Important limitations include:

1. The final contig can depend on whether right or left extension is performed first.
2. Removing used reads can shorten the final contig but prevents the same read from being used multiple times.
3. Similar or duplicate reads may still be selected even after one read identifier is removed.
4. The scoring function was designed heuristically and was not fully optimized.
5. The current FASTA parser assumes one sequence line per FASTA record.
6. The method does not perform consensus correction across all reads aligned to the contig.

Future improvements could include testing both extension orders, bidirectional extension, branch tracking, consensus correction, stricter duplicate read handling, score optimization, or a graph based representation of alternative paths.

## Current Recommended Settings

Based on parameter experiments, the strict final assembler was run with:

```text
k_right = 7
n_right = 4
k_left = 7
n_left = 4
mRate = 0.1
```

Exploratory experiments also tested asymmetric right and left window settings. Those experiments were useful for understanding algorithm behavior, but the final strict assembler removes previously used reads across extension steps.
To obtain the longest contig, allowing the reuse of reads for left extension with the following settings:
```text
k_right = 7
n_right = 6
k_left = 7
n_left = 5
mRate = 0.1
```

## Author

Nicholas Cordaro