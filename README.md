# Contig Builder
## Overview
Python project for building a larger contiguous sequence from an initial seed query sequence and set of reads.
The current design builds the sequence from kmer lookup tables. Currently, the repo is in its experimental stage and contained within a jupyter notebook used to test the lookup table behavior.

## Current Status
Completed so far:

FASTA loading utilities
Prefix and suffix helpers with window shift arg
Reverse complement generation
End kmer lookup table construction
Kmer length experiments
Window depth experiments

In progress:

best_seq design and testing
Contig construction logic
Relative location tracking
Final contig coordinate conversion

## Installation
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt

## Notes
Repo is still being organized. Current implementation is exploratory