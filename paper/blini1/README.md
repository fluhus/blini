# Paper-related material

This directory contains publication-related files that are
not part of the core program.
This is mostly the manuscripts, plots and benchmarking code.

## Directory structure

**Subdirectories:**

* `gentestdata`: Generates the simulated viral dataset.
* `gentestdatabig`: Generates the simulated bacterial dataset.
* `simul`: Library for simulating data.
* `testclust`: Tests the clustering results.
* `testsearch`: Tests the search results.
* `results`: Results from running the comparisons.
  * `*.txt`: Manually filled files with results.
  * `*.png`: Figures created from the text files.

**Files:**

* `paper.md`: The JOSS manuscript text.
* `paper2.md`: The Arxiv manuscript text.
* `paper.bib`: Paper bibliography.
* `template.tex`: Pandoc template for `paper2.md`.
* `plots.py`: Generates plots for the paper using `results/*.txt`.
* `run_comparisons.sh`: Commands for running the benchmarked tools.

## Running benchmarks

### Required data & code

* [Viral dataset](https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/) (genomic.fna)
* [Bacterial dataset](http://segatalab.cibio.unitn.it/data/Pasolli_et_al.html) (4930 SGBs)
  * Need to extract the files and concatenate them into one file
    `representatives.fa`.
* Tools compared against:
  [mmseqs2](https://github.com/soedinglab/MMseqs2/) and
  [sourmash](https://github.com/sourmash-bio/sourmash/)
* Python helper scripts for plotting:
  [myplot](https://github.com/fluhus/lab-utils/blob/master/myplot.py) and
  [cmbfig](https://github.com/fluhus/lab-utils/blob/master/cmbfig.py)
* Optional: a snapshot of the generated test data it available on
  [Zenodo](https://zenodo.org/records/17904728).

### How to run

1. Set the variables in `run_comparisons.sh`.
2. Add myplot's directory to `PYTHONPATH`.
3. Add MMSeqs and Sourmash to `PATH`.
4. Run the commands in `run_comparisons.sh` from the project's root directory,
   preferably one by one to make sure everything runs as expected.

A temporary directory `testdata/fasta` will be created:
* `vir_*.fa`: Sequences for viral search benchmark.
* `mut_*.fa`: Mutated sequences for viral search benchmark.
* `big.fa`: Sequences for bacterial search benchmark.
* `clust_frag.fa`: Sequences for fragment-clustering benchmark.
* `clust_snps.fa`: Sequences for mutated-clustering benchmark.
