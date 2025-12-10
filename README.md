# Blini

Lightweight nucleotide sequence searching and dereplication.

## Requirements

None.
[Download](https://github.com/fluhus/blini/releases)
and get started!

## Usage (basic)

### Searching

With both `-q` and `-r` set, Blini looks up the query entries in the
reference entries.
The reference may either be a fasta or a pre-sketched index.

```sh
blini -q query.fasta -r reference.fasta -o output.csv
# Or
blini -q query.fasta -r reference.blini -o output.csv
```

### Sketching

With only `-r` set, Blini pre-sketches the given reference for use
in search operations.
This makes lookup operations quicker.

```sh
blini -r reference.fasta -o reference.blini
```

### Clustering

With only `-q` set, Blini dereplicates (clusters) the query set.

```sh
blini -q input.fasta -o output_prefix
```

The outputs are a fasta file with the representatives,
and a JSON file with the cluster assignments.

### Other options

* `-h` display help on the available flags.
* `-c` for searching, calculate containment of query in the reference
  rather than full match.
* `-m` for searching and clustering,
  minimal similarity for a match.
* `-s` scale; use 1/s of kmers for similarity.
* `-u` for search, include unmatched queries in the output.

## Usage (advanced)

### Choosing the scale value (`-s`)

**Scale should be at most 1/25 the length of the sequnces analyzed.**

*Scale* is the k-mer subsampling ratio.
A scale of 100 means that 1/100 of the k-mers are used for distance calculations.
Doubling the scale halves RAM and CPU usage, but also loses some accuracy.
For accurate results, sketches of size 25 and above are needed.
This means that the scale needs to be up to `sequence length / 25`.

The default scale of 100 is effective for sequences of length 2500 and above.
For sequences of length 1000, for example, the scale needs to be at most 40.

### Parallelizing reference sketching

Sketch files (`.blini`) can be concatenated
if they were created using the same scale.
This is equivalent to having the different original datasets sketched
together in one run.
Therefore, big reference datasets can be broken down and sketched in parallel.

## Limitations

* Blini supports nucleotide sequences only.
  Amino-acids are currently not supported.
* Blini runs on a single file with sequences,
  where each sequence is a separate species.
  Support for multiple files and multiple sequences per species
  will be added in the future.
* No multi-threading at the moment.
  Still fast, innit?

## Testing

The
[testdata](https://github.com/fluhus/blini/tree/main/testdata)
directory contains mock data for testing Blini.
The tests run automatically with each push to this repository,
but you can also run them locally.
See the directory's README for further instructions.

## Join the conversation

Got a question? Feedback? Found a bug?

Feel free to
[open an issue](https://github.com/fluhus/blini/issues/new/choose)
or
[start a discussion](https://github.com/fluhus/blini/discussions/new/choose).
