# Blini

Lightweight nucleotide sequence clustering and searching.

## Requirements

[FastANI](https://github.com/ParBLiSS/FastANI/releases/)
in the system PATH.
FastANI for Windows is available under Releases.

## Usage

### Clustering

```
blini -i input.fasta -o output.json [optional flags]
```

### Searching

```
blini -i query.fasta -r reference.fasta -o output.json [optional flags]
```

### Output

Blini creates two output files: one by fasta name and one by number,
which is the 0-based index of the sequence in the input file.

### Help

```
blini -h
```

## Limitations

* Blini supports nucleotide sequences only. Amino-acids are not supported.
* Blini runs on a single file with sequences,
  where each sequence is a separate species.
  Support for multiple files and multiple sequences per species
  will be added in the future.
