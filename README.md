# Blini

Lightweight nucleotide sequence clustering and searching.

## Requirements

[fastANI](https://github.com/ParBLiSS/FastANI/releases/)
in the system PATH.

## Usage

### Clustering

```
blini -i input.fasta -o output.json [additional flags]
```

### Searching

```
blini -i query.fasta -r reference.fasta -o output.json [additional flags]
```

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
