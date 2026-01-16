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

<details>
<summary>Output file format</summary>

```
similarity,query,reference
99%,Query sequence 1,Influenza virus
97%,Query sequence 2,Pepper mottle virus
...
```

The first row is always `similarity,query,reference`
and then a row for each match between a query sequence and a reference
sequence that passed the similarity threshold.
There can be several matches per query, and several matches per reference.

</details>

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

<details>
<summary>JSON file format</summary>

```json
{
  "byName": [
    ["Coronavirus 1", "Coronavirus 2"],
    ["Influenza virus 1", "Influenza virus 2"]
  ],
  "byNumber": [
    [1, 4],
    [3, 2]
  ]
}
```

The `byName` value holds the names of sequences in each cluster,
with each cluster's representative first.
The `byNumber` value is the same, with the index of each sequence
in the input.

</details>

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

## Examples

<details>
<summary>Clustering coronavirus species</summary>

Let's download some coronavirus genomes from
[NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=694002&assembly_level=3:3).
We concatenate the sequences into one file, `cov.fa`.

The file should look something like this:

```
>OZ067591.1 Severe acute respiratory syndrome coronavirus 2
CTTTCGATCTCTTGTAGATCTGTTC...
>MW719567.1 Sarbecovirus RhGB01, complete genome
GGAGGATATCACCTGCGGATAAAAG...
>EF424622.1 Giraffe coronavirus US/OH3-TC/2006, complete genome
CGTGCGTGCATCCCGCTTCACTGAT...
```

Run clustering with 99% minimal similarity:

```
blini -q cov.fa -m 0.99 -o cov_clust
```

The output files are `cov_clust.json` and `cov_clust.fasta`.

Inside `cov_clust.json` we find (truncated for clarity):

```json
{
  "byName": [
    [
      "KC164505.2 Betacoronavirus England 1, complete genome",
      "NC_019843.3 Middle East respiratory syndrome-related coronavirus...",
      "NC_038294.1 Betacoronavirus England 1 isolate H123990006, complete genome",
      "KC667074.1 Human betacoronavirus 2c England-Qatar/2012, complete genome"
    ],
    [
      "DQ084199.1 bat SARS coronavirus HKU3-2, complete genome",
      "GQ153540.1 Bat SARS coronavirus HKU3-5, complete genome",
      "GQ153545.1 Bat SARS coronavirus HKU3-10, complete genome",
      "GQ153548.1 Bat SARS coronavirus HKU3-13, complete genome"
    ],
    [
      "NC_017083.1 Rabbit coronavirus HKU14, complete genome",
      "JN874559.1 Rabbit coronavirus HKU14 strain HKU14-1, complete genome"
    ],
    [
      "KF636752.1 Bat Hp-betacoronavirus/Zhejiang2013, complete genome",
      "NC_025217.1 Bat Hp-betacoronavirus/Zhejiang2013, complete genome"
    ]
  ]
}
```

Inside `cov_clust.fasta` we find the first out of each cluster in the JSON file:

```
>KC164505.2 Betacoronavirus England 1, complete genome
...
>DQ084199.1 bat SARS coronavirus HKU3-2, complete genome
...
>NC_017083.1 Rabbit coronavirus HKU14, complete genome
...
>KF636752.1 Bat Hp-betacoronavirus/Zhejiang2013, complete genome
...
```

</details>

<details>
<summary>Identifying mysterious bacteria</summary>

Let's download some bacterial contigs from the
[Segata Lab](http://segatalab.cibio.unitn.it/data/Pasolli_et_al.html).
We'll use
[proGenomes v4](https://progenomes.embl.de/download.cgi)
as the reference.

Pre-sketch the reference - not mandatory but recommended for big datasets.
We will use a scale of 200 rather than the default 100
so the big reference can fit within a personal computer's RAM.
This one-time preprocessing might take about an hour.

```
blini -s 200 -r pg4_genomes_representatives.fna.gz -o pg4.blini
```

Run search with default minimal similarity (90%).
We use `-c` because we expect the query contigs to be subsequences of the
reference whole genomes.

```
blini -q segata_contigs.fa -r pg4.blini -o out.csv -c
```

Voilà! Inside `out.csv` we find (truncated for clarity):

```
similarity,query,reference
99%,gnl|X|NBFFLGAN_1,"Lactobacillus salivarius UCC118, complete genome"
94%,gnl|X|BPBGINNO_1,"Granulicatella elegans ATCC 700633 genomic scaffold supercont2.1, whole genome shotgun sequence"
90%,gnl|X|LPGKDJLJ_1,"Enterobacter bugandensis strain UENF-21GII scaf_2_1030401, whole genome shotgun sequence"
90%,gnl|X|LPGKDJLJ_1,"Lelliottia nimipressuralis strain 51 GCID_CRU_0002_NODE_1.ctg_1, whole genome shotgun sequence"
91%,gnl|X|LPGKDJLJ_1,"Enterobacter bugandensis strain EBG2 NODE_1_length_2688700_cov_52.224477, whole genome shotgun sequence"
98%,gnl|X|LPGKDJLJ_1,"Enterobacter sp. M4-VN DNA, sequence04, whole genome shotgun sequence"
```

</details>

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
