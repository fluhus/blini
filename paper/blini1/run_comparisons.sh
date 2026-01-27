# Runs benchmarks on blini, sourmash and mmseqs.

datadir=
cmbfigdir=

outdir=testdata/output
fastadir=testdata/fasta

## DATA GENERATION

# Generate simulated dataset.
go run ./paper/gentestdata $datadir/viral.1.1.genomic.fna.gz
go run ./paper/gentestdatabig $datadir/representatives.fa


# PRE-SKETCH REFERENCE DATASET

mkdir -p $outdir
blini -r $datadir/viral.1.1.genomic.fna.gz -o $outdir/viral.blini
sourmash sketch dna --singleton $datadir/viral.1.1.genomic.fna.gz -o $outdir/viral.sm
mmseqs createdb $datadir/viral.1.1.genomic.fna.gz $outdir/viral.mm


## SEARCH TASKS

# Search genomes with blini.
time blini -q $fastadir/vir_all.fa -r $outdir/viral.blini -o $outdir/blini_vir.csv
time blini -q $fastadir/mut_all.fa -r $outdir/viral.blini -o $outdir/blini_mut.csv

# Searches a single file with sourmash.
function smsearch {
  sourmash sketch dna --singleton "$1" -o tmp.sm &> /dev/null &&
  sourmash search -q --ignore-abundance -o "$2" tmp.sm $outdir/viral.sm &&
  rm tmp.sm
}

# Search genomes with sourmash.
time find $fastadir/ -name 'vir_[0-9]*' |
while read f; do
  smsearch $f $outdir/sm_$(basename $f .fa).csv
done

time find $fastadir/ -name 'mut_[0-9]*' |
while read f; do
  smsearch $f $outdir/sm_$(basename $f .fa).csv
done

# Search genomes with mmseqs.
time mmseqs easy-search -v 1 --search-type 3 --threads 1 --min-seq-id 0.9 $fastadir/vir_all.fa $outdir/viral.mm $outdir/mms tmp
time mmseqs easy-search -v 1 --search-type 3 --threads 1 --min-seq-id 0.9 $fastadir/mut_all.fa $outdir/viral.mm $outdir/mmsm tmp

# Evaluate searches.
go run ./paper/testsearch


## CLUSTERING

# Blini cluster.
for s in 25 50 100 200; do
  reset
  for i in 1 2 3 4 5; do
    time blini -q $fastadir/clust_snps.fa -o $outdir/blini_$s -s $s -c -m 0.97
  done
  read -p "Done $s"
done

# MMseqs cluster.
for t in 1 4; do
  reset
  for i in 1 2 3 4 5; do
    rm -fr tmp &&
    time \
      mmseqs easy-linclust -v 1 --threads $t \
      $fastadir/clust_snps.fa $outdir/mm tmp \
      --min-seq-id 0.97 --seq-id-mode 1 --cov-mode 1 &&
    rm -fr tmp
  done
  read -p "Done $t"
done

# Evaluate clustering.
go run ./paper/testclust


## PLOT COMBINING

# Search
python $cmbfigdir/cmbfig.py \
  -c 3 -wt 0.66 -o paper/results/search.png -i \
  $outdir/search_found.png \
  $outdir/search_others.png \
  $outdir/search_time.png

# Cluster frag
python $cmbfigdir/cmbfig.py \
  -c 2 -wt 0.66 -o paper/results/clust_frag.png -i \
  $outdir/clust_frag_nclust.png \
  $outdir/clust_frag_ari.png \
  $outdir/clust_frag_time.png \
  $outdir/clust_frag_mem.png

# Cluster snps
python $cmbfigdir/cmbfig.py \
  -c 2 -wt 0.66 -o paper/results/clust_snps.png -i \
  $outdir/clust_snps_nclust.png \
  $outdir/clust_snps_ari.png \
  $outdir/clust_snps_time.png \
  $outdir/clust_snps_mem.png
