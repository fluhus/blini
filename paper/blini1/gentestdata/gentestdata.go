// Picks out random genomes from an input reference.
package main

import (
	"fmt"
	"math/rand/v2"
	"os"
	"regexp"
	"slices"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/blini/paper/simul"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/gnum"
	"github.com/fluhus/gostuff/reservoir"
)

const (
	nSampled   = 100
	nFragments = 300
	nMutants   = 100
	minFragLen = 1000
	minSeqLen  = 10000

	outDir = "testdata/fasta"
)

func main() {
	fmt.Println("Sampling sequences")
	seqs, err := sampleSequences(os.Args[1])
	if err != nil {
		panic(err)
	}

	fmt.Println("Generating test dataset")
	mustMkdir(outDir)
	f, err := aio.Create(outDir + "/vir_all.fa")
	if err != nil {
		panic(err)
	}
	defer f.Close()
	fm, err := aio.Create(outDir + "/mut_all.fa")
	if err != nil {
		panic(err)
	}
	defer fm.Close()
	fcf, err := aio.Create(outDir + "/clust_frag.fa")
	if err != nil {
		panic(err)
	}
	defer fcf.Close()
	fcs, err := aio.Create(outDir + "/clust_snps.fa")
	if err != nil {
		panic(err)
	}
	defer fcs.Close()

	nameRE := regexp.MustCompile(`^\S+`)
	var lens []int
	for i, fa := range seqs {
		// Output for searching.
		fa.Name = nameRE.Find(fa.Name)
		txt, _ := fa.MarshalText()
		file := fmt.Sprintf("%s/vir_%d.fa", outDir, i+1)
		os.WriteFile(file, txt, 0o600)
		f.Write(txt)

		fa2 := &fasta.Fasta{Name: fa.Name, Sequence: simul.MutSeqPerc(fa.Sequence, 1)}
		txt, _ = fa2.MarshalText()
		file = fmt.Sprintf("%s/mut_%d.fa", outDir, i+1)
		os.WriteFile(file, txt, 0o600)
		fm.Write(txt)

		// Output for clustering.
		fa3 := &fasta.Fasta{Name: append(fa.Name, []byte("..0")...), Sequence: fa.Sequence}
		fa3.Write(fcf)
		fa3.Write(fcs)
		lens = append(lens, len(fa3.Sequence))
		for i := range nFragments {
			if len(fa.Sequence) < minFragLen {
				panic(fmt.Sprintf("seq too short: %d", len(fa.Sequence)))
			}
			ln := minFragLen + rand.IntN(len(fa.Sequence)+1-minFragLen)
			fa4 := &fasta.Fasta{
				Name:     fmt.Append(fa.Name, "..", i+1),
				Sequence: simul.RandSubseq(fa.Sequence, ln),
			}
			lens = append(lens, len(fa4.Sequence))
			fa4.Write(fcf)
		}
		for i := range nMutants {
			fa4 := &fasta.Fasta{
				Name:     fmt.Append(fa.Name, "..", i+1),
				Sequence: simul.MutSeqPerc(fa.Sequence, 1),
			}
			fa4.Write(fcs)
		}
	}
	slices.Sort(lens)
	fmt.Println("Fragment lengths:", gnum.NQuantiles(lens, 10))
}

// Reads a subsample of the sequences.
func sampleSequences(inFile string) ([]*fasta.Fasta, error) {
	r := reservoir.New[*fasta.Fasta](nSampled)
	for fa, err := range fasta.File(inFile) {
		if err != nil {
			return nil, err
		}
		if len(fa.Sequence) < minSeqLen {
			continue
		}
		r.Add(fa)
	}
	if len(r.Elements) < nSampled {
		return nil, fmt.Errorf("reservoir too small: %v", len(r.Elements))
	}
	return r.Elements, nil
}

// Calls mkdir and panics if it fails.
func mustMkdir(path string) {
	err := os.MkdirAll(path, 0o744)
	if err != nil {
		panic(err)
	}
}
