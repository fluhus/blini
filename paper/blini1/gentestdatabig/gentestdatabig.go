// Generates the big test dataset from the bacterial reference.
package main

import (
	"fmt"
	"math/rand/v2"
	"os"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/blini/paper/simul"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/reservoir"
)

const (
	nSequences    = 100000
	nSequencesRes = 10
	seqLen        = 10000

	outDir = "testdata/fasta"
)

func main() {
	os.MkdirAll(outDir, 0o744)

	fmt.Println("Counting sequences")
	refFile := os.Args[1]
	n := 0
	for fa, err := range fasta.File(refFile) {
		if err != nil {
			panic(err)
		}
		if len(fa.Sequence) < seqLen {
			continue
		}
		n++
	}
	fmt.Println(n)

	dist := createSequenceDistribution(n)
	res := reservoir.New[*fasta.Fasta](nSequencesRes)

	fmt.Println("Sampling")
	fout, err := aio.Create(outDir + "/big.fa")
	if err != nil {
		panic(err)
	}
	defer fout.Close()

	i := 0
	for fa, err := range fasta.File(refFile) {
		if err != nil {
			panic(err)
		}
		if len(fa.Sequence) < seqLen {
			continue
		}
		nsubseqs := len(fa.Sequence) - seqLen + 1
		for range dist[i] {
			ii := rand.IntN(nsubseqs)
			fa1 := &fasta.Fasta{
				Name:     fa.Name,
				Sequence: simul.MutSeq(fa.Sequence[ii:ii+seqLen], seqLen/1000),
			}
			fa1.Write(fout)
			res.Add(fa1)
		}
		i++
	}

	for i, fa := range res.Elements {
		txt, _ := fa.MarshalText()
		f := fmt.Sprintf("%s/big_%d.fa", outDir, i+1)
		os.WriteFile(f, txt, 0o644)
	}
}

func createSequenceDistribution(n int) []int {
	result := make([]int, n)
	for range nSequences {
		result[rand.IntN(len(result))]++
	}
	return result
}
