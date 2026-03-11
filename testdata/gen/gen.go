// Generates test data for Blini's CI tests.
package main

import (
	"fmt"
	"math/rand/v2"
	"slices"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/blini/paper/simul"
	"github.com/fluhus/gostuff/aio"
)

const (
	outDir = "testdata"
)

func main() {
	if err := genClustData(); err != nil {
		panic(err)
	}
	if err := genSearchData(); err != nil {
		panic(err)
	}
	if err := genDumpData(); err != nil {
		panic(err)
	}
}

// Generates dummy data for testing the clustering function.
func genClustData() error {
	const (
		nClusters = 3
		nElements = 3
		minLen    = 5000
		maxLen    = 10000
	)

	var fas []*fasta.Fasta
	for i := range nClusters {
		// Create cluster center.
		fa := &fasta.Fasta{
			Name:     fmt.Appendf(nil, "ref%d", i+1),
			Sequence: simul.RandSeq(minLen + rand.IntN(maxLen-minLen+1)),
		}
		fas = append(fas, fa)

		// Create elements.
		for j := range nElements {
			ssLen := len(fa.Sequence)/2 + rand.IntN(len(fa.Sequence)/2) - 1
			ss := &fasta.Fasta{
				Name:     fmt.Appendf(nil, "ref%d.%d", i+1, j+1),
				Sequence: simul.RandSubseq(fa.Sequence, ssLen),
			}
			ss.Sequence = simul.MutSeqPerc(ss.Sequence, 1)
			fas = append(fas, ss)
		}
	}

	shuffle(fas)

	fout, err := aio.Create(outDir + "/clust.fa.zst")
	if err != nil {
		return err
	}
	defer fout.Close()
	for _, fa := range fas {
		if err := fa.Write(fout); err != nil {
			return err
		}
	}
	return nil
}

// Generates dummy data for testing the search function.
func genSearchData() error {
	const (
		nRefs    = 3
		nQueries = 3
		minLen   = 5000
		maxLen   = 10000
	)

	var refs []*fasta.Fasta
	var queries []*fasta.Fasta
	for i := range nRefs {
		// Create reference sequence.
		fa := &fasta.Fasta{
			Name:     fmt.Appendf(nil, "ref%d", i+1),
			Sequence: simul.RandSeq(minLen + rand.IntN(maxLen-minLen+1)),
		}
		refs = append(refs, fa)

		// Create queries.
		for j := range nQueries {
			qLen := len(fa.Sequence)/2 + rand.IntN(len(fa.Sequence)/2) - 1
			q := &fasta.Fasta{
				Name:     fmt.Appendf(nil, "query%d.%d", i+1, j+1),
				Sequence: simul.RandSubseq(fa.Sequence, qLen),
			}
			q.Sequence = simul.MutSeqPerc(q.Sequence, 2)
			queries = append(queries, q)
		}
	}

	// Add some queries that match nothing.
	for i := range nQueries {
		// Half the length of a ref sequence.
		ssLen := (minLen + rand.IntN(maxLen-minLen+1)) / 2
		q := &fasta.Fasta{
			Name:     fmt.Appendf(nil, "query%d.%d", nRefs+1, i+1),
			Sequence: simul.RandSeq(ssLen),
		}
		queries = append(queries, q)
	}

	shuffle(refs)
	shuffle(queries)

	rout, err := aio.Create(outDir + "/refs.fa.zst")
	if err != nil {
		return err
	}
	defer rout.Close()
	for _, fa := range refs {
		if err := fa.Write(rout); err != nil {
			return err
		}
	}

	qout, err := aio.Create(outDir + "/queries.fa.zst")
	if err != nil {
		return err
	}
	defer qout.Close()
	for _, fa := range queries {
		if err := fa.Write(qout); err != nil {
			return err
		}
	}
	return nil
}

// Generates dummy data for testing the distance-dump function.
func genDumpData() error {
	const (
		seqLen  = 10000
		mutPerc = 2
		mutNucs = seqLen * mutPerc / 100
	)
	s := newSNPer(seqLen)
	a := simul.RandSeq(seqLen)
	b := s.mutate(a, mutNucs)
	c := s.mutate(a, mutNucs)
	d := s.mutate(b, mutNucs)

	out, err := aio.Create(outDir + "/dump.fa.zst")
	if err != nil {
		return err
	}
	defer out.Close()
	if err := (&fasta.Fasta{[]byte("a"), a}).Write(out); err != nil {
		return err
	}
	if err := (&fasta.Fasta{[]byte("b"), b}).Write(out); err != nil {
		return err
	}
	if err := (&fasta.Fasta{[]byte("c"), c}).Write(out); err != nil {
		return err
	}
	if err := (&fasta.Fasta{[]byte("d"), d}).Write(out); err != nil {
		return err
	}
	return nil
}

// Shuffles a slice.
func shuffle[T any](a []T) {
	rand.Shuffle(len(a), func(i, j int) {
		a[i], a[j] = a[j], a[i]
	})
}

// Creates mutually exclusive SNPs.
type snper struct {
	p []int
}

func newSNPer(n int) *snper {
	return &snper{rand.Perm(n)}
}

func (s *snper) mutate(b []byte, n int) []byte {
	b = slices.Clone(b)
	for _, i := range s.p[:n] {
		b[i] = simul.MutNuc(b[i])
	}
	s.p = s.p[n:]
	return b
}
