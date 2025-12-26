// Reusable code for the blini algorithm.
package libblini

import (
	"iter"
	"slices"

	"github.com/fluhus/blini/sketching"
)

const (
	kmerLen     = 21
	idxScale    = 4
	indexSuffix = ".blini" // Suffix of pre-sketched files.
)

// Holds sketches, metadata and index.
type Dataset struct {
	sketches [][]uint64       // Frac min hash sketches.
	lens     []int            // Sequence lengths.
	names    []string         // Sequence names.
	scale    uint64           // Kmer selection scale.
	idx      *sketching.Index // Search index.
}

type SearchResult struct {
	I          int     // Reference serial number, 0-based.
	Name       string  // Reference name.
	Similarity float64 // Between 0-1.
}

func (d *Dataset) Search(seq []byte, contn bool) iter.Seq[SearchResult] {
	return d.SearchSketch(
		sketching.Sketch(seq, kmerLen, d.scale), len(seq), contn)
}

func (d *Dataset) SearchSketch(s []uint64, ln int, contn bool) iter.Seq[SearchResult] {
	return func(yield func(SearchResult) bool) {
		for _, i := range d.idx.Search(s) {
			sim := similarity(s, d.sketches[i], ln, d.lens[i], contn)
			sr := SearchResult{
				I:          i,
				Name:       d.names[i],
				Similarity: sim,
			}
			if !yield(sr) {
				return
			}
		}
	}
}

func (d *Dataset) CleanIndex() {
	d.idx.Clean()
}

func (d *Dataset) Name(i int) string {
	return d.names[i]
}

func (d *Dataset) Sketch(i int) []uint64 {
	return slices.Clone(d.sketches[i])
}

func (d *Dataset) SequenceLen(i int) int {
	return d.lens[i]
}

func (d *Dataset) Len() int {
	return len(d.sketches)
}
