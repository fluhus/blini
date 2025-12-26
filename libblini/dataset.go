// Reusable code for the blini algorithm.
package libblini

import (
	"iter"
	"slices"

	"github.com/fluhus/blini/sketching"
)

const (
	kmerLen     = 21       // K-mer length.
	idxScale    = 4        // What ratio of the hashes is used for the index.
	indexSuffix = ".blini" // Suffix of pre-sketched files.
)

// Dataset holds sketches, metadata and index.
type Dataset struct {
	sketches [][]uint64       // Frac min hash sketches.
	lens     []int            // Sequence lengths.
	names    []string         // Sequence names.
	scale    uint64           // Kmer selection scale.
	idx      *sketching.Index // Search index.
}

// SearchResult is a single match between a query and a reference sequence.
type SearchResult struct {
	I          int     // Reference serial number, 0-based.
	Name       string  // Reference name.
	Similarity float64 // Between 0-1.
}

// Search sketches seq and and looks it up in the dataset.
func (d *Dataset) Search(seq []byte, contn bool) iter.Seq[SearchResult] {
	return d.SearchSketch(
		sketching.Sketch(seq, kmerLen, d.scale), len(seq), contn)
}

// SearchSketch looks up a sketch in the dataset.
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

// CleanIndex removes singletons from the index.
// Used in clustering operations to reduce memory consumption.
func (d *Dataset) CleanIndex() {
	d.idx.Clean()
}

// Name returns the name of the i'th sequence.
func (d *Dataset) Name(i int) string {
	return d.names[i]
}

// Sketch returns a clone of the i'th sketch.
func (d *Dataset) Sketch(i int) []uint64 {
	return slices.Clone(d.sketches[i])
}

// SequenceLen returns the length of the i'th sequence.
func (d *Dataset) SequenceLen(i int) int {
	return d.lens[i]
}

// Len returns the number of sketches in the dataset.
func (d *Dataset) Len() int {
	return len(d.sketches)
}
