// Reusable code for the blini algorithm.
package libblini

import (
	"iter"

	"github.com/fluhus/blini/sketching"
	"golang.org/x/exp/constraints"
)

const (
	kmerLen     = 21       // K-mer length.
	idxScale    = 4        // What ratio of the hashes is used for the index.
	indexSuffix = ".blini" // Suffix of pre-sketched files.
)

// Dataset holds sketches, metadata and index.
type Dataset[T constraints.Unsigned] struct {
	sketches [][]T               // Frac min hash sketches.
	lens     []int               // Sequence lengths.
	names    []string            // Sequence names.
	scale    uint64              // Kmer selection scale.
	idx      *sketching.Index[T] // Search index.
}

// SearchResult is a single match between a query and a reference sequence.
type SearchResult struct {
	I          int     // Reference serial number, 0-based.
	Name       string  // Reference name.
	Similarity float64 // Between 0-1.
}

// Search sketches seq and and looks it up in the dataset.
func (d *Dataset[T]) Search(seq []byte, contn bool) iter.Seq[SearchResult] {
	s := &Sketch[T]{
		h:     sketching.Sketch[T](seq, kmerLen, d.scale),
		ln:    len(seq),
		scale: d.scale,
	}
	return d.SearchSketch(s, contn)
}

// SearchSketch looks up a sketch in the dataset.
func (d *Dataset[T]) SearchSketch(s *Sketch[T], contn bool) iter.Seq[SearchResult] {
	if d.idx == nil {
		panic("this dataset was created with no search index")
	}
	return func(yield func(SearchResult) bool) {
		for _, i := range d.idx.Search(s.h) {
			sim := similarity(s.h, d.sketches[i], s.ln, d.lens[i], contn)
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
func (d *Dataset[T]) CleanIndex() {
	d.idx.Clean()
}

// Sketch returns a clone of the i'th sketch.
func (d *Dataset[T]) Sketch(i int) *Sketch[T] {
	return &Sketch[T]{
		h:     d.sketches[i],
		ln:    d.lens[i],
		name:  d.names[i],
		scale: d.scale,
	}
}

// Len returns the number of sketches in the dataset.
func (d *Dataset[T]) Len() int {
	return len(d.sketches)
}
