package libblini

import (
	"iter"
	"slices"

	"golang.org/x/exp/constraints"
)

// Sketch represents a single sequence that has been sketched.
type Sketch[T constraints.Unsigned] struct {
	h     []T    // Sketch hashes.
	ln    int    // Sequence length.
	name  string // Sequence name.
	scale uint64 // Kmer selection scale.
}

// Hashes returns an iterator of this sketch's elements.
func (s *Sketch[T]) Hashes() iter.Seq[T] {
	return slices.Values(s.h)
}

// Len returns the length of the original sequence.
func (s *Sketch[T]) Len() int {
	return s.ln
}

// Name returns the name of the original sequence.
func (s *Sketch[T]) Name() string {
	return s.name
}

// Scale returns the sketch's scale
// (one out of how many hashes are used).
func (s *Sketch[T]) Scale() uint64 {
	return s.scale
}

// Similarity returns an estimation of average nucleotide identity (ANI)
// between the sequence represented by this sketch and the other.
func (s *Sketch[T]) Similarity(other *Sketch[T], contn bool) float64 {
	return similarity(s.h, other.h, s.ln, other.ln, contn)
}
