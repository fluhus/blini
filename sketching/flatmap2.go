package sketching

import (
	"cmp"
	"iter"
	"slices"
)

// A multi-map backed by sorted slices of of keys and values.
type flatmap2[K cmp.Ordered, V any] struct {
	k []K
	v []V
}

// Inserts v to the values of k.
func (f *flatmap2[K, V]) put(k K, v V) {
	f.k = append(f.k, k)
	f.v = append(f.v, v)
}

// Returns an iterator over the values of k.
func (f *flatmap2[K, V]) get(k K) iter.Seq[V] {
	return func(yield func(V) bool) {
		i, _ := slices.BinarySearch(f.k, k)
		for j := i; j < len(f.k); j++ {
			if f.k[j] != k {
				break
			}
			if !yield(f.v[j]) {
				break
			}
		}
	}
}

// Sorts the slice for binary search.
func (f *flatmap2[K, V]) finalize() {
	sort2(f.k, f.v)
}

// No-op.
func (f *flatmap2[K, V]) clean() {}

var _ hashIndex[int, int] = &flatmap2[int, int]{}
