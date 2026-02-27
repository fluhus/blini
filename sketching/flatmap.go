package sketching

import (
	"cmp"
	"iter"
	"slices"
)

// A single key-value pair in the flat map.
type flatmapEntry[K cmp.Ordered, V any] struct {
	k K
	v V
}

// Compares two entries by key.
func cmpEntries[K cmp.Ordered, V any](a, b flatmapEntry[K, V]) int {
	return cmp.Compare(a.k, b.k)
}

// A multi-map backed by a sorted slice of key-value tuples.
type flatmap[K cmp.Ordered, V any] struct {
	s []flatmapEntry[K, V]
}

// Inserts v to the values of k.
func (f *flatmap[K, V]) put(k K, v V) {
	f.s = append(f.s, flatmapEntry[K, V]{k, v})
}

// Returns an iterator over the values of k.
func (f flatmap[K, V]) get(k K) iter.Seq[V] {
	return func(yield func(V) bool) {
		i, _ := slices.BinarySearchFunc(f.s, flatmapEntry[K, V]{k: k}, cmpEntries)
		for _, e := range f.s[i:] {
			if e.k != k {
				break
			}
			if !yield(e.v) {
				break
			}
		}
	}
}

// Sorts the slice for binary search.
func (f *flatmap[K, V]) finalize() {
	slices.SortFunc(f.s, cmpEntries)
}

// No-op.
func (f *flatmap[K, V]) clean() {}

var _ hashIndex[int, int] = &flatmap[int, int]{}
