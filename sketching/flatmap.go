package sketching

import (
	"iter"
	"slices"

	"github.com/fluhus/gostuff/snm"
	"golang.org/x/exp/constraints"
)

const (
	flatmapSize = 4096
)

// A multi-map backed by sorted slices of of keys and values.
type flatmap[K constraints.Unsigned, V any] struct {
	k [][]K
	v [][]V
}

// Returns a newly allocated flatmap.
func newFlatmap[K constraints.Unsigned, V any]() *flatmap[K, V] {
	return &flatmap[K, V]{
		k: make([][]K, flatmapSize),
		v: make([][]V, flatmapSize),
	}
}

// Inserts v to the values of k.
func (f *flatmap[K, V]) put(k K, v V) {
	a := k % K(len(f.k))
	f.k[a] = append(f.k[a], k)
	f.v[a] = append(f.v[a], v)
}

// Returns an iterator over the values of k.
func (f *flatmap[K, V]) get(k K) iter.Seq[V] {
	return func(yield func(V) bool) {
		a := k % K(len(f.k))
		kk := f.k[a]
		vv := f.v[a]
		i, _ := slices.BinarySearch(kk, k)
		for j := i; j < len(kk); j++ {
			if kk[j] != k {
				break
			}
			if !yield(vv[j]) {
				break
			}
		}
	}
}

// Sorts the slice for binary search.
func (f *flatmap[K, V]) finalize() {
	for i := range f.k {
		sort2(f.k[i], f.v[i])

		// Reallocate because dynamically grown slices can theoretically
		// take up to twice their length.
		f.k[i] = snm.TightClone(f.k[i])
		f.v[i] = snm.TightClone(f.v[i])
	}
}

var _ hashIndex[uint, int] = (*flatmap[uint, int])(nil)
