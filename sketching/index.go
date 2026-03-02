package sketching

import (
	"cmp"
	"fmt"
	"iter"
	"slices"

	"github.com/fluhus/gostuff/sets"
	"golang.org/x/exp/constraints"
	"golang.org/x/exp/maps"
)

type idType = uint32

// Index allows quick lookups for sketches.
type Index[T constraints.Unsigned] struct {
	idx   hashIndex[T, idType]
	scale int
}

// NewIndex returns a new index that stores 1/scale of hashes.
func NewIndex[T constraints.Unsigned](scale int) *Index[T] {
	return &Index[T]{
		idx:   newFlatmap[T, idType](),
		scale: scale,
	}
}

// Add adds the given sketch with the given serial number.
func (idx *Index[T]) Add(s []T, i int) {
	if !slices.IsSorted(s) {
		panic(fmt.Sprintf("unsorted sketch: %v", s))
	}
	upto := (len(s) + idx.scale - 1) / idx.scale
	for _, x := range s[:upto] {
		// BUG(fluhus): Consider making put receive a slice.
		idx.idx.put(x, idType(i))
	}
}

// Search returns serial numbers of sketches that share hashes with
// the given sketch.
func (idx *Index[T]) Search(s []T) []int {
	if minIndexHitRatio > 1 {
		cnt := map[idType]int{}
		for _, x := range s {
			for i := range idx.idx.get(x) {
				cnt[i]++
			}
		}
		var result []int
		r := minIndexHitRatio * idx.scale
		min := (len(s) + r - 1) / r
		for k, v := range cnt {
			if v >= min {
				result = append(result, int(k))
			}
		}
		return result
	}
	set := sets.Set[int]{}
	for _, x := range s {
		for i := range idx.idx.get(x) {
			set.Add(int(i))
		}
	}
	return maps.Keys(set)
}

// Finalize runs final processing after adding data and before using the index.
func (idx *Index[T]) Finalize() {
	idx.idx.finalize()
}

// A common interface for the index data structures.
// Used for testing different indexes.
type hashIndex[K cmp.Ordered, V any] interface {
	put(K, V)          // Adds a hash to an ID.
	get(K) iter.Seq[V] // Returns an iterator of IDs for a hash.
	finalize()         // Finalize index construction before use.
}
