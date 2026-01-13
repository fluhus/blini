package sketching

import (
	"fmt"
	"slices"

	"github.com/fluhus/gostuff/sets"
	"golang.org/x/exp/constraints"
	"golang.org/x/exp/maps"
)

type idType = uint32

// Index allows quick lookups for sketches.
type Index[T constraints.Unsigned] struct {
	idx   svmap[T, idType]
	scale int
}

// NewIndex returns a new index that stores 1/scale of hashes.
func NewIndex[T constraints.Unsigned](scale int) *Index[T] {
	return &Index[T]{
		idx:   newSVMap[T, idType](),
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
		idx.idx.put(x, idType(i))
	}
}

// Search returns serial numbers of sketches that share hashes with
// the given sketch.
func (idx *Index[T]) Search(s []T) []int {
	set := sets.Set[int]{}
	for _, x := range s {
		for i := range idx.idx.get(x) {
			set.Add(int(i))
		}
	}
	return maps.Keys(set)
}

// Clean removes keys with only one element.
// Use only for clustering.
func (idx *Index[T]) Clean() {
	n1 := len(idx.idx.singles)
	n2 := len(idx.idx.slices)
	idx.idx.clearSingles()
	idx.idx.singles = maps.Clone(idx.idx.singles) // Reduce memory footprint.
	fmt.Printf("Cleaning: %d ==> %d (%.0f%%)\n",
		n1, n2, float64(n2)/float64(n1)*100)
}
