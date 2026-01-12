package sketching

import (
	"fmt"
	"math"

	"github.com/fluhus/gostuff/sets"
	"golang.org/x/exp/maps"
)

// Index allows quick lookups for sketches.
type Index struct {
	idx svmap[uint64, int]
	mx  uint64
}

// NewIndex returns a new index that stores 1/scale of hashes.
func NewIndex(scale uint64) *Index {
	return &Index{
		idx: newSVMap[uint64, int](),
		mx:  math.MaxUint64 / scale,
	}
}

// Add adds the given sketch with the given serial number.
func (idx *Index) Add(s []uint64, i int) {
	for _, x := range s {
		if x > idx.mx {
			break
		}
		idx.idx.put(x, i)
	}
}

// Search returns serial numbers of sketches that share hashes with
// the given sketch.
func (idx *Index) Search(s []uint64) []int {
	set := sets.Set[int]{}
	for _, x := range s {
		if x > idx.mx {
			break
		}
		for i := range idx.idx.get(x) {
			set.Add(i)
		}
	}
	return maps.Keys(set)
}

// Clean removes keys with only one element.
// Use only for clustering.
func (idx *Index) Clean() {
	n1 := len(idx.idx.singles)
	n2 := len(idx.idx.slices)
	idx.idx.clearSingles()
	idx.idx.singles = maps.Clone(idx.idx.singles) // Reduce memory footprint.
	fmt.Printf("Cleaning: %d ==> %d (%.0f%%)\n",
		n1, n2, float64(n2)/float64(n1)*100)
}
