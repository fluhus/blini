// Package frindex implements a search index on MinHash sketches.
package frindex

import (
	"github.com/fluhus/gostuff/minhash"
	"github.com/fluhus/gostuff/sets"
	"golang.org/x/exp/constraints"
	"golang.org/x/exp/maps"
)

// Index is a search index on MinHash sketches.
type Index[H constraints.Integer] struct {
	k int         // Number of values to look at
	n int         // Min common elements
	i int         // Serial of the next addition
	m map[H][]int // Hashes map
}

// New returns a new index that looks for k matches out of the
// n minimal elements.
func New[H constraints.Integer](k, n int) *Index[H] {
	return &Index[H]{k, n, 0, map[H][]int{}}
}

// Add adds the given sketches to the index.
func (idx *Index[H]) Add(mhs ...*minhash.MinHash[H]) {
	for _, mh := range mhs {
		v := mh.View()
		v = v[max(len(v)-idx.k, 0):]
		for _, x := range v {
			idx.m[x] = append(idx.m[x], idx.i)
		}
		idx.i++
	}
}

// Query returns the serial numbers of sketches that share
// at least k elements with the given one.
func (idx *Index[H]) Query(mh *minhash.MinHash[H]) []int {
	if idx.n == 1 { // Optimization for n=1.
		found := make(sets.Set[int], idx.i*idx.k/len(idx.m)*2)
		v := mh.View()
		v = v[max(len(v)-idx.k, 0):]
		for _, x := range v {
			found.Add(idx.m[x]...)
		}
		return maps.Keys(found)
	}
	cnt := make(map[int]int, idx.i*idx.k/len(idx.m)*2)
	v := mh.View()
	v = v[max(len(v)-idx.k, 0):]
	for _, x := range v {
		for _, f := range idx.m[x] {
			cnt[f]++
		}
	}
	var result []int
	for k, v := range cnt {
		if v >= idx.n {
			result = append(result, k)
		}
	}
	return result
}
