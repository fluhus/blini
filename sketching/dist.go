// Package sketching provides sequence sketching and
// distance calculation functionality.
package sketching

import (
	"cmp"

	"github.com/fluhus/biostuff/mash"
	"github.com/fluhus/gostuff/sets"
)

const (
	// Punish containment for possible mismatches.
	compensatingCont = true

	// Add a ghost count to the union part of Jaccard calculations.
	// This "punishes" low kmer counts and makes them look less similar.
	ghostUnion = 1
)

// Jaccard returns the Jaccard similarity between a and b.
func Jaccard[T cmp.Ordered](a, b []T) float64 {
	i := sets.SortedIntersectionLen(a, b)
	u := len(a) + len(b) - i
	u += ghostUnion
	return float64(i) / float64(u)
}

// Containment returns a Jaccard-like similarity for the containment
// of a in b.
func Containment[T cmp.Ordered](a, b []T) float64 {
	i := sets.SortedIntersectionLen(a, b)
	u := len(a)
	if compensatingCont {
		u += len(a) - i
	}
	u += ghostUnion
	return float64(i) / float64(u)
}

// MyDist returns a Mash distance with compensation for length
// difference.
func MyDist[T cmp.Ordered](a, b []T, alen, blen int, k int) float64 {
	if alen > blen { // a should be the smaller.
		a, b = b, a
		alen, blen = blen, alen
	}
	r := float64(alen) / float64(blen)
	d := mash.FromJaccard(Containment(a, b), k)
	return d*r + (1 - r)
}
