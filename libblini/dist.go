package libblini

import (
	"github.com/fluhus/biostuff/mash"
	"github.com/fluhus/blini/sketching"
)

const (
	useMyDist = true // Use a new experiemental distance func.
)

func similarity(s1, s2 []uint64, l1, l2 int, contn bool) float64 {
	if useMyDist {
		return 1 - myDist(s1, s2, l1, l2, contn)
	} else {
		return 1 - mash.FromJaccard(jaccard(s1, s2, contn), kmerLen)
	}
}

// Returns the jaccard/containment jaccard.
func jaccard(a, b []uint64, contn bool) float64 {
	if contn {
		return sketching.Containment(a, b)
	} else {
		return sketching.Jaccard(a, b)
	}
}

// Returns a specialized distance.
func myDist(s1, s2 []uint64, l1, l2 int, contn bool) float64 {
	if contn {
		return mash.FromJaccard(sketching.Containment(s1, s2), kmerLen)
	} else {
		return sketching.MyDist(s1, s2, l1, l2, kmerLen)
	}
}
