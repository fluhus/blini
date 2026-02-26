package libblini

import (
	"cmp"
	"slices"
	"testing"
)

func TestSortedPerm(t *testing.T) {
	input := []string{"a", "g", "d", "b"}
	want := []int{0, 3, 2, 1}
	got := sortedPerm(len(input), func(i1, i2 int) int {
		return cmp.Compare(input[i1], input[i2])
	})
	if !slices.Equal(got, want) {
		t.Fatalf("sortedPerm(%q)=%v, want %v", input, got, want)
	}
}

func TestBabiScores(t *testing.T) {
	db := &Dataset[uint64]{}
	db.sketches = [][]uint64{
		{1, 4, 2},
		{2},
		{1, 3, 4},
		{5, 1, 3, 2},
	}
	want := []int{5, 2, 4, 5}
	got := db.babiScores()
	if !slices.Equal(got, want) {
		t.Fatalf("babiScores(%v)=%v, want %v", db.sketches, got, want)
	}
}
