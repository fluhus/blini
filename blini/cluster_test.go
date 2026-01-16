package main

import (
	"cmp"
	"iter"
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
	input := [][]hashType{
		{1, 4, 2},
		{2},
		{1, 3, 4},
		{5, 1, 3, 2},
	}
	want := []int{5, 2, 4, 5}
	got := babiScores(func() iter.Seq[iter.Seq[hashType]] {
		return func(yield func(iter.Seq[hashType]) bool) {
			for _, x := range input {
				if !yield(slices.Values(x)) {
					break
				}
			}
		}
	})
	if !slices.Equal(got, want) {
		t.Fatalf("babiScores(%v)=%v, want %v", input, got, want)
	}
}
