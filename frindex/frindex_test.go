package frindex

import (
	"slices"
	"testing"

	"github.com/fluhus/gostuff/minhash"
	"github.com/fluhus/gostuff/snm"
)

func TestIndex_n2(t *testing.T) {
	input := [][]uint64{
		{1, 3, 5, 7, 9},
		{1, 3, 4, 7, 9},
		{3, 5, 7, 9, 11},
		{0, 1, 2, 7, 9},
	}
	mhs := snm.SliceToSlice(input, func(u []uint64) *minhash.MinHash[uint64] {
		mh := minhash.New[uint64](5)
		for _, x := range u {
			mh.Push(x)
		}
		return mh.Frozen()
	})
	want := [][]int{{0, 1, 2}, {0, 1}, {0, 2}, {3}}
	idx := New[uint64](3, 2)
	for _, mh := range mhs {
		idx.Add(mh)
	}
	got := snm.SliceToSlice(mhs, func(mh *minhash.MinHash[uint64]) []int {
		return snm.Sorted(idx.Query(mh))
	})
	if !slices.EqualFunc(want, got, slices.Equal) {
		t.Fatalf("New(%v).Query()=%v, want %v", input, got, want)
	}
}

func TestIndex_n1(t *testing.T) {
	input := [][]uint64{
		{1, 3, 5},
		{2, 4, 6},
		{1, 2, 3},
		{3, 5, 10},
		{11, 12, 13},
	}
	mhs := snm.SliceToSlice(input, func(u []uint64) *minhash.MinHash[uint64] {
		mh := minhash.New[uint64](3)
		for _, x := range u {
			mh.Push(x)
		}
		return mh.Frozen()
	})
	want := [][]int{{0, 2, 3}, {1, 2}, {0, 1, 2, 3}, {0, 2, 3}, {4}}
	idx := New[uint64](3, 1)
	for _, mh := range mhs {
		idx.Add(mh)
	}
	got := snm.SliceToSlice(mhs, func(mh *minhash.MinHash[uint64]) []int {
		return snm.Sorted(idx.Query(mh))
	})
	if !slices.EqualFunc(want, got, slices.Equal) {
		t.Fatalf("New(%v).Query()=%v, want %v", input, got, want)
	}
}
