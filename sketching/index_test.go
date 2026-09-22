package sketching

import (
	"maps"
	"slices"
	"testing"
)

func TestHashIndex(t *testing.T) {
	// New index structures can be added here.
	idxs := []hashIndex[uint, uint]{
		newFlatmap[uint, uint](),
	}

	for _, idx := range idxs {
		idx.put(5, 31)
		idx.put(5+4096, 111)
		idx.put(3, 55)
		idx.put(5, 12)
		idx.put(6, 90)
		idx.put(5, 90)
		idx.finalize()

		tests := []struct {
			k    uint
			want []uint
		}{
			{3, []uint{55}},
			{5, []uint{31, 12, 90}},
			{6, []uint{90}},
			{7, nil},
			{5 + 4096, []uint{111}},
		}
		for _, test := range tests {
			if got := slices.Collect(idx.get(test.k)); !slices.Equal(got, test.want) {
				t.Errorf("(%T).get(%d)=%d, want %d", idx, test.k, got, test.want)
			}
		}

		wantNVals := map[uint]int{3: 1, 5: 3, 6: 1, 5 + 4096: 1}
		gotNVals := map[uint]int{}
		for k, n := range idx.nvals() {
			if _, ok := gotNVals[k]; ok {
				t.Errorf("(%T).nvals(): duplicate key: %v", idx, k)
			}
			gotNVals[k] = n
		}
		if !maps.Equal(gotNVals, wantNVals) {
			t.Errorf("(%T).nvals()=%v, want %v", idx, gotNVals, wantNVals)
		}
	}
}
