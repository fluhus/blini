package sketching

import (
	"slices"
	"testing"
)

func TestHashIndex(t *testing.T) {
	idxs := []hashIndex[uint, uint]{
		newSVMap[uint, uint](),
		&flatmap[uint, uint]{},
		&flatmap2[uint, uint]{},
		newFlatmap3[uint, uint](),
	}

	for _, idx := range idxs {
		idx.put(5, 31)
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
		}
		for _, test := range tests {
			if got := slices.Collect(idx.get(test.k)); !slices.Equal(got, test.want) {
				t.Errorf("(%T).get(%d)=%d, want %d", idx, test.k, got, test.want)
			}
		}
	}
}
