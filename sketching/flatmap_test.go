package sketching

import (
	"slices"
	"testing"
)

func TestFlatMap(t *testing.T) {
	s := &flatmap[int, int]{}
	s.put(5, 31)
	s.put(3, 55)
	s.put(5, 12)
	s.put(6, 90)
	s.put(5, 90)
	s.finalize()

	tests := []struct {
		k    int
		want []int
	}{
		{3, []int{55}},
		{5, []int{31, 12, 90}},
		{6, []int{90}},
		{7, nil},
	}
	for _, test := range tests {
		if got := slices.Collect(s.get(test.k)); !slices.Equal(got, test.want) {
			t.Errorf("get(%d)=%d, want %d", test.k, got, test.want)
		}
	}
}
