package sketching

import (
	"slices"
	"testing"
)

func TestSVMapClean(t *testing.T) {
	s := newSVMap[int, int]()
	s.put(5, 31)
	s.put(3, 55)
	s.put(5, 12)
	s.put(6, 90)
	s.put(5, 90)

	tests := []struct {
		k    int
		want []int
	}{
		{3, nil},
		{5, []int{31, 12, 90}},
		{6, nil},
		{7, nil},
	}
	s.clean()
	for _, test := range tests {
		if got := slices.Collect(s.get(test.k)); !slices.Equal(got, test.want) {
			t.Errorf("get(%d)=%d, want %d", test.k, got, test.want)
		}
	}
}
