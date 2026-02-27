package sketching

import (
	"slices"
	"testing"
)

func TestSVMap(t *testing.T) {
	s := newSVMap[int, int]()
	s.put(5, 31)
	s.put(3, 55)
	s.put(5, 12)
	s.put(6, 90)
	s.put(5, 90)

	tests := []struct {
		k         int
		want      []int
		wantClear []int
	}{
		{3, []int{55}, nil},
		{5, []int{31, 12, 90}, []int{31, 12, 90}},
		{6, []int{90}, nil},
		{7, nil, nil},
	}
	for _, test := range tests {
		if got := slices.Collect(s.get(test.k)); !slices.Equal(got, test.want) {
			t.Errorf("get(%d)=%d, want %d", test.k, got, test.want)
		}
	}
	s.clean()
	for _, test := range tests {
		if got := slices.Collect(s.get(test.k)); !slices.Equal(got, test.wantClear) {
			t.Errorf("get(%d)=%d, want %d", test.k, got, test.wantClear)
		}
	}
}
