package sketching

import (
	"slices"
	"testing"
)

func TestSort2(t *testing.T) {
	a := []int{5, 3, 7, 1, 9, 6}
	b := []string{"b", "g", "e", "f", "v", "s"}
	wanta := []int{1, 3, 5, 6, 7, 9}
	wantb := []string{"f", "g", "b", "s", "e", "v"}
	gota := slices.Clone(a)
	gotb := slices.Clone(b)
	sort2(gota, gotb)
	if !slices.Equal(gota, wanta) {
		t.Errorf("sort2(%v)=%v, want %v", a, gota, wanta)
	}
	if !slices.Equal(gotb, wantb) {
		t.Errorf("sort2(%v)=%v, want %v", b, gotb, wantb)
	}
}
