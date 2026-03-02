package sketching

import (
	"cmp"
	"sort"
)

// Sorts a and permutes b in the same way.
func sort2[T cmp.Ordered, S any](a []T, b []S) {
	sort.Sort(&sort2t[T, S]{a, b})
}

type sort2t[T cmp.Ordered, S any] struct {
	a []T
	b []S
}

func (s *sort2t[T, S]) Len() int {
	return len(s.a)
}

func (s *sort2t[T, S]) Less(i, j int) bool {
	return s.a[i] < s.a[j]
}

func (s *sort2t[T, S]) Swap(i, j int) {
	s.a[i], s.a[j] = s.a[j], s.a[i]
	s.b[i], s.b[j] = s.b[j], s.b[i]
}
