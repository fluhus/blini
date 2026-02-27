package sketching

import (
	"fmt"
	"iter"
	"maps"
	"os"
)

// A map with slice values, optimized for cases where most
// values are singletons.
type svmap[K comparable, V any] struct {
	singles map[K]V
	slices  map[K][]V
}

// Creates a new empty svmap.
func newSVMap[K comparable, V any]() *svmap[K, V] {
	return &svmap[K, V]{
		singles: map[K]V{},
		slices:  map[K][]V{},
	}
}

// Appends v to the values of k.
func (s *svmap[K, V]) put(k K, v V) {
	if _, ok := s.singles[k]; ok {
		s.slices[k] = append(s.slices[k], v)
	} else {
		s.singles[k] = v
	}
}

// Yields the elements associated with k.
func (s *svmap[K, V]) get(k K) iter.Seq[V] {
	return func(yield func(V) bool) {
		if v, ok := s.singles[k]; ok {
			if !yield(v) {
				return
			}
			for _, v := range s.slices[k] {
				if !yield(v) {
					return
				}
			}
		}
	}
}

// Removes keys with one value.
func (s *svmap[K, V]) clean() {
	n1 := len(s.singles)
	n2 := len(s.slices)
	for k := range s.singles {
		if s.slices[k] == nil {
			delete(s.singles, k)
		}
	}
	s.singles = maps.Clone(s.singles) // Reduce memory footprint.
	fmt.Fprintf(os.Stderr, "Cleaning: %d ==> %d (%.0f%%)\n",
		n1, n2, float64(n2)/float64(n1)*100)
}

// No-op.
func (s *svmap[K, V]) finalize() {}

var _ hashIndex[int, int] = &svmap[int, int]{}
