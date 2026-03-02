package libblini

import (
	"cmp"
	"fmt"
	"math"
	"os"
	"slices"
	"time"

	"github.com/fluhus/gostuff/ptimer"
	"github.com/fluhus/gostuff/snm"
)

// Cluster clusters this dataset's sketches,
// with minSim as the minimal similarity between each element and the cluster's
// representative.
// contn is whether containment should be checked rather than full match.
func (d *Dataset[T]) Cluster(minSim float64, contn bool) [][]int {
	perm := sortedPerm(d.Len(), func(a, b int) int {
		return cmp.Compare(d.Sketch(b).Len(), d.Sketch(a).Len())
	})
	if babiClustering {
		perm = d.babiSort()
	}

	friends := 0
	var clusters [][]int
	pt := ptimer.NewFunc(func(i int) string {
		return fmt.Sprintf("%d (%dc %df)", i, len(clusters), friends/i)
	})
	done := make([]bool, d.Len())
	var ignored []int // Keeping ignored IDs for verification.
	for _, i := range perm {
		if done[i] {
			pt.Inc()
			continue
		}
		if d.IsIgnored(i) {
			ignored = append(ignored, i)
			pt.Inc()
			continue
		}
		done[i] = true
		s := d.Sketch(i)

		// Create cluster.
		c := []int{i}
		for sr := range d.SearchSketch(s, contn) {
			friends++
			if done[sr.I] {
				continue
			}
			if sr.Similarity < minSim {
				continue
			}
			c = append(c, sr.I)
			done[sr.I] = true
		}
		clusters = append(clusters, c)
		pt.Inc()
	}
	pt.Done()

	// This is an assertion. Failure means a bug in the code.
	checkClusterAssignment(append(clusters, ignored), d.Len())

	if secondAssn {
		fmt.Fprintln(os.Stderr, "Improving assignments")
		clusters = d.improveAssignments(clusters, contn)
		checkClusterAssignment(append(clusters, ignored), d.Len())
	}

	// Sort clusters for deterministic output.
	for _, c := range clusters {
		slices.Sort(c[1:]) // First element is the representative.
	}
	slices.SortFunc(clusters, func(a, b []int) int {
		return cmp.Compare(a[0], b[0])
	})

	return clusters
}

// Returns the indexes of slice elements if they were sorted.
func sortedPerm(n int, cmp func(int, int) int) []int {
	return snm.SortedFunc(
		snm.Slice(n, func(i int) int { return i }), cmp)
}

// Returns a sorted permutation according to each sketch's babi score.
func (d *Dataset[T]) babiSort() []int {
	scores := d.babiScores()
	return sortedPerm(d.Len(), func(a, b int) int {
		return cmp.Compare(scores[b], scores[a])
	})
}

// Returns the babi score of each element, by the order of iteration.
func (d *Dataset[T]) babiScores() []int {
	t := time.Now()
	cnt := map[T]int{}
	n := 0
	for _, sketch := range d.sketches {
		n++
		for _, h := range sketch {
			cnt[h]++
		}
	}
	if babiPrints {
		fmt.Fprintln(os.Stderr, "Babi counting took", time.Since(t))
	}

	if babiIgnoreCommon {
		t := time.Now()
		before := len(cnt)
		thrsh := int(math.Round(math.Pow(float64(n), 0.9)))
		for k, v := range cnt {
			if v > thrsh {
				delete(cnt, k)
			}
		}
		if babiPrints {
			diff := before - len(cnt)
			fmt.Fprintf(os.Stderr,
				"Babi filtering took %v, threshold %v, removed %v/%v (%.1f%%)\n",
				time.Since(t), thrsh, diff, before,
				float64(diff)/float64(before)*100)
		}
	}

	t = time.Now()
	scores := make([]int, 0, n)
	for _, sketch := range d.sketches {
		s := 0
		for _, h := range sketch {
			s += cnt[h] - 1
		}
		scores = append(scores, s)
	}
	if babiPrints {
		fmt.Fprintln(os.Stderr, "Babi scoring took", time.Since(t))
	}
	return scores
}

// Checks that the clusters include all the numbers from 0 to n-1
// and with no repetitions.
func checkClusterAssignment(clusters [][]int, n int) {
	all := snm.Sorted(slices.Concat(clusters...))
	if len(all) != n {
		panic(fmt.Sprintf("bad number of elements: %v, want %v",
			len(all), n))
	}
	for i, x := range all {
		if x != i {
			panic(fmt.Sprintf("bad element: %v, want %v", x, i))
		}
	}
}

// Reassigns each non-representative element to its closest representative.
func (d *Dataset[T]) improveAssignments(clusters [][]int, contn bool) [][]int {
	// Keep only reps in the search index.
	d.reindex(snm.SliceToSlice(clusters, func(c []int) int { return c[0] }))

	m := map[int][]int{}
	original := map[int]int{}
	for _, c := range clusters {
		m[c[0]] = []int{}
		for _, x := range c {
			original[x] = c[0]
		}
	}
	pt := ptimer.New()
	for i := range d.Len() {
		pt.Inc()
		if d.IsIgnored(i) {
			continue
		}
		// Skip representatives.
		if m[i] != nil {
			continue
		}

		best, besti := -1.0, original[i]
		for c := range d.SearchSketch(d.Sketch(i), contn) {
			if c.Similarity > best {
				best = c.Similarity
				besti = c.I
			}
		}
		m[besti] = append(m[besti], i)
	}
	pt.Done()
	var newClusters [][]int
	for k, v := range m {
		newClusters = append(newClusters, append([]int{k}, v...))
	}
	return newClusters
}
