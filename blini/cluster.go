// Clustering logic.

package main

import (
	"cmp"
	"fmt"
	"iter"
	"math"
	"slices"
	"time"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/blini/libblini"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/ptimer"
	"github.com/fluhus/gostuff/snm"
)

// Main function for clustering operation.
func mainCluster() error {
	fmt.Println("--------------------")
	fmt.Println("CLUSTERING OPERATION")
	fmt.Println("--------------------")
	fmt.Println("Scale:", *scale)
	fmt.Println("Min sim:", *minSim)

	if *unmatched {
		return fmt.Errorf("flag -u is for search, not for clustering")
	}

	db, err := libblini.ReadDataset[hashType](*qFile, *scale, true)
	if err != nil {
		return err
	}
	db.CleanIndex()

	fmt.Println("Clustering")
	perm := sortedPerm(db.Len(), func(a, b int) int {
		return cmp.Compare(db.Sketch(b).Len(), db.Sketch(a).Len())
	})
	if babiClustering {
		perm = babiSort(db)
	}

	friends := 0
	var clusters [][]int
	pt := ptimer.NewFunc(func(i int) string {
		return fmt.Sprintf("%d (%dc %df)", i, len(clusters), friends/i)
	})
	done := make([]bool, db.Len())
	for _, i := range perm {
		if done[i] {
			pt.Inc()
			continue
		}
		done[i] = true
		s := db.Sketch(i)

		// Create cluster.
		c := []int{i}
		for sr := range db.SearchSketch(s, *contn) {
			if done[sr.I] {
				continue
			}
			if sr.Similarity < *minSim {
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
	checkClusterAssignment(clusters, db.Len())

	if secondAssn {
		fmt.Println("Improving assignments")
		clusters = improveAssignments(db, clusters)
		checkClusterAssignment(clusters, db.Len())
	}

	// Sort clusters for deterministic output.
	for _, c := range clusters {
		slices.Sort(c[1:]) // First element is the representative.
	}
	slices.SortFunc(clusters, func(a, b []int) int {
		return cmp.Compare(a[0], b[0])
	})

	if *oFile != "" {
		fmt.Println("Generating output")

		// Create clusters by names.
		byName := snm.SliceToSlice(clusters, func(c []int) []string {
			return snm.SliceToSlice(c, func(i int) string {
				return db.Sketch(i).Name()
			})
		})

		// JSON output.
		output := map[string]any{
			"byNumber": clusters,
			"byName":   byName,
		}
		if err := jio.Write(*oFile+".json", output); err != nil {
			return err
		}

		// Fasta output.
		reps := snm.SliceToSlice(clusters, func(a []int) int { return a[0] })
		fout, err := aio.Create(*oFile + ".fasta")
		if err != nil {
			return err
		}
		defer fout.Close()
		i := -1
		for fa, err := range fasta.File(*qFile) {
			if err != nil {
				return err
			}
			if len(reps) == 0 {
				break
			}
			i++
			if i == reps[0] {
				if err := fa.Write(fout); err != nil {
					return err
				}
				reps = reps[1:]
			}
		}
		fmt.Println("Cluster assignment:    ", *oFile+".json")
		fmt.Println("Representative genomes:", *oFile+".fasta")
	} else {
		fmt.Println("No output")
	}

	return nil
}

// Returns the indexes of slice elements if they were sorted.
func sortedPerm(n int, cmp func(int, int) int) []int {
	return snm.SortedFunc(
		snm.Slice(n, func(i int) int { return i }), cmp)
}

// Returns a sorted permutation according to each sketch's babi score.
func babiSort(db *libblini.Dataset[hashType]) []int {
	scores := babiScores(func() iter.Seq[iter.Seq[hashType]] {
		return func(yield func(iter.Seq[hashType]) bool) {
			for i := range db.Len() {
				if !yield(db.Sketch(i).Hashes()) {
					break
				}
			}
		}
	})
	return sortedPerm(db.Len(), func(a, b int) int {
		return cmp.Compare(scores[b], scores[a])
	})
}

// Returns the babi score of each element, by the order of iteration.
func babiScores(f func() iter.Seq[iter.Seq[hashType]]) []int {
	t := time.Now()
	cnt := map[hashType]int{}
	n := 0
	for sketch := range f() {
		n++
		for h := range sketch {
			cnt[h]++
		}
	}
	if babiPrints {
		fmt.Println("Babi counting took", time.Since(t))
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
			fmt.Printf("Babi filtering took %v, threshold %v, removed %v/%v (%.1f%%)\n",
				time.Since(t), thrsh, diff, before,
				float64(diff)/float64(before)*100)
		}
	}

	t = time.Now()
	scores := make([]int, 0, n)
	for sketch := range f() {
		s := 0
		for h := range sketch {
			s += cnt[h] - 1
		}
		scores = append(scores, s)
	}
	if babiPrints {
		fmt.Println("Babi scoring took", time.Since(t))
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
func improveAssignments(db *libblini.Dataset[hashType], clusters [][]int) [][]int {
	// Keep only reps in the search index.
	db.Reindex(snm.SliceToSlice(clusters, func(c []int) int { return c[0] }))

	m := map[int][]int{}
	for _, c := range clusters {
		m[c[0]] = []int{}
	}
	pt := ptimer.New()
	for i := range db.Len() {
		pt.Inc()
		// Skip representatives.
		if m[i] != nil {
			continue
		}

		best, besti := -1.0, -1
		for c := range db.SearchSketch(db.Sketch(i), *contn) {
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
