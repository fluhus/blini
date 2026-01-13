// Clustering logic.

package main

import (
	"cmp"
	"fmt"
	"maps"
	"slices"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/blini/libblini"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/ptimer"
	"github.com/fluhus/gostuff/sets"
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

	db, err := libblini.ReadDataset[hashType](*qFile, *scale)
	if err != nil {
		return err
	}
	db.CleanIndex()

	fmt.Println("Clustering")
	perm := sortedPerm(db.Len(), func(a, b int) int {
		return cmp.Compare(db.SequenceLen(b), db.SequenceLen(a))
	})
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
		for sr := range db.SearchSketch(s, db.SequenceLen(i), *contn) {
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

	{ // TODO(fluhus): Organize this code.
		alli := sets.Set[int]{}
		lens := 0
		for _, c := range clusters {
			alli.Add(c...)
			lens += len(c)
		}
		if lens != len(alli) {
			fmt.Println("lens!=alli:", lens, len(alli))
		}
		if !maps.Equal(alli, sets.Of(perm...)) {
			fmt.Println("bad alli")
		}
	}

	// Sort clusters for deterministic output.
	for _, c := range clusters {
		slices.Sort(c[1:]) // First element is the representative.
	}
	slices.SortFunc(clusters, func(a, b []int) int {
		return cmp.Compare(a[0], b[0])
	})

	// Create clusters by names.
	byName := snm.SliceToSlice(clusters, func(c []int) []string {
		return snm.SliceToSlice(c, func(i int) string {
			return db.Name(i)
		})
	})

	if *oFile != "" {
		fmt.Println("Generating output")
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
