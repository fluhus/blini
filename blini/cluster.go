// Clustering logic.

package main

import (
	"fmt"
	"os"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/blini/libblini"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/snm"
)

// Main function for clustering operation.
func mainCluster() error {
	fmt.Fprintln(os.Stderr, "--------------------")
	fmt.Fprintln(os.Stderr, "CLUSTERING OPERATION")
	fmt.Fprintln(os.Stderr, "--------------------")
	fmt.Fprintln(os.Stderr, "Scale:", *scale)
	fmt.Fprintln(os.Stderr, "Min sim:", *minSim)

	if *unmatched {
		return fmt.Errorf("flag -u is for search, not for clustering")
	}

	db, err := libblini.ReadDataset[hashType](*qFile, *scale, true, ignoreShort)
	if err != nil {
		return err
	}

	fmt.Fprintln(os.Stderr, "Clustering")
	clusters := db.Cluster(*minSim, *contn)

	if *oFile != "" {
		fmt.Fprintln(os.Stderr, "Generating output")

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
		fmt.Fprintln(os.Stderr, "Cluster assignment:    ", *oFile+".json")
		fmt.Fprintln(os.Stderr, "Representative genomes:", *oFile+".fasta")
	} else {
		fmt.Fprintln(os.Stderr, "No output")
	}

	return nil
}
