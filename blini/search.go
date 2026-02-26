// Search logic.

package main

import (
	"encoding/csv"
	"fmt"
	"io"
	"os"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/blini/libblini"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/ptimer"
)

const (
	unmatchedRef = "(unmatched)" // The "reference" value of an unmatched query.
)

// Main function for search operation.
func mainSearch() error {
	fmt.Fprintln(os.Stderr, "----------------")
	fmt.Fprintln(os.Stderr, "SEARCH OPERATION")
	fmt.Fprintln(os.Stderr, "----------------")
	db, err := libblini.ReadDataset[hashType](*rFile, *scale, true, ignoreShort)
	if err != nil {
		return err
	}

	fmt.Fprintln(os.Stderr, "Searching")
	var fout io.Writer
	if *oFile != "" {
		f, err := aio.Create(*oFile)
		if err != nil {
			return err
		}
		defer f.Close()
		fout = f
	} else {
		fout = io.Discard
	}
	out := csv.NewWriter(fout)
	defer out.Flush()

	out.Write([]string{"similarity", "query", "reference"})

	var matches int
	pt := ptimer.NewFunc(func(i int) string {
		return fmt.Sprintf("%d (%d matches)", i, matches)
	})
	for fa, err := range fasta.File(*qFile) {
		if err != nil {
			return err
		}
		found := false
		for sr := range db.Search(fa.Sequence, *contn) {
			if sr.Similarity >= *minSim {
				matches++
				output := []string{
					fmt.Sprintf("%.0f%%", sr.Similarity*100),
					string(fa.Name),
					sr.Name,
				}
				out.Write(output)
				found = true
			}
		}
		if !found && *unmatched { // Report unmatched query.
			output := []string{"0%", string(fa.Name), unmatchedRef}
			out.Write(output)
		}
		pt.Inc()
	}
	pt.Done()

	return nil
}
