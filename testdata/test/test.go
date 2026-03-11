// Tests the outputs of running Blini on the test data.
package main

import (
	"cmp"
	"fmt"
	"os"
	"regexp"
	"slices"
	"strconv"
	"strings"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/gostuff/csvx"
	"github.com/fluhus/gostuff/gnum"
	"github.com/fluhus/gostuff/iterx"
	"github.com/fluhus/gostuff/jio"
	"golang.org/x/exp/maps"
)

const (
	dataDir    = "testdata"
	resultsDir = "tmp"
)

func main() {
	if err := testClusters(); err != nil {
		panic(err)
	}
	if err := testClustersSeqs(); err != nil {
		panic(err)
	}
	if err := testSearch(resultsDir+"/search.csv", false); err != nil {
		panic(err)
	}
	if err := testSearch(resultsDir+"/search_u.csv", true); err != nil {
		panic(err)
	}
	if err := testSearch(resultsDir+"/search2.csv", false); err != nil {
		panic(err)
	}
	if err := testSearch(resultsDir+"/search2_u.csv", true); err != nil {
		panic(err)
	}
	if err := testDump(); err != nil {
		panic(err)
	}
	fmt.Println("OK!")
}

// Tests that all sequences went to the right cluster.
func testClusters() error {
	var result struct {
		ByName [][]string
	}
	if err := jio.Read(resultsDir+"/clust.json", &result); err != nil {
		return err
	}
	got := result.ByName
	sortClusters(got)
	want := [][]string{
		{"ref1", "ref1.1", "ref1.2", "ref1.3"},
		{"ref2", "ref2.1", "ref2.2", "ref2.3"},
		{"ref3", "ref3.1", "ref3.2", "ref3.3"},
	}
	if !slices.EqualFunc(got, want, slices.Equal) {
		return fmt.Errorf("wrong clustering result: got %v, want %v",
			got, want)
	}
	return nil
}

// Tests that the output sequences are the same as input.
func testClustersSeqs() error {
	// Read original sequences.
	want := map[string]string{}
	for fa, err := range fasta.File(dataDir + "/clust.fa.zst") {
		if err != nil {
			return err
		}
		want[string(fa.Name)] = string(fa.Sequence)
	}

	// Test result sequences.
	namePattern := regexp.MustCompile(`^ref\d$`)
	for fa, err := range fasta.File(resultsDir + "/clust.fasta") {
		if err != nil {
			return err
		}
		if !namePattern.Match(fa.Name) {
			return fmt.Errorf("unexpected name: %q", fa.Name)
		}
		if string(fa.Sequence) != want[string(fa.Name)] {
			return fmt.Errorf("mismatching sequence for %q", fa.Name)
		}
	}

	return nil
}

// Tests that each query was matched with its reference.
func testSearch(file string, unmatched bool) error {
	const wantSimilarity = 98

	type entry struct {
		Similarity, Query, Reference string
	}

	queryRefMap := map[string]string{
		"query1.1": "ref1",
		"query1.2": "ref1",
		"query1.3": "ref1",
		"query2.1": "ref2",
		"query2.2": "ref2",
		"query2.3": "ref2",
		"query3.1": "ref3",
		"query3.2": "ref3",
		"query3.3": "ref3",
	}
	unmatchedPattern := regexp.MustCompile(`^query4\.\d$`)

	wantQueries := maps.Keys(queryRefMap)
	if unmatched {
		wantQueries = append(wantQueries, "query4.1", "query4.2", "query4.3")
	}
	var gotQueries []string

	similarities := map[string]bool{ // Accepted similarity values.
		fmt.Sprint(wantSimilarity, "%"):   true,
		fmt.Sprint(wantSimilarity+1, "%"): true,
		fmt.Sprint(wantSimilarity-1, "%"): true,
	}

	for row, err := range csvx.DecodeFileHeader[entry](file) {
		if err != nil {
			return err
		}
		gotQueries = append(gotQueries, row.Query)
		if unmatched {
			if unmatchedPattern.MatchString(row.Query) &&
				row.Reference == "(unmatched)" && row.Similarity == "0%" {
				continue
			}
		}
		if ref, ok := queryRefMap[row.Query]; !ok || row.Reference != ref {
			return fmt.Errorf("unexpected match: %q %q", row.Query, row.Reference)
		}
		if !similarities[row.Similarity] {
			return fmt.Errorf("unexpected similarity: %s", row.Similarity)
		}
	}

	slices.Sort(wantQueries)
	slices.Sort(gotQueries)
	if !slices.Equal(gotQueries, wantQueries) {
		return fmt.Errorf("unexpected queries: %q", gotQueries)
	}

	return nil
}

func testDump() error {
	const tolerance = 0.0021
	const file = resultsDir + "/dump.txt"

	if _, err := os.Stat(file); err != nil {
		// File missing is ok for blini1.
		if strings.HasSuffix(err.Error(), "no such file or directory") {
			fmt.Println("Dump file missing. Continuing...")
			return nil
		}
		return err
	}

	want := []float64{0.02, 0.02, 0.04, 0.04, 0.02, 0.06}
	var got []float64
	for line, err := range iterx.LinesFile(file) {
		if err != nil {
			return err
		}
		x, err := strconv.ParseFloat(line, 64)
		if err != nil {
			return err
		}
		got = append(got, x)
	}
	if len(got) != len(want) {
		return fmt.Errorf("bad result length: %d, want %v", len(got), len(want))
	}
	for i := range got {
		if gnum.Diff(got[i], want[i]) > tolerance {
			return fmt.Errorf("got[%d]=%f, want %f", i, got[i], want[i])
		}
	}
	return nil
}

// Sorts a clustering assignment to make clusterings comparable.
func sortClusters[T cmp.Ordered](a [][]T) {
	for _, x := range a {
		// Skip first element to keep cluster representative first.
		slices.Sort(x[1:])
	}
	slices.SortFunc(a, func(x, y []T) int {
		return cmp.Compare(x[0], y[0])
	})
}
