// Tests clustering results of blini and mmseqs.
package main

import (
	"fmt"
	"path/filepath"
	"regexp"

	"github.com/fluhus/gostuff/clustering"
	"github.com/fluhus/gostuff/csvx"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/sets"
	"github.com/fluhus/gostuff/snm"
	"golang.org/x/exp/maps"
)

const (
	outDir = "testdata/output"
)

func main() {
	fmt.Println("Testing mmseqs clusters")
	groups, err := readMMSClusters(outDir + "/mm_cluster.tsv")
	if err != nil {
		panic(err)
	}
	fmt.Println(len(groups), "clusters")

	gotTags := toTags(groups)
	wantTags := toTags(toWantGroups(groups))
	fmt.Println("ARI:", clustering.AdjustedRandIndex(gotTags, wantTags))

	fmt.Println("Testing blini clusters")
	bliniFiles, _ := filepath.Glob(outDir + "/blini_*.json")
	for _, f := range bliniFiles {
		bln := struct {
			ByName [][]string
		}{}
		if err := jio.Read(f, &bln); err != nil {
			panic(err)
		}
		groups := bln.ByName

		gotTags := toTags(groups)
		wantTags := toTags(toWantGroups(groups))
		scale := f[10 : len(f)-5]
		fmt.Printf("[%s] ARI: %f (%d clusters)\n", scale,
			clustering.AdjustedRandIndex(gotTags, wantTags), len(groups))
	}
}

// Returns the clusters from MMseq's TSV output.
func readMMSClusters(file string) ([][]string, error) {
	pairs := sets.Set[[2]string]{}
	for line, err := range csvx.File(file, csvx.TSV) {
		if err != nil {
			return nil, err
		}
		pairs.Add([2]string{line[0], line[1]})
	}
	m := map[string][]string{}
	for pair := range pairs {
		m[pair[0]] = append(m[pair[0]], pair[1])
	}
	return maps.Values(m), nil
}

// Returns the cluster index number for each element.
func toTags(groups [][]string) []int {
	m := map[string]int{}
	for i, g := range groups {
		for _, s := range g {
			m[s] = i
		}
	}
	return snm.SliceToSlice(snm.Sorted(maps.Keys(m)), func(k string) int {
		return m[k]
	})
}

// Returns the "real" grouping.
func toWantGroups(groups [][]string) [][]string {
	re := regexp.MustCompile(`^(.*)\.\.\d+`)
	m := map[string][]string{}
	for _, g := range groups {
		for _, s := range g {
			k := re.FindStringSubmatch(s)[1]
			if k == "" {
				panic(fmt.Sprintf("could not match element: %q", s))
			}
			m[k] = append(m[k], s)
		}
	}
	return maps.Values(m)
}
