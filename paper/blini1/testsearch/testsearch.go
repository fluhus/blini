// Tests search results of blini and sourmash.
package main

import (
	"fmt"
	"path/filepath"
	"regexp"

	"github.com/fluhus/gostuff/csvx"
	"github.com/fluhus/gostuff/sets"
	"golang.org/x/exp/maps"
)

const (
	outDir = "testdata/output"
)

func main() {
	e, err := readBlini(outDir + "/blini_vir.csv")
	if err != nil {
		panic(err)
	}
	fmt.Print("Blini: ")
	fmt.Println(findMatches(e))

	e, err = readBlini(outDir + "/blini_mut.csv")
	if err != nil {
		panic(err)
	}
	fmt.Print("Blini (mut): ")
	fmt.Println(findMatches(e))

	e, err = readSourmash(outDir + "/sm_vir_*.csv")
	if err != nil {
		panic(err)
	}
	fmt.Print("Sourmash: ")
	fmt.Println(findMatches(e))

	e, err = readSourmash(outDir + "/sm_mut_*.csv")
	if err != nil {
		panic(err)
	}
	fmt.Print("Sourmash (mut): ")
	fmt.Println(findMatches(e))

	e, err = readMMSeqs(outDir + "/mms")
	if err != nil {
		panic(err)
	}
	fmt.Print("MMseqs: ")
	fmt.Println(findMatches(e))

	e, err = readMMSeqs(outDir + "/mmsm")
	if err != nil {
		panic(err)
	}
	fmt.Print("MMseqs (mut): ")
	fmt.Println(findMatches(e))
}

// Reads matches produced by Blini.
func readBlini(file string) ([][2]string, error) {
	type entry struct {
		Query, Reference string
	}
	var result [][2]string
	for e, err := range csvx.DecodeFileHeader[entry](file) {
		if err != nil {
			return nil, err
		}
		result = append(result, [2]string{
			e.Query,
			firstWord(e.Reference)},
		)
	}
	return result, nil
}

// Reads matches produced by Sourmash.
func readSourmash(glob string) ([][2]string, error) {
	type entry struct {
		Query     string `csvdec:"query_name"`
		Reference string `csvdec:"name"`
	}
	var result [][2]string
	files, err := filepath.Glob(glob)
	if err != nil {
		return nil, err
	}
	for _, file := range files {
		for e, err := range csvx.DecodeFileHeader[entry](file) {
			if err != nil {
				return nil, err
			}
			result = append(result, [2]string{
				e.Query,
				firstWord(e.Reference)},
			)
		}
	}
	return result, nil
}

func readMMSeqs(file string) ([][2]string, error) {
	type entry struct {
		Query     string
		Reference string
	}
	result := sets.Set[[2]string]{}
	for e, err := range csvx.DecodeFile[entry](file, csvx.TSV) {
		if err != nil {
			return nil, err
		}
		result.Add([2]string{e.Query, e.Reference})
	}
	return maps.Keys(result), nil
}

// Returns the number of matches (query=reference)
// and the number of mismatches.
func findMatches(entries [][2]string) (int, int) {
	matches, mismatches := 0, 0
	for _, e := range entries {
		if e[0] == e[1] {
			matches++
		} else {
			mismatches++
		}
	}
	return matches, mismatches
}

var splitter = regexp.MustCompile(`^\S+`)

// Returns the first non-space substring.
func firstWord(s string) string {
	return splitter.FindString(s)
}
