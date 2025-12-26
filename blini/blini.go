package main

import (
	"flag"
	"fmt"
	"os"
	"runtime/debug"
)

/*
TODO
- Make ptimer output to stdout
- Tests for sketching and for common
- Grouping by file, by regex?
*/

const (
	unmatchedRef = "(unmatched)" // The "reference" value of an unmatched query.
)

var (
	qFile     = flag.String("q", "", "Query file")
	rFile     = flag.String("r", "", "Reference file")
	oFile     = flag.String("o", "", "Output file or prefix")
	contn     = flag.Bool("c", false, "Use containment rather than full match")
	minSim    = flag.Float64("m", 0.9, "Minimum similarity for match")
	scale     = flag.Uint64("s", 100, "Use 1/`scale` of the kmers")
	unmatched = flag.Bool("u", false, "Include unmatched queries in search output")

	version = "development version"
)

func main() {
	flag.Parse()
	debug.SetGCPercent(20)

	var err error
	if *qFile != "" && *rFile != "" {
		err = mainSearch()
	} else if *qFile != "" {
		err = mainCluster()
	} else if *rFile != "" {
		err = mainSketch()
	} else {
		fmt.Printf("Blini (%s)\n\n", version)
		fmt.Println("Please select -q for clustering, -r for sketching,",
			"or both for searching.")
		flag.PrintDefaults()
		os.Exit(1)
	}
	if err != nil {
		fmt.Println("ERROR:", err)
		os.Exit(2)
	}
}
