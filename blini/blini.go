package main

import (
	"flag"
	"fmt"
	"os"
	"reflect"
	"runtime/debug"
)

/*
TODO
- Tests for sketching and for common
- Grouping by file, by regex?
*/

var (
	qFile     = flag.String("q", "", "Query file")
	rFile     = flag.String("r", "", "Reference file")
	oFile     = flag.String("o", "", "Output file or prefix")
	contn     = flag.Bool("c", false, "Use containment rather than full match")
	minSim    = flag.Float64("m", 0.9, "Minimum similarity for match")
	scale     = flag.Uint64("s", defaultScale, "Use 1/`scale` of the kmers")
	unmatched = flag.Bool("u", false, "Include unmatched queries in search output")

	dump        bool // An optional flag.
	ignoreShort bool // An optional flag.

	version = "development version"
)

func main() {
	flag.Parse()
	debug.SetGCPercent(20)

	if experimental {
		fmt.Fprintln(os.Stderr, "Hash type:", reflect.TypeFor[hashType]())
	}

	var err error
	if *qFile != "" && *rFile != "" {
		err = mainSearch()
	} else if (withDump && dump) && *qFile != "" {
		err = mainDump()
	} else if *qFile != "" {
		err = mainCluster()
	} else if *rFile != "" {
		err = mainSketch()
	} else {
		fmt.Fprintf(os.Stderr, "Blini (%s)\n\n", version)
		fmt.Fprintln(os.Stderr, "Please select -q for clustering, -r for sketching,",
			"or both for searching.")
		flag.PrintDefaults()
		os.Exit(1)
	}
	if err != nil {
		fmt.Fprintln(os.Stderr, "ERROR:", err)
		os.Exit(2)
	}
}

func init() {
	if withDump {
		flag.BoolVar(&dump, "d", false,
			"With -q, instead of clustering print pairwise distances into the output file")
	}
	if withIgnoreShort {
		flag.BoolVar(&ignoreShort, "ignore-too-short", false,
			"Ignore sequences that are too short, rather than err")
	}
}
