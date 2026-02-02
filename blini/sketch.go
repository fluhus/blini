// Sketching logic.

package main

import (
	"fmt"
	"os"

	"github.com/fluhus/blini/libblini"
)

// Main function for sketching operation.
func mainSketch() error {
	fmt.Fprintln(os.Stderr, "----------------")
	fmt.Fprintln(os.Stderr, "SKETCH OPERATION")
	fmt.Fprintln(os.Stderr, "----------------")
	fmt.Fprintln(os.Stderr, "Scale:", *scale)

	if *unmatched {
		return fmt.Errorf("flag -u is for search, not for sketching")
	}

	return libblini.CreateSketchFile[hashType](*rFile, *oFile, *scale, ignoreShort)
}
