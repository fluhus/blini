// Sketching logic.

package main

import (
	"fmt"

	"github.com/fluhus/blini/libblini"
)

// Main function for sketching operation.
func mainSketch() error {
	fmt.Println("----------------")
	fmt.Println("SKETCH OPERATION")
	fmt.Println("----------------")
	fmt.Println("Scale:", *scale)

	if *unmatched {
		return fmt.Errorf("flag -u is for search, not for sketching")
	}

	return libblini.CreateSketchFile[hashType](*rFile, *oFile, *scale, ignoreShort)
}
