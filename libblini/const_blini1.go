//go:build !blinix && !blini2

package libblini

const (
	minSketchSize = 0 // Sketches shorter than this will trigger an error.

	babiClustering   = false // Sort by babi scores rather than by length.
	babiPrints       = false // Print timing of babi stages.
	babiIgnoreCommon = true  // Ignore too common elements.
	secondAssn       = false // Run a second round of cluster assignments.
)
