//go:build blini2

package libblini

const (
	minSketchSize = 25 // Sketches shorter than this will trigger an error.

	babiClustering   = true  // Sort by babi scores rather than by length.
	babiPrints       = false // Print timing of babi stages.
	babiIgnoreCommon = true  // Ignore too common elements.
	secondAssn       = true  // Run a second round of cluster assignments.
)
