//go:build !blini2 && !blinix

package main

const (
	experimental = false // Is this an experimental build.
	defaultScale = 100   // Default value for the scale flag.
	withDump     = false // Is the distance dump feature enabled.

	babiClustering   = false // Sort by babi scores rather than by length.
	babiPrints       = true  // Print timing of babi stages.
	babiIgnoreCommon = true  // Ignore too common elements.
	secondAssn       = false // Run a second round of cluster assignments.
)

// Size of hashes used here.
type hashType = uint64
