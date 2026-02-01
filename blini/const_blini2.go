//go:build blini2

package main

const (
	experimental = false // Is this an experimental build.
	defaultScale = 40    // Default value for the scale flag.
	withDump     = true  // Is the distance dump feature enabled.

	babiClustering   = true  // Sort by babi scores rather than by length.
	babiPrints       = false // Print timing of babi stages.
	babiIgnoreCommon = true  // Ignore too common elements.
	secondAssn       = true  // Run a second round of cluster assignments.
)

// Size of hashes used here.
type hashType = uint32
