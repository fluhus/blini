//go:build blinix

package main

const (
	experimental    = true // Is this an experimental build.
	defaultScale    = 40   // Default value for the scale flag.
	withDump        = true // Is the distance dump feature enabled.
	withIgnoreShort = true // Enable ignore too short sequences feature.

	babiClustering   = true // Sort by babi scores rather than by length.
	babiPrints       = true // Print timing of babi stages.
	babiIgnoreCommon = true // Ignore too common elements.
	secondAssn       = true // Run a second round of cluster assignments.
)

// Size of hashes used here.
type hashType = uint32
