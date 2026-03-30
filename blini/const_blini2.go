//go:build !blini1 && !blinix

package main

const (
	experimental    = false // Is this an experimental build.
	defaultScale    = 40    // Default value for the scale flag.
	withDump        = true  // Is the distance dump feature enabled.
	withIgnoreShort = true  // Enable ignore too short sequences feature.
)

// Size of hashes used here.
type hashType = uint32
