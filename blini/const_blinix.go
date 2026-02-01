//go:build blinix

package main

const (
	// Is this an experimental build.
	experimental = true

	// Default value for the scale flag.
	defaultScale = 40

	// Is the distance dump feature enabled.
	withDump = true
)

// Size of hashes used here.
type hashType = uint32
