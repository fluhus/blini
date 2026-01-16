//go:build blexp

package main

const (
	// Is this an experimental build.
	experimental = true

	// Default value for the scale flag.
	defaultScale = 40
)

// Size of hashes used here.
type hashType = uint32
