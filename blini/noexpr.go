//go:build !blexp

package main

const (
	// Is this an experimental build.
	experimental = false

	// Default value for the scale flag.
	defaultScale = 100
)

// Size of hashes used here.
type hashType = uint64
