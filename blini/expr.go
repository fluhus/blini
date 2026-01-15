//go:build blexp

package main

const (
	// Is this an experimental build.
	expr = true

	defaultScale = 40
)

// Size of hashes used here.
type hashType = uint32
