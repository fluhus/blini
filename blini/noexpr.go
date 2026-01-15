//go:build !blexp

package main

const (
	// Is this an experimental build.
	expr = false

	defaultScale = 100
)

// Size of hashes used here.
type hashType = uint64
