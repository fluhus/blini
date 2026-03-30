//go:build blini1

package main

const (
	experimental    = false // Is this an experimental build.
	defaultScale    = 100   // Default value for the scale flag.
	withDump        = false // Enable distance dump feature.
	withIgnoreShort = false // Enable ignore too short sequences feature.
)

// Size of hashes used here.
type hashType = uint64
