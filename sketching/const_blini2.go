//go:build !blini1 && !blinix

package sketching

const (
	// Require at least 1/x of hashes to match for an index hit.
	minIndexHitRatio = 100
)
