//go:build !blini2 && !blinix

package sketching

const (
	// Require at least 1/x of hashes to match for an index hit.
	minIndexHitRatio = 1

	// Which type of index to use.
	indexType = useSVMapIndex
)
