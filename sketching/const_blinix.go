//go:build blinix

package sketching

const (
	// Require at least 1/x of hashes to match for an index hit.
	minIndexHitRatio = 100

	// Which type of index to use.
	indexType = useFlatMap3Index
)
