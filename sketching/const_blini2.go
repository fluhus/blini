//go:build blini2

package sketching

const (
	// Require at least 1/x of hashes to match for an index hit.
	minIndexHitRatio = 100

	// Use flat map rather than SV map for the index.
	useFlatMapIndex = false
)
