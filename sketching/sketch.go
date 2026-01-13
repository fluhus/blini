package sketching

import (
	"bytes"
	"math"

	"github.com/fluhus/biostuff/sequtil"
	"github.com/fluhus/gostuff/hashx"
	"github.com/fluhus/gostuff/sets"
	"github.com/fluhus/gostuff/snm"
	"golang.org/x/exp/constraints"
	"golang.org/x/exp/maps"
)

// Sketch returns a sketch with 1/scale kmer hashes/
func Sketch[T constraints.Unsigned](seq []byte, k int, scale uint64) []T {
	seq = bytes.ToUpper(seq)
	hashes := make(sets.Set[T], len(seq)/int(scale))
	mx := math.MaxUint64 / scale
	for sseq := range sequtil.SubsequencesWith(seq, "atcgATCG") {
		for s := range sequtil.CanonicalSubsequences(sseq, k) {
			h := hashx.Bytes(s)
			if h > mx {
				continue
			}
			hashes.Add(T(h))
		}
	}
	return snm.Sorted(maps.Keys(hashes))
}
