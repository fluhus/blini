package main

import (
	"cmp"
	"flag"
	"fmt"
	"os"
	"os/signal"
	"regexp"
	"runtime"
	"runtime/debug"
	"slices"
	"strconv"
	"sync/atomic"
	"time"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/biostuff/mash"
	"github.com/fluhus/blini/fastani"
	"github.com/fluhus/blini/frindex"
	"github.com/fluhus/gostuff/gnum"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/minhash"
	"github.com/fluhus/gostuff/ptimer"
	"github.com/fluhus/gostuff/sets"
	"github.com/fluhus/gostuff/snm"
	"golang.org/x/exp/maps"
)

const (
	sslen       = 10000 // Length of subsequences for min-hashing.
	ssHalfSteps = false // Whether to use x2 subsequences.
	short       = 0     // Stop after [short] sequences, for debugging.
)

var (
	inFile     = flag.String("i", "", "Input fasta")
	inRef      = flag.String("r", "", "Reference to search in")
	outFile    = flag.String("o", "", "Output file prefix")
	sketchSize = flag.Int("n", 20, "Sketch size")
	minFriends = flag.Int("m", 1, "Minimum common hashes for testing")
	minANIPerc = flag.Float64("p", 90, "Minimum percent identity for clustering")
	seed       = flag.Uint("s", 0, "Hash seed (default 0)")
	tolerate   = flag.Bool("t", false, "Tolerate 1 missing part")
)

func main() {
	if err := myLinClust(); err != nil {
		fmt.Fprintln(os.Stderr, "ERROR:", err)
		os.Exit(2)
	}
}

// Runs clustering on input sequences.
func myLinClust() error {
	flag.Parse()
	mash.Seed = uint32(*seed)
	debug.SetGCPercent(33)

	fmt.Println("Perc:", *minANIPerc)
	fmt.Println("Sketch size:", *sketchSize)
	fmt.Println("Min friends:", *minFriends)
	fmt.Println("Input:", *inFile)
	if *inRef != "" {
		fmt.Println("Reference:", *inRef)
	}
	fmt.Println("Output:", *outFile)
	if sslen > 0 {
		fmt.Println("Using subseqs:", sslen)
	}
	fmt.Println("Tolerate:", *tolerate)

	fmt.Print("Checking fastANI... ")
	if err := fastani.Check(); err != nil {
		return err
	}
	fmt.Println("OK")

	if *inRef != "" {
		return mySearch()
	}

	start := time.Now()

	fmt.Println("Reading sequences")
	pt := ptimer.New()
	var seqs []*fasta.Fasta
	for fa, err := range fasta.IterFile(*inFile) {
		if err != nil {
			return err
		}
		seqs = append(seqs, fa)
		pt.Inc()
		if short > 0 && len(seqs) >= short {
			break
		}
	}
	pt.Done()

	fmt.Println("Sorting")
	pt = ptimer.New()
	perm := snm.Slice(len(seqs), func(i int) int { return i })
	names := snm.SliceToSlice(seqs, func(f *fasta.Fasta) string {
		return string(f.Name)
	})
	slices.SortFunc(perm, func(a, b int) int {
		return cmp.Compare(len(seqs[b].Sequence), len(seqs[a].Sequence))
	})
	merp := flipPerm(perm)
	pt.Done()

	fmt.Println("Indexing")
	pt = ptimer.New()
	idx := frindex.New[uint64](*sketchSize, *minFriends)
	var mhss [][]*minhash.MinHash[uint64]
	var imh []int
	for i, fa := range seqs {
		if sslen == 0 {
			mh := mash.Sequences(*sketchSize, 21, fa.Sequence)
			mhss = append(mhss, []*minhash.MinHash[uint64]{mh})
			idx.Add(mh)
		} else {
			bla := ssmash(fa.Sequence, sslen)
			for _, mh := range bla {
				imh = append(imh, i)
				idx.Add(mh)
			}
			mhss = append(mhss, bla)
		}
		pt.Inc()
	}
	pt.Done()

	setInterruptHandler() // For cleaning up nicely in case of interruption.

	dir, err := os.MkdirTemp("", "blini-")
	if err != nil {
		return err
	}
	defer os.RemoveAll(dir)
	fmt.Println("Writing temp files into:", dir)
	pt = ptimer.New()
	for i, fa := range seqs {
		if interrupted.Load() {
			return fmt.Errorf("interrupted")
		}
		f := fmt.Sprintf("%s/%d.fa", dir, i)
		txt, _ := fa.MarshalText()
		if err := os.WriteFile(f, txt, 0o644); err != nil {
			return nil
		}
		pt.Inc()
	}
	pt.Done()

	runtime.GC()

	fmt.Println("Clustering")
	n, nn := 0, 0
	var clusters [][]int
	pt = ptimer.NewFunc(func(i int) string {
		return fmt.Sprintf("%d (%dc %.0ff)",
			i, len(clusters), float64(n)/float64(nn))
	})
	for _, i := range perm {
		if interrupted.Load() {
			return fmt.Errorf("interrupted")
		}
		if mhss[i] == nil {
			pt.Inc()
			continue
		}
		nn++
		c := []int{i}
		var friends []int
		if sslen == 0 {
			friends = idx.Query(mhss[i][0])
		} else {
			frs := sets.Set[int]{}
			for _, mh := range mhss[i] {
				for _, fr := range idx.Query(mh) {
					frs.Add(imh[fr])
				}
			}
			friends = maps.Keys(frs)
		}
		friends = snm.FilterSlice(friends, func(j int) bool {
			return merp[j] > merp[i] && mhss[j] != nil // Keep only shorter sequences.
		})
		n += len(friends)
		if len(friends) == 0 {
			clusters = append(clusters, c)
			pt.Inc()
			continue
		}
		r := []string{fmt.Sprintf("%s/%d.fa", dir, i)}
		q := snm.SliceToSlice(friends, func(i int) string {
			return fmt.Sprintf("%s/%d.fa", dir, i)
		})
		for e, err := range fastani.CompareFiles(q, r) {
			if err != nil {
				return err
			}
			p := float64(e.PercID) * float64(e.PartsAligned) / float64(e.PartsTotal)
			if *tolerate {
				p = float64(e.PercID) * float64(e.PartsAligned) / float64(max(e.PartsTotal-1, e.PartsAligned))
			}
			if p >= *minANIPerc {
				n := faNumber(e.Query)
				if n == -1 {
					return fmt.Errorf("bad query name: %s, %v", e.Query, e)
				}
				mhss[n] = nil
				c = append(c, n)
			}
		}
		clusters = append(clusters, c)
		pt.Inc()
	}
	pt.Done()
	if err := jio.Save(*outFile+".bynumber.json", clusters); err != nil {
		return err
	}
	byName := snm.SliceToSlice(clusters, func(a []int) []string {
		return snm.SliceToSlice(a, func(i int) string { return names[i] })
	})
	if err := jio.Save(*outFile+".byname.json", byName); err != nil {
		return err
	}

	fmt.Println("Took", time.Since(start))
	return nil
}

// Searches query sequences in a reference dataset.
func mySearch() error {
	start := time.Now()

	fmt.Println("Reading sequences")
	pt := ptimer.New()
	var seqs []*fasta.Fasta
	var rnames []string
	var lens []int
	for fa, err := range fasta.IterFile(*inRef) {
		if err != nil {
			return err
		}
		seqs = append(seqs, fa)
		rnames = append(rnames, string(fa.Name))
		lens = append(lens, len(fa.Sequence))
		pt.Inc()
	}
	pt.Done()

	fmt.Println("Indexing")
	pt = ptimer.New()
	idx := frindex.New[uint64](*sketchSize, *minFriends)
	var imh []int
	for i, fa := range seqs {
		if sslen == 0 {
			mh := mash.Sequences(*sketchSize, 21, fa.Sequence)
			idx.Add(mh)
		} else {
			bla := ssmash(fa.Sequence, sslen)
			for _, mh := range bla {
				imh = append(imh, i)
				idx.Add(mh)
			}
		}
		pt.Inc()
	}
	pt.Done()

	setInterruptHandler() // For cleaning up nicely in case of interruption.

	dir, err := os.MkdirTemp("", "blini-")
	if err != nil {
		return err
	}
	defer os.RemoveAll(dir)
	fmt.Println("Writing temp files into:", dir)
	pt = ptimer.New()
	for i, fa := range seqs {
		if interrupted.Load() {
			return fmt.Errorf("interrupted")
		}
		f := fmt.Sprintf("%s/%d.fa", dir, i)
		txt, _ := fa.MarshalText()
		if err := os.WriteFile(f, txt, 0o644); err != nil {
			return nil
		}
		pt.Inc()
	}
	pt.Done()

	runtime.GC()

	fmt.Println("Searching")
	n, nn := 0, 0
	matches := map[int]map[int]float64{}
	var qnames []string
	pt = ptimer.NewFunc(func(i int) string {
		return fmt.Sprintf("%d (%dm %.1ff)",
			i, len(matches), float64(n)/float64(nn))
	})
	for fa, err := range fasta.IterFile(*inFile) {
		if err != nil {
			return err
		}
		if interrupted.Load() {
			return fmt.Errorf("interrupted")
		}
		if short > 0 && pt.N >= short {
			break
		}
		nn++
		qnames = append(qnames, string(fa.Name))
		var friends []int
		if sslen == 0 {
			friends = idx.Query(mash.Sequences(*sketchSize, 21, fa.Sequence))
		} else {
			frs := sets.Set[int]{}
			for _, mh := range ssmash(fa.Sequence, sslen) {
				for _, fr := range idx.Query(mh) {
					frs.Add(imh[fr])
				}
			}
			friends = maps.Keys(frs)
		}
		if len(friends) == 0 {
			pt.Inc()
			continue
		}
		n += len(friends)
		var shorter, longer []int
		for _, fr := range friends {
			if lens[fr] < len(fa.Sequence) {
				shorter = append(shorter, fr)
			} else {
				longer = append(longer, fr)
			}
		}
		faf := fmt.Sprintf("%s/%d.q.fa", dir, pt.N)
		{
			txt, _ := fa.MarshalText()
			if err := os.WriteFile(faf, txt, 0o644); err != nil {
				return err
			}
		}
		c := map[int]float64{}
		if len(shorter) > 0 {
			r := []string{faf}
			q := snm.SliceToSlice(friends, func(i int) string {
				return fmt.Sprintf("%s/%d.fa", dir, i)
			})
			for e, err := range fastani.CompareFiles(q, r) {
				if err != nil {
					return err
				}
				p := float64(e.PercID) * float64(e.PartsAligned) / float64(e.PartsTotal)
				if p >= *minANIPerc {
					c[faNumber(e.Query)] = p
				}
			}
		}
		if len(longer) > 0 {
			q := []string{faf}
			r := snm.SliceToSlice(friends, func(i int) string {
				return fmt.Sprintf("%s/%d.fa", dir, i)
			})
			for e, err := range fastani.CompareFiles(q, r) {
				if err != nil {
					return err
				}
				p := float64(e.PercID) * float64(e.PartsAligned) / float64(e.PartsTotal)
				if *tolerate {
					p = float64(e.PercID) * float64(e.PartsAligned) / float64(max(e.PartsTotal-1, e.PartsAligned))
				}
				if p >= *minANIPerc {
					c[faNumber(e.Ref)] = p
				}
			}
		}
		if err := os.Remove(faf); err != nil {
			return err
		}
		if len(c) > 0 {
			matches[nn-1] = c
		}
		pt.Inc()
	}
	pt.Done()

	if err := jio.Save(*outFile+".bynumber.json", matches); err != nil {
		return err
	}
	byName := snm.MapToMap(matches, func(i int, m map[int]float64) (string, map[string]float64) {
		return qnames[i], snm.MapToMap(m, func(i int, f float64) (string, float64) {
			return rnames[i], f
		})
	})
	if err := jio.Save(*outFile+".byname.json", byName); err != nil {
		return err
	}
	fmt.Println("Took", time.Since(start))
	return nil
}

// Indicates that an interrupt signal was sent.
var interrupted atomic.Bool

// Connects interrupt signals to the indicator.
func setInterruptHandler() {
	c := make(chan os.Signal, 1)
	signal.Notify(c, os.Interrupt)
	go func() {
		<-c
		interrupted.Store(true)
		<-c
		signal.Reset(os.Interrupt)
	}()
}

// Extracts number from fasta name.
var faNumberRE = regexp.MustCompile(`\d+\.fa$`)

// Returns a fasta file's number, or -1 if failed to find.
func faNumber(s string) int {
	m := faNumberRE.FindString(s)
	if len(m) < 4 {
		return -1
	}
	i, err := strconv.Atoi(m[:len(m)-3])
	if err != nil {
		return -1
	}
	return i
}

// Min-hashes a sequence's parts.
func ssmash(b []byte, sslen int) []*minhash.MinHash[uint64] {
	n := len(b)
	nss := max(gnum.Idiv(n, sslen), 1)
	nss2 := nss
	if ssHalfSteps {
		nss *= 2
		nss2 = nss - 1
	}
	var bla []*minhash.MinHash[uint64]
	for i := range nss2 {
		start := i * n / nss
		end := (i + 1) * n / nss
		if ssHalfSteps {
			end = (i + 2) * n / nss
		}
		mh := mash.Sequences(*sketchSize, 21, b[start:end])
		bla = append(bla, mh)
	}
	return bla
}

func flipPerm(a []int) []int {
	b := make([]int, len(a))
	for i, v := range a {
		b[v] = i
	}
	return b
}
