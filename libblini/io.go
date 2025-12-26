package libblini

import (
	"fmt"
	"io"
	"iter"
	"strings"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/blini/sketching"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/bnry"
	"github.com/fluhus/gostuff/ptimer"
)

type sketch struct {
	h     []uint64 // Sketch hashes.
	ln    int      // Sequence length.
	name  string   // Sequence name.
	scale uint64   // Kmer selection scale.
}

func ReadDataset(file string, scale uint64) (*Dataset, error) {
	var d *Dataset
	var err error
	if strings.HasSuffix(file, indexSuffix) {
		fmt.Println("Reading prepared sketches")
		d, err = collectSketches(readSketches(file))
	} else {
		fmt.Println("Sketching reference sequences")
		d, err = collectSketches(sketchFile(file, scale))
	}
	if err != nil {
		return nil, err
	}
	fmt.Println("Indexing")
	d.index()
	return d, nil
}

func CreateSketchFile(inFile, outFile string, scale uint64) error {
	var out io.Writer
	if outFile == "" {
		fmt.Println("No output")
		out = io.Discard
	} else {
		if !strings.HasSuffix(outFile, indexSuffix) {
			outFile += indexSuffix
		}
		fmt.Println("Saving to:", outFile)
		f, err := aio.Create(outFile)
		if err != nil {
			return err
		}
		defer f.Close()
		out = f
	}

	fmt.Println("Sketching sequences")
	pt := ptimer.New()
	for s, err := range sketchFile(inFile, scale) {
		if err != nil {
			return err
		}
		if err := bnry.Write(out, s.h, s.ln, s.name, s.scale); err != nil {
			return err
		}
		pt.Inc()
	}
	pt.Done()
	return nil
}

// Sketches an input fasta file and iterates over the sketches.
func sketchFile(file string, scale uint64) iter.Seq2[sketch, error] {
	return func(yield func(sketch, error) bool) {
		for fa, err := range fasta.File(file) {
			if err != nil {
				yield(sketch{}, err)
				return
			}
			var e sketch
			e.h = sketching.Sketch(fa.Sequence, kmerLen, scale)
			e.ln = len(fa.Sequence)
			e.name = string(fa.Name)
			e.scale = scale
			if !yield(e, nil) {
				return
			}
		}
	}
}

// Iterates over sketches in a file.
func readSketches(file string) iter.Seq2[sketch, error] {
	return func(yield func(sketch, error) bool) {
		f, err := aio.Open(file)
		if err != nil {
			yield(sketch{}, err)
			return
		}
		defer f.Close()

		for {
			var e sketch
			err := bnry.Read(f, &e.h, &e.ln, &e.name, &e.scale)
			if err != nil {
				if err == io.EOF {
					return
				}
				yield(sketch{}, err)
				return
			}
			if !yield(e, nil) {
				return
			}
		}
	}
}

// Collects sketches from an iterator,
// validating that their scales are the same.
func collectSketches(seq iter.Seq2[sketch, error]) (*Dataset, error) {
	skch := &Dataset{}
	first := true
	pt := ptimer.New()
	for s, err := range seq {
		if err != nil {
			return skch, err
		}
		if first {
			skch.scale = s.scale
			first = false
		} else {
			if s.scale != skch.scale {
				return skch, fmt.Errorf("mismatching scales: %d, %d",
					skch.scale, s.scale)
			}
		}
		skch.sketches = append(skch.sketches, s.h)
		skch.lens = append(skch.lens, s.ln)
		skch.names = append(skch.names, s.name)
		pt.Inc()
	}
	pt.Done()
	return skch, nil
}

// Builds the search index for this dataset.
func (d *Dataset) index() {
	d.idx = sketching.NewIndex(d.scale * idxScale)
	pt := ptimer.New()
	for i, s := range d.sketches {
		d.idx.Add(s, i)
		pt.Inc()
	}
	pt.Done()
}
