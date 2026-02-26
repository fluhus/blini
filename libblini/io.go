package libblini

import (
	"fmt"
	"io"
	"iter"
	"os"
	"strings"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/blini/sketching"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/bnry"
	"github.com/fluhus/gostuff/ptimer"
	"github.com/fluhus/gostuff/sets"
	"golang.org/x/exp/constraints"
)

// ReadDataset reads a sketch dataset from a file.
// Reads pre-sketched data if the file ends with .blini,
// otherwise treats the data as fasta and sketches it.
func ReadDataset[T constraints.Unsigned](
	file string, scale uint64, index bool, ignoreShort bool) (*Dataset[T], error) {
	var d *Dataset[T]
	var err error
	if strings.HasSuffix(file, indexSuffix) {
		fmt.Fprintln(os.Stderr, "Reading prepared sketches")
		d, err = collectSketches(readSketches[T](file), index, ignoreShort)
	} else {
		fmt.Fprintln(os.Stderr, "Sketching reference sequences")
		d, err = collectSketches(sketchFile[T](file, scale), index, ignoreShort)
	}
	if err != nil {
		return nil, err
	}
	return d, nil
}

// CreateSketchFile sketches a fasta file and outputs the sketches into a file.
func CreateSketchFile[T constraints.Unsigned](
	inFile, outFile string, scale uint64, ignoreShort bool) error {
	var out io.Writer
	if outFile == "" {
		fmt.Fprintln(os.Stderr, "No output")
		out = io.Discard
	} else {
		if !strings.HasSuffix(outFile, indexSuffix) {
			outFile += indexSuffix
		}
		fmt.Fprintln(os.Stderr, "Saving to:", outFile)
		f, err := aio.Create(outFile)
		if err != nil {
			return err
		}
		defer f.Close()
		out = f
	}

	fmt.Fprintln(os.Stderr, "Sketching sequences")
	ignored := 0
	pt := ptimer.New()
	if ignoreShort {
		pt = ptimer.NewFunc(func(i int) string {
			return fmt.Sprintf("%d (ignored %d)", i, ignored)
		})
	}
	for s, err := range sketchFile[T](inFile, scale) {
		if err != nil {
			return err
		}
		if len(s.h) < minSketchSize {
			if ignoreShort {
				ignored++
				pt.Inc()
				continue
			}
			return fmt.Errorf(
				"sketch has %v elements, but %v are needed for accurate results: %q",
				len(s.h), minSketchSize, s.name)
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
func sketchFile[T constraints.Unsigned](
	file string, scale uint64) iter.Seq2[Sketch[T], error] {
	return func(yield func(Sketch[T], error) bool) {
		for fa, err := range fasta.File(file) {
			if err != nil {
				yield(Sketch[T]{}, err)
				return
			}
			var e Sketch[T]
			e.h = sketching.Sketch[T](fa.Sequence, kmerLen, scale)
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
func readSketches[T constraints.Unsigned](file string) iter.Seq2[Sketch[T], error] {
	return func(yield func(Sketch[T], error) bool) {
		f, err := aio.Open(file)
		if err != nil {
			yield(Sketch[T]{}, err)
			return
		}
		defer f.Close()

		for {
			var e Sketch[T]
			err := bnry.Read(f, &e.h, &e.ln, &e.name, &e.scale)
			if err != nil {
				if err == io.EOF {
					return
				}
				yield(Sketch[T]{}, err)
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
func collectSketches[T constraints.Unsigned](
	seq iter.Seq2[Sketch[T], error], index bool, ignoreShort bool) (*Dataset[T], error) {
	d := &Dataset[T]{
		ignored: sets.Set[int]{},
	}
	if index {
		d.idx = sketching.NewIndex[T](idxScale)
	}
	pt := ptimer.New()
	if ignoreShort {
		pt = ptimer.NewFunc(func(i int) string {
			return fmt.Sprintf("%d (ignored %d)", i, len(d.ignored))
		})
	}
	for s, err := range seq {
		if err != nil {
			return d, err
		}
		if len(d.sketches) == 0 {
			d.scale = s.scale
		} else {
			if s.scale != d.scale {
				return d, fmt.Errorf("mismatching scales: %d, %d",
					d.scale, s.scale)
			}
		}
		if len(s.h) < minSketchSize {
			if ignoreShort {
				// Add an empty placeholder to keep indexes the same.
				d.sketches = append(d.sketches, nil)
				d.lens = append(d.lens, 0)
				d.names = append(d.names, "")
				d.ignored.Add(len(d.sketches) - 1)
				pt.Inc()
				continue
			}
			return nil, fmt.Errorf(
				"sketch has %v elements, but %v are needed for accurate results: %q",
				len(s.h), minSketchSize, s.name)
		}
		d.sketches = append(d.sketches, s.h)
		d.lens = append(d.lens, s.ln)
		d.names = append(d.names, s.name)
		if index {
			d.idx.Add(s.h, len(d.sketches)-1)
		}
		pt.Inc()
	}
	pt.Done()
	return d, nil
}
