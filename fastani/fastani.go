package fastani

import (
	"bytes"
	"fmt"
	"iter"
	"os"
	"os/exec"
	"path/filepath"
	"strings"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/csvdec/v2"
)

const (
	fastani = "fastANI"
)

// CompareSeqs runs FastANI on the given fastas.
func CompareSeqs(query []*fasta.Fasta, ref []*fasta.Fasta) iter.Seq2[Entry, error] {
	return func(yield func(Entry, error) bool) {
		dir, err := os.MkdirTemp("/dev/shm", "fastani-")
		if err != nil {
			yield(Entry{}, err)
			return
		}
		qfile := filepath.Join(dir, "q")
		rfile := filepath.Join(dir, "r")
		ofile := filepath.Join(dir, "o")
		defer os.RemoveAll(dir)

		if err := writeFas(query, qfile); err != nil {
			yield(Entry{}, err)
			return
		}
		if err := writeFas(ref, rfile); err != nil {
			yield(Entry{}, err)
			return
		}

		cmd := exec.Command(fastani, "-q", qfile, "-r", rfile, "-o", ofile)
		out := bytes.NewBuffer(nil)
		cmd.Stdout = out
		cmd.Stderr = out
		if err := cmd.Run(); err != nil {
			yield(Entry{}, fmt.Errorf("%v %v %s", err, cmd.Args, out.Bytes()))
			return
		}

		for e, err := range readANI(ofile) {
			if !yield(e, err) {
				break
			}
		}
	}
}

// CompareFiles runs FastANI on the given fasta file names.
func CompareFiles(query []string, ref []string) iter.Seq2[Entry, error] {
	return func(yield func(Entry, error) bool) {
		dir, err := os.MkdirTemp("/dev/shm", "fastani-")
		if err != nil {
			yield(Entry{}, err)
			return
		}
		qfile := filepath.Join(dir, "q")
		rfile := filepath.Join(dir, "r")
		ofile := filepath.Join(dir, "o")
		defer os.RemoveAll(dir)

		qdata := strings.Join(query, "\n")
		rdata := strings.Join(ref, "\n")
		if err := os.WriteFile(qfile, []byte(qdata), 0o644); err != nil {
			yield(Entry{}, err)
			return
		}
		if err := os.WriteFile(rfile, []byte(rdata), 0o644); err != nil {
			yield(Entry{}, err)
			return
		}

		cmd := exec.Command(fastani, "--ql", qfile, "--rl", rfile, "-o", ofile)
		out := bytes.NewBuffer(nil)
		cmd.Stdout = out
		cmd.Stderr = out
		if err := cmd.Run(); err != nil {
			yield(Entry{}, fmt.Errorf("%v %v %s", err, cmd.Args, out.Bytes()))
			return
		}

		for e, err := range readANI(ofile) {
			if !yield(e, err) {
				break
			}
		}
	}
}

// Reads FastANI entries from its output file.
func readANI(file string) iter.Seq2[Entry, error] {
	return func(yield func(Entry, error) bool) {
		f, err := aio.Open(file)
		if err != nil {
			yield(Entry{}, err)
			return
		}
		defer f.Close()
		r := csvdec.New(f)
		r.Comma = '\t'
		for e, err := range csvdec.Iter[Entry](r) {
			if !yield(e, err) {
				break
			}
		}
	}
}

// Writes fasta sequences to an output file.
func writeFas(fas []*fasta.Fasta, file string) error {
	f, err := aio.Create(file)
	if err != nil {
		return err
	}
	defer f.Close()
	for _, fa := range fas {
		txt, _ := fa.MarshalText()
		if _, err := f.Write(txt); err != nil {
			return err
		}
	}
	return nil
}

// Entry is a single FastANI comparison.
type Entry struct {
	Query        string
	Ref          string
	PercID       float64
	PartsAligned int
	PartsTotal   int
}

// Check checks that fastANI can be run successfully.
func Check() error {
	cmd := exec.Command(fastani, "-h")
	err := cmd.Run()
	if err != nil {
		return fmt.Errorf("could not run fastANI: %w", err)
	}
	return nil
}
