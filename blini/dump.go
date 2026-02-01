package main

import (
	"fmt"

	"github.com/fluhus/blini/libblini"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/ptimer"
)

// Main function for dump operation.
func mainDump() error {
	fmt.Println("-----------------------")
	fmt.Println("DISTANCE DUMP OPERATION")
	fmt.Println("-----------------------")
	fmt.Println("Scale:", *scale)

	if *unmatched {
		return fmt.Errorf("flag -u is for search, not for dump")
	}
	if *oFile == "" {
		return fmt.Errorf("no output file given")
	}

	db, err := libblini.ReadDataset[hashType](*qFile, *scale, false, ignoreShort)
	if err != nil {
		return err
	}

	fmt.Println("Dumping distances")
	fout, err := aio.Create(*oFile)
	if err != nil {
		return err
	}
	defer fout.Close()

	n := db.Len()
	nd := n * (n - 1) / 2
	pt := ptimer.NewMessage(fmt.Sprint("{}/", nd))
	for i := range n {
		for jj := range n - i - 1 {
			j := jj + i + 1
			d := 1 - db.Sketch(i).Similarity(db.Sketch(j), false)
			fmt.Fprintln(fout, d)
			pt.Inc()
		}
	}
	pt.Done()

	return nil
}
