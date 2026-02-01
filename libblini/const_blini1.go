//go:build !blinix && !blini2

package libblini

const (
	minSketchSize = 0 // Sketches shorter than this will trigger an error.
)
