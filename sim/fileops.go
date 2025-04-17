package sim

import (
	"fmt"
	"os"
)

type FileOps struct {
	Wigfile *os.File
	Bedfile *os.File
}

func writeWigfileHeader(outfile *os.File, trackname string) error {
	// Compose .wig header with track definition line
	header := fmt.Sprintf("track type=wiggle_0 name=\"%s\" visibility=full autoscale=off color=50,150,255 priority=10\n", trackname)

	// Write header to file
	_, err := outfile.WriteString(header)
	if err != nil {
		return fmt.Errorf("error writing wig header: %v", err)
	}

	return nil
}

func writeBedfileHeader(outfile *os.File, trackname string) error {
	// Compose bed header with track definition line
	header := fmt.Sprintf("track name=rLooper description=\"%s\" useScore=1\n", trackname)

	// Write header to file
	_, err := outfile.WriteString(header)
	if err != nil {
		return fmt.Errorf("error writing bed header: %v", err)
	}

	return nil
}