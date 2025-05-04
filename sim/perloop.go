package sim

import (
	"fmt"
	"os"

	"golooper/config"
)

func SimulationA(config *config.Config) error {
	infile, err := os.Open(config.InfileName)
	if err != nil {
		return fmt.Errorf("error opening input file: %v", err)
	}
	defer infile.Close()

	// Create all output files using FileOps
	outFiles, err := CreateOutputFiles(config) // TODO: move fileops use to first place results are written
	if err != nil {
		return fmt.Errorf("error creating output files: %v", err)
	}
	defer outFiles.Close()

	// Test reading first few bytes of input file to validate it's accessible
	testBuf := make([]byte, 100)
	n, err := infile.Read(testBuf)
	if err != nil {
		return fmt.Errorf("error reading from input file: %v", err)
	}
	if n == 0 {
		return fmt.Errorf("input file is empty")
	}

	// Reset file pointer to beginning after test read
	_, err = infile.Seek(0, 0)
	if err != nil {
		return fmt.Errorf("error resetting file position: %v", err)
	}

	// Test writing to output files
	testStr := "Test output\n"
	if _, err := outFiles.BasePairProbWig.WriteString(testStr); err != nil {
		return fmt.Errorf("error writing to base pair prob wig: %v", err)
	}
	if _, err := outFiles.AverageEnergyWig.WriteString(testStr); err != nil {
		return fmt.Errorf("error writing to average energy wig: %v", err)
	}
	if _, err := outFiles.MinFreeEnergyWig.WriteString(testStr); err != nil {
		return fmt.Errorf("error writing to min free energy wig: %v", err)
	}
	if _, err := outFiles.BasePairProbBed.WriteString(testStr); err != nil {
		return fmt.Errorf("error writing to base pair prob bed: %v", err)
	}
	if _, err := outFiles.MinFreeEnergyBed.WriteString(testStr); err != nil {
		return fmt.Errorf("error writing to min free energy bed: %v", err)
	}
	if _, err := outFiles.ExtendedBasePairProbWig.WriteString(testStr); err != nil {
		return fmt.Errorf("error writing to extended base pair prob wig: %v", err)
	}

	// TODO: Add the rest of the simulation logic here

	return nil
}
