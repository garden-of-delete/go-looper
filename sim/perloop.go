package sim

import (
	"fmt"
	"os"

	"golooper/cmd"
)

func SimulationA(config *cmd.Config) error {

	// Open output files
	basePairProbWig, err := os.Create(config.OutfileName + "_bpprob.wig")
	if err != nil {
		return fmt.Errorf("error creating bpprob wig file: %v", err)
	}
	defer basePairProbWig.Close()

	averageEnergyWig, err := os.Create(config.OutfileName + "_avgG.wig")
	if err != nil {
		return fmt.Errorf("error creating avgG wig file: %v", err)
	}
	defer averageEnergyWig.Close()

	minFreeEnergyWig, err := os.Create(config.OutfileName + "_mfe.wig")
	if err != nil {
		return fmt.Errorf("error creating mfe wig file: %v", err)
	}
	defer minFreeEnergyWig.Close()

	basePairProbBed, err := os.Create(config.OutfileName + "_bpprob.bed")
	if err != nil {
		return fmt.Errorf("error creating bpprob bed file: %v", err)
	}
	defer basePairProbBed.Close()

	minFreeEnergyBed, err := os.Create(config.OutfileName + "_mfe.bed")
	if err != nil {
		return fmt.Errorf("error creating mfe bed file: %v", err)
	}
	defer minFreeEnergyBed.Close()

	extendedBasePairProbWig, err := os.Create(config.OutfileName + "_extbpprob.wig")
	if err != nil {
		return fmt.Errorf("error creating extbpprob wig file: %v", err)
	}
	defer extendedBasePairProbWig.Close()

	// Write headers to wig files
	if err := writeWigfileHeader(basePairProbWig, "Base pair probability"); err != nil {
		return fmt.Errorf("error writing bpprob wig header: %v", err)
	}

	if err := writeWigfileHeader(averageEnergyWig, "Average free energy"); err != nil {
		return fmt.Errorf("error writing avgG wig header: %v", err)
	}

	if err := writeWigfileHeader(minFreeEnergyWig, "Minimum free energy"); err != nil {
		return fmt.Errorf("error writing mfe wig header: %v", err)
	}

	if err := writeWigfileHeader(extendedBasePairProbWig, "Extended base pair probability"); err != nil {
		return fmt.Errorf("error writing extbpprob wig header: %v", err)
	}

	// Write headers to bed files
	if err := writeBedfileHeader(basePairProbBed, "Base pair probability"); err != nil {
		return fmt.Errorf("error writing bpprob bed header: %v", err)
	}

	if err := writeBedfileHeader(minFreeEnergyBed, "Minimum free energy"); err != nil {
		return fmt.Errorf("error writing mfe bed header: %v", err)
	}

	return nil
}
