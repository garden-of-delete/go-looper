package sim

import (
	"fmt"
	"golooper/config"
	"os"
	"path/filepath"
)

func createFileWithDir(path string) (*os.File, error) {
	dir := filepath.Dir(path)
	if err := os.MkdirAll(dir, 0755); err != nil {
		return nil, fmt.Errorf("error creating directory %s: %v", dir, err)
	}
	return os.Create(path)
}

func createOutputFile(path string, headerFunc func(*os.File) error, headerName string) (*os.File, error) {
	file, err := createFileWithDir(path)
	if err != nil {
		return nil, fmt.Errorf("error creating %s file at %s: %v", headerName, path, err)
	}
	if err := headerFunc(file); err != nil {
		file.Close()
		return nil, fmt.Errorf("error writing %s header: %v", headerName, err)
	}
	return file, nil
}

type FileOps struct {
	BasePairProbWig         *os.File
	AverageEnergyWig        *os.File
	MinFreeEnergyWig        *os.File
	BasePairProbBed         *os.File
	MinFreeEnergyBed        *os.File
	ExtendedBasePairProbWig *os.File
}

func writeWigfileHeader(outfile *os.File, trackname string) error {
	header := fmt.Sprintf("track type=wiggle_0 name=\"%s\" visibility=full autoscale=off color=50,150,255 priority=10\n", trackname)
	if _, err := outfile.WriteString(header); err != nil {
		return fmt.Errorf("error writing wig header: %v", err)
	}
	return nil
}

func writeBedfileHeader(outfile *os.File, trackname string) error {
	header := fmt.Sprintf("track name=rLooper description=\"%s\" useScore=1\n", trackname)
	if _, err := outfile.WriteString(header); err != nil {
		return fmt.Errorf("error writing bed header: %v", err)
	}
	return nil
}

func (f *FileOps) Close() error {
	var errs []error

	if f.BasePairProbWig != nil {
		if err := f.BasePairProbWig.Close(); err != nil {
			errs = append(errs, fmt.Errorf("error closing base pair prob wig: %v", err))
		}
	}
	if f.AverageEnergyWig != nil {
		if err := f.AverageEnergyWig.Close(); err != nil {
			errs = append(errs, fmt.Errorf("error closing average energy wig: %v", err))
		}
	}
	if f.MinFreeEnergyWig != nil {
		if err := f.MinFreeEnergyWig.Close(); err != nil {
			errs = append(errs, fmt.Errorf("error closing min free energy wig: %v", err))
		}
	}
	if f.BasePairProbBed != nil {
		if err := f.BasePairProbBed.Close(); err != nil {
			errs = append(errs, fmt.Errorf("error closing base pair prob bed: %v", err))
		}
	}
	if f.MinFreeEnergyBed != nil {
		if err := f.MinFreeEnergyBed.Close(); err != nil {
			errs = append(errs, fmt.Errorf("error closing min free energy bed: %v", err))
		}
	}
	if f.ExtendedBasePairProbWig != nil {
		if err := f.ExtendedBasePairProbWig.Close(); err != nil {
			errs = append(errs, fmt.Errorf("error closing extended base pair prob wig: %v", err))
		}
	}

	if len(errs) > 0 {
		return fmt.Errorf("errors closing files: %v", errs)
	}
	return nil
}

func CreateOutputFiles(config *config.Config) (*FileOps, error) {
	fileOps := &FileOps{}
	basePath := filepath.Join(filepath.Dir(config.OutfileName), filepath.Base(config.OutfileName))
	base := filepath.Base(config.OutfileName)

	cleanup := func() error {
		return fileOps.Close()
	}

	var err error

	fileOps.BasePairProbWig, err = createOutputFile(
		basePath+"_bpprob.wig",
		func(f *os.File) error { return writeWigfileHeader(f, "signal1_"+base) },
		"base pair probability wig",
	)
	if err != nil {
		cleanup()
		return nil, err
	}

	fileOps.AverageEnergyWig, err = createOutputFile(
		basePath+"_avgG.wig",
		func(f *os.File) error { return writeWigfileHeader(f, "signal2_"+base) },
		"average energy wig",
	)
	if err != nil {
		cleanup()
		return nil, err
	}

	fileOps.MinFreeEnergyWig, err = createOutputFile(
		basePath+"_mfe.wig",
		func(f *os.File) error { return writeWigfileHeader(f, "signal3_"+base) },
		"minimum free energy wig",
	)
	if err != nil {
		cleanup()
		return nil, err
	}

	fileOps.BasePairProbBed, err = createOutputFile(
		basePath+"_bpprob.bed",
		func(f *os.File) error { return writeBedfileHeader(f, "signal1_peaks_"+base) },
		"base pair probability bed",
	)
	if err != nil {
		cleanup()
		return nil, err
	}

	fileOps.MinFreeEnergyBed, err = createOutputFile(
		basePath+"_mfe.bed",
		func(f *os.File) error { return writeBedfileHeader(f, "signal3_peaks_"+base) },
		"minimum free energy bed",
	)
	if err != nil {
		cleanup()
		return nil, err
	}

	fileOps.ExtendedBasePairProbWig, err = createOutputFile(
		basePath+"_extbpprob.wig",
		func(f *os.File) error { return writeWigfileHeader(f, "signal4_"+base) },
		"extended base pair probability wig",
	)
	if err != nil {
		cleanup()
		return nil, err
	}

	return fileOps, nil
}
