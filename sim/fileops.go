package sim

import (
	"bufio"
	"fmt"
	"golooper/config"
	"io"
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

type bufferedFile struct {
	*bufio.Writer
	f *os.File
}

func newBufferedFile(f *os.File) *bufferedFile {
	return &bufferedFile{Writer: bufio.NewWriter(f), f: f}
}

func (b *bufferedFile) Close() error {
	if err := b.Flush(); err != nil {
		b.f.Close()
		return fmt.Errorf("error flushing %s: %v", b.f.Name(), err)
	}
	if err := b.f.Close(); err != nil {
		return fmt.Errorf("error closing %s: %v", b.f.Name(), err)
	}
	return nil
}

func createOutputFile(path string, headerFunc func(io.Writer) error, headerName string) (*bufferedFile, error) {
	file, err := createFileWithDir(path)
	if err != nil {
		return nil, fmt.Errorf("error creating %s file at %s: %v", headerName, path, err)
	}
	buffered := newBufferedFile(file)
	if err := headerFunc(buffered); err != nil {
		buffered.Close()
		return nil, fmt.Errorf("error writing %s header: %v", headerName, err)
	}
	return buffered, nil
}

type FileOps struct {
	BasePairProbWig         *bufferedFile
	AverageEnergyWig        *bufferedFile
	MinFreeEnergyWig        *bufferedFile
	BasePairProbBed         *bufferedFile
	MinFreeEnergyBed        *bufferedFile
	ExtendedBasePairProbWig *bufferedFile
}

func writeWigfileHeader(w io.Writer, trackname string) error {
	header := fmt.Sprintf("track type=wiggle_0 name=\"%s\" visibility=full autoscale=off color=50,150,255 priority=10\n", trackname)
	if _, err := io.WriteString(w, header); err != nil {
		return fmt.Errorf("error writing wig header: %v", err)
	}
	return nil
}

func writeBedfileHeader(w io.Writer, trackname string) error {
	header := fmt.Sprintf("track name=rLooper description=\"%s\" useScore=1\n", trackname)
	if _, err := io.WriteString(w, header); err != nil {
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
		func(w io.Writer) error { return writeWigfileHeader(w, "signal1_"+base) },
		"base pair probability wig",
	)
	if err != nil {
		cleanup()
		return nil, err
	}

	fileOps.AverageEnergyWig, err = createOutputFile(
		basePath+"_avgG.wig",
		func(w io.Writer) error { return writeWigfileHeader(w, "signal2_"+base) },
		"average energy wig",
	)
	if err != nil {
		cleanup()
		return nil, err
	}

	fileOps.MinFreeEnergyWig, err = createOutputFile(
		basePath+"_mfe.wig",
		func(w io.Writer) error { return writeWigfileHeader(w, "signal3_"+base) },
		"minimum free energy wig",
	)
	if err != nil {
		cleanup()
		return nil, err
	}

	fileOps.BasePairProbBed, err = createOutputFile(
		basePath+"_bpprob.bed",
		func(w io.Writer) error { return writeBedfileHeader(w, "signal1_peaks_"+base) },
		"base pair probability bed",
	)
	if err != nil {
		cleanup()
		return nil, err
	}

	fileOps.MinFreeEnergyBed, err = createOutputFile(
		basePath+"_mfe.bed",
		func(w io.Writer) error { return writeBedfileHeader(w, "signal3_peaks_"+base) },
		"minimum free energy bed",
	)
	if err != nil {
		cleanup()
		return nil, err
	}

	fileOps.ExtendedBasePairProbWig, err = createOutputFile(
		basePath+"_extbpprob.wig",
		func(w io.Writer) error { return writeWigfileHeader(w, "signal4_"+base) },
		"extended base pair probability wig",
	)
	if err != nil {
		cleanup()
		return nil, err
	}

	return fileOps, nil
}
