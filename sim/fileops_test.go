package sim

import (
	"bufio"
	"fmt"
	"golooper/config"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestFileOps(t *testing.T) {
	tempDir, err := os.MkdirTemp("", "golooper_test_*")
	if err != nil {
		t.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tempDir)

	testConfig := &config.Config{
		OutfileName: filepath.Join(tempDir, "test_output"),
	}
	base := filepath.Base(testConfig.OutfileName)

	fileOps, err := CreateOutputFiles(testConfig)
	if err != nil {
		t.Fatalf("Failed to create output files: %v", err)
	}

	expectedFiles := []string{
		"_bpprob.wig",
		"_avgG.wig",
		"_mfe.wig",
		"_bpprob.bed",
		"_mfe.bed",
		"_extbpprob.wig",
	}

	for _, suffix := range expectedFiles {
		filePath := testConfig.OutfileName + suffix
		if _, err := os.Stat(filePath); os.IsNotExist(err) {
			t.Errorf("Expected file %s was not created", filePath)
		}
	}

	checkWigHeader(t, fileOps.BasePairProbWig, "signal1_"+base)
	checkWigHeader(t, fileOps.AverageEnergyWig, "signal2_"+base)
	checkWigHeader(t, fileOps.MinFreeEnergyWig, "signal3_"+base)
	checkWigHeader(t, fileOps.ExtendedBasePairProbWig, "signal4_"+base)

	checkBedHeader(t, fileOps.BasePairProbBed, "signal1_peaks_"+base)
	checkBedHeader(t, fileOps.MinFreeEnergyBed, "signal3_peaks_"+base)

	if err := fileOps.Close(); err != nil {
		t.Errorf("Failed to close files: %v", err)
	}

	for _, suffix := range expectedFiles {
		filePath := testConfig.OutfileName + suffix
		file, err := os.OpenFile(filePath, os.O_RDONLY, 0)
		if err != nil {
			t.Errorf("Failed to open file %s for verification: %v", filePath, err)
			continue
		}
		file.Close()
	}
}

func checkWigHeader(t *testing.T, file *os.File, expectedName string) {
	t.Helper()

	if _, err := file.Seek(0, 0); err != nil {
		t.Errorf("Failed to seek to beginning of file: %v", err)
		return
	}

	reader := bufio.NewReader(file)
	header, err := reader.ReadString('\n')
	if err != nil {
		t.Errorf("Failed to read wig header: %v", err)
		return
	}

	expectedHeader := fmt.Sprintf("track type=wiggle_0 name=\"%s\" visibility=full autoscale=off color=50,150,255 priority=10\n", expectedName)
	if !strings.Contains(header, expectedHeader) {
		t.Errorf("Wig header mismatch. Got: %s, Expected: %s", header, expectedHeader)
	}
}

func checkBedHeader(t *testing.T, file *os.File, expectedName string) {
	t.Helper()

	if _, err := file.Seek(0, 0); err != nil {
		t.Errorf("Failed to seek to beginning of file: %v", err)
		return
	}

	reader := bufio.NewReader(file)
	header, err := reader.ReadString('\n')
	if err != nil {
		t.Errorf("Failed to read bed header: %v", err)
		return
	}

	expectedHeader := fmt.Sprintf("track name=rLooper description=\"%s\" useScore=1\n", expectedName)
	if !strings.Contains(header, expectedHeader) {
		t.Errorf("Bed header mismatch. Got: %s, Expected: %s", header, expectedHeader)
	}
}
