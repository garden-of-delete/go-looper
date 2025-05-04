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
	// Create a temporary directory for test files
	tempDir, err := os.MkdirTemp("", "golooper_test_*")
	if err != nil {
		t.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tempDir)

	// Create a test config
	testConfig := &config.Config{
		OutfileName: filepath.Join(tempDir, "test_output"),
	}

	// Test file creation
	fileOps, err := CreateOutputFiles(testConfig)
	if err != nil {
		t.Fatalf("Failed to create output files: %v", err)
	}

	// Verify all files were created
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

	// Test wig file headers
	checkWigHeader(t, fileOps.BasePairProbWig, "Base Pair Probability")
	checkWigHeader(t, fileOps.AverageEnergyWig, "Average Energy")
	checkWigHeader(t, fileOps.MinFreeEnergyWig, "Minimum Free Energy")
	checkWigHeader(t, fileOps.ExtendedBasePairProbWig, "Extended Base Pair Probability")

	// Test bed file headers
	checkBedHeader(t, fileOps.BasePairProbBed, "Base Pair Probability")
	checkBedHeader(t, fileOps.MinFreeEnergyBed, "Minimum Free Energy")

	// Test file closing
	if err := fileOps.Close(); err != nil {
		t.Errorf("Failed to close files: %v", err)
	}

	// Verify files are closed by trying to read them
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

	// Reset file position to beginning
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

	// Reset file position to beginning
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
