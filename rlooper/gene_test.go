package rlooper

import (
	"os"
	"path/filepath"
	"sync"
	"testing"
)

func TestComputeStructuresSerial(t *testing.T) {
	// Get the absolute path to the project root
	wd, err := os.Getwd()
	if err != nil {
		t.Fatalf("Failed to get working directory: %v", err)
	}
	// Navigate up to the project root
	projectRoot := filepath.Dir(wd)

	gene := NewGene(filepath.Join(projectRoot, "res/gattaca.fa"))
	model := NewParamsReasonableDefaults()
	minLoopLength := 2

	result := gene.computeStructuresSerial(&model, minLoopLength, false)
	result2 := gene.computeStructuresSerial(&model, minLoopLength, true)

	if len(result) != 21 {
		t.Errorf("Expected 21 structures, got %d", len(result))
	}
	if len(result2) != 42 {
		t.Errorf("Expected 42 structures, got %d", len(result2))
	}
	// TODO: compare energetics for a structure against manual calculation
}

func TestComputeStructuresConcurrent(t *testing.T) {
	// Get the absolute path to the project root
	wd, err := os.Getwd()
	if err != nil {
		t.Fatalf("Failed to get working directory: %v", err)
	}
	// Navigate up to the project root
	projectRoot := filepath.Dir(wd)
	gene := NewGene(filepath.Join(projectRoot, "res/gattaca.fa"))
	model := NewParamsReasonableDefaults()
	minLoopLength := 2
	ec := &ExecutionContext{
		NumThreads: 2,
		WaitGroup:  &sync.WaitGroup{},
	}

	result := gene.computeStructuresConcurrent(ec, &model, minLoopLength, false)
	if len(result) != 21 {
		t.Errorf("Expected 21 structures, got %d", len(result))
	}
}
