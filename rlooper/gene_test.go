package rlooper

import (
	"fmt"
	"os"
	"path/filepath"
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
	fmt.Println(len(result2))
	// TODO: compare energetics for a structure against manual calculation


	// if len(result2) != 21 {
	// 	t.Errorf("Expected 21 structures, got %d", len(result2))
	// }
}
