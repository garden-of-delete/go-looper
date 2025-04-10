package main

import (
	"fmt"
	"go-looper/src"
)

func main() {

	//gene := src.NewGene("res/pfc53_full.fa")
	gene := src.NewGene("res/gattaca.fa")
	// rloopModel := src.NewParamsReasonableDefaults() // best set of experimentally validated params
	// structures := gene.ComputeStructures(&rloopModel)

	windows := src.FromLinearWindows(gene.Sequence, 2)
	fmt.Println(windows)

	// for _, w := range windows {
	// 	w.Start
	// }

	fmt.Println(gene.GeneName)
}
