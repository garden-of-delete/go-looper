package main

import (
	"fmt"
	rlooper "go-looper/rlooper"
)

func main() {

	//gene := src.NewGene("res/pfc53_full.fa")
	gene := rlooper.NewGene("res/gattaca.fa")
	// rloopModel := src.NewParamsReasonableDefaults() // best set of experimentally validated params
	// structures := gene.ComputeStructures(&rloopModel)

	windows := rlooper.FromLinearWindows(gene.Sequence, 2)
	fmt.Println(windows)

	// for _, w := range windows {
	// 	w.Start
	// }

	fmt.Println(gene.GeneName)
}
