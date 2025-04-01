package main

import (
	"fmt"
	"go-looper/genome"
)

func main() {

	gene := genome.NewGene("res/pfc53_full.fa")
	windower := genome.NewWindower(gene.Sequence, genome.WithMinWindowSize(2))
	windower

	fmt.Println(gene.GeneName)
}
