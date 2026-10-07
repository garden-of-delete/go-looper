package cmd

import (
	"fmt"
	"os"

	"golooper/sim"

	"github.com/spf13/cobra"
)

var ensembleCmd = &cobra.Command{
	Use:   "ensemble",
	Short: "Run the full ensemble simulation on nucleic acid sequences",
	Long: `Run the full ensemble simulation on nucleic acid sequences to compute
R-loop formation energetics and signals, writing wig/bed tracks and optional
dump, top-structure and residual reports.`,
	PreRunE: func(cmd *cobra.Command, args []string) error {
		return requireInputOutput()
	},
	Run: func(cmd *cobra.Command, args []string) {
		if err := sim.RunEnsembleSim(&cfg, threads, os.Stdout); err != nil {
			fmt.Printf("Error running simulation: %v\n", err)
			os.Exit(1)
		}
	},
}

func init() {
	rootCmd.AddCommand(ensembleCmd)
}
