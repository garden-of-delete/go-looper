package cmd

import (
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
	SilenceUsage: true,
	PreRunE: func(cmd *cobra.Command, args []string) error {
		return requireInputOutput()
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		return sim.RunEnsembleSim(&cfg, threads, os.Stdout)
	},
}

func init() {
	rootCmd.AddCommand(ensembleCmd)
}
