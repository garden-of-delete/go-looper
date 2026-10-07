package cmd

import (
	"fmt"
	"os"

	"golooper/sim"

	"github.com/spf13/cobra"
)

var perloopCmd = &cobra.Command{
	Use:   "perloop",
	Short: "Run perloop analysis on nucleic acid sequences",
	Long: `Run perloop analysis on nucleic acid sequences to analyze R-loop formation
and energetics. This command takes an input file containing sequence data and
generates output at the specified path.`,
	PreRunE: func(cmd *cobra.Command, args []string) error {
		return requireInputOutput()
	},
	Run: func(cmd *cobra.Command, args []string) {
		if err := sim.SimulationA(&cfg); err != nil {
			fmt.Printf("Error running simulation: %v\n", err)
			os.Exit(1)
		}
	},
}

func init() {
	rootCmd.AddCommand(perloopCmd)
}
