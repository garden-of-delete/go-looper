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
	Run: func(cmd *cobra.Command, args []string) {
		// Validate required flags
		if cfg.InfileName == "" {
			fmt.Println("Error: input file is required")
			os.Exit(1)
		}
		if cfg.OutfileName == "" {
			fmt.Println("Error: output file is required")
			os.Exit(1)
		}

		// Run the simulation
		if err := sim.SimulationA(&cfg); err != nil {
			fmt.Printf("Error running simulation: %v\n", err)
			os.Exit(1)
		}
	},
}

func init() {
	rootCmd.AddCommand(perloopCmd)

	// Add flags specific to perloop command
	perloopCmd.Flags().BoolP("verbose", "v", false, "Enable verbose output")
}
