package cmd

import (
	"fmt"
	"os"

	"github.com/spf13/cobra"
)

var perloopCmd = &cobra.Command{
	Use:   "perloop [input_file] [output_path]",
	Short: "Run perloop analysis on nucleic acid sequences",
	Long: `Run perloop analysis on nucleic acid sequences to analyze R-loop formation
and energetics. This command takes an input file containing sequence data and
generates output at the specified path.`,
	Args: cobra.ExactArgs(2),
	Run: func(cmd *cobra.Command, args []string) {
		inputFile := args[0]
		outputPath := args[1]

		// Basic validation
		if _, err := os.Stat(inputFile); os.IsNotExist(err) {
			fmt.Printf("Error: input file %s does not exist\n", inputFile)
			os.Exit(1)
		}

		// TODO: Implement the actual perloop analysis logic here
		fmt.Printf("Running perloop analysis on %s, output will be saved to %s\n", inputFile, outputPath)
	},
}

func init() {
	rootCmd.AddCommand(perloopCmd)

	// Add flags specific to perloop command
	perloopCmd.Flags().BoolP("verbose", "v", false, "Enable verbose output")
	perloopCmd.Flags().StringP("config", "c", "", "Path to configuration file")
}
