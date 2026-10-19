package cmd

import (
	"fmt"

	"golooper/sim"

	"github.com/spf13/cobra"
)

var sigmaMin float64
var sigmaMax float64
var sigmaStep float64

var perloopCmd = &cobra.Command{
	Use:   "perloop",
	Short: "Run perloop analysis on nucleic acid sequences",
	Long: `Run perloop analysis on nucleic acid sequences to analyze R-loop formation
and energetics. This command takes an input file containing sequence data and
generates output at the specified path.`,
	SilenceUsage: true,
	PreRunE: func(cmd *cobra.Command, args []string) error {
		if err := requireInputOutput(); err != nil {
			return err
		}
		return validateSigmaRange(sigmaMin, sigmaMax, sigmaStep)
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		return sim.RunPerloop(&cfg, sigmaMin, sigmaMax, sigmaStep, threads)
	},
}

func validateSigmaRange(min, max, step float64) error {
	if step <= 0 {
		return fmt.Errorf("--sigma-step must be > 0, got %g", step)
	}
	if min > max {
		return fmt.Errorf("--sigma-min (%g) must be <= --sigma-max (%g)", min, max)
	}
	return nil
}

func init() {
	perloopCmd.Flags().Float64Var(&sigmaMin, "sigma-min", 0, "lower bound of the superhelicity sweep")
	perloopCmd.Flags().Float64Var(&sigmaMax, "sigma-max", 0, "upper bound of the superhelicity sweep")
	perloopCmd.Flags().Float64Var(&sigmaStep, "sigma-step", 0.01, "step size of the superhelicity sweep")
	perloopCmd.MarkFlagRequired("sigma-min")
	perloopCmd.MarkFlagRequired("sigma-max")
	rootCmd.AddCommand(perloopCmd)
}
