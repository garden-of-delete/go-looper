package cmd

import (
	"fmt"
	"strconv"

	"golooper/config"

	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

// cfg holds the global configuration
var cfg config.Config

// temporary variables to hold flag values
var nucleationFreeEnergy float64
var superhelicityDomain string
var superhelicalDensity float64
var minRLoopLength int
var reverse bool
var complement bool
var unconstrained bool
var homopolymer float64
var infilename string
var outfilename string

var rootCmd = &cobra.Command{
	Use:   "golooper [subcommand] [input_file] [output_path] [options]",
	Short: "golooper is a CLI application for running biophysical simulations on nucleic acid energetics.",
	Long: `Go Looper is a CLI application that provides various commands
for running biophysical simulations on nucleic acid energetics with an emphasis on R-loops (genomic DNA/RNA hybrids).`,
	PersistentPreRun: func(cmd *cobra.Command, args []string) {
		cfg.InfileName = infilename
		cfg.OutfileName = outfilename

		// Set the pointers in the config only if the flags are set
		cmd.Flags().Visit(func(f *pflag.Flag) {
			switch f.Name {
			case "a":
				cfg.NucleationFreeEnergy = &nucleationFreeEnergy
			case "N":
				if superhelicityDomain == "auto" {
					cfg.AutoDomainSize = true
				} else {
					value, err := strconv.Atoi(superhelicityDomain)
					if err == nil {
						cfg.SuperhelicityDomain = &value
					}
				}
			case "sigma":
				cfg.SuperhelicalDensity = &superhelicalDensity
			case "minlength":
				cfg.MinRLoopLength = &minRLoopLength
			case "reverse":
				cfg.Reverse = &reverse
			case "complement":
				cfg.Complement = &complement
			case "unconstrained":
				cfg.Unconstrained = &unconstrained
			case "homopolymer":
				cfg.Homopolymer = &homopolymer
			}
		})
	},
}

// GetConfig returns the current configuration
func GetConfig() config.Config {
	return cfg
}

// Execute adds all child commands to the root command and sets flags appropriately.
func Execute() error {
	return rootCmd.Execute()
}

func init() {
	initFlags()

	rootCmd.AddCommand(&cobra.Command{
		Use:   "show-config",
		Short: "Display the current configuration values",
		Long: `Display the current configuration values that will be used for the simulation.
Values not explicitly set via command line flags will use the model's default values.`,
		Run: func(cmd *cobra.Command, args []string) {
			fmt.Println("Current Configuration:")
			fmt.Println("---------------------")
			if cfg.NucleationFreeEnergy == nil {
				fmt.Println("Nucleation Free Energy (--a): not set (will use model default)")
			} else {
				fmt.Printf("Nucleation Free Energy (--a): %.2f Kcal/mol\n", *cfg.NucleationFreeEnergy)
			}
			if cfg.AutoDomainSize {
				fmt.Println("Superhelicity Domain (--N): auto (will be determined automatically)")
			} else if cfg.SuperhelicityDomain == nil {
				fmt.Println("Superhelicity Domain (--N): not set (will use model default)")
			} else {
				fmt.Printf("Superhelicity Domain (--N): %d nucleotides\n", *cfg.SuperhelicityDomain)
			}
			if cfg.SuperhelicalDensity == nil {
				fmt.Println("Superhelical Density (--sigma): not set (will use model default)")
			} else {
				fmt.Printf("Superhelical Density (--sigma): %.3f (%.1f%%)\n", *cfg.SuperhelicalDensity, *cfg.SuperhelicalDensity*100)
			}
			if cfg.MinRLoopLength == nil {
				fmt.Println("Minimum R-loop Length (--minlength): not set (will use model default)")
			} else {
				fmt.Printf("Minimum R-loop Length (--minlength): %d nucleotides\n", *cfg.MinRLoopLength)
			}
			if cfg.Reverse == nil {
				fmt.Println("Reverse Direction (--reverse): not set (will use model default)")
			} else {
				fmt.Printf("Reverse Direction (--reverse): %v\n", *cfg.Reverse)
			}
			if cfg.Complement == nil {
				fmt.Println("Complement Strand (--complement): not set (will use model default)")
			} else {
				fmt.Printf("Complement Strand (--complement): %v\n", *cfg.Complement)
			}
			if cfg.Unconstrained == nil {
				fmt.Println("Unconstrained Model (--unconstrained): not set (will use model default)")
			} else {
				fmt.Printf("Unconstrained Model (--unconstrained): %v\n", *cfg.Unconstrained)
			}
			if cfg.Homopolymer == nil {
				fmt.Println("Homopolymer Energy (--homopolymer): not set (will use model default)")
			} else {
				fmt.Printf("Homopolymer Energy (--homopolymer): %.2f Kcal/mol\n", *cfg.Homopolymer)
			}
			fmt.Printf("Invert Output (--invert): %v\n", cfg.Invert)
			fmt.Printf("Dump Calculations (--dump): %v\n", cfg.Dump)
			fmt.Printf("Circular Sequence (--circular): %v\n", cfg.Circular)
			fmt.Printf("Calculate Residuals (--residuals): %v\n", cfg.Residuals)
			fmt.Printf("Local Average Energy (--local-average-energy): %v\n", cfg.LocalAverageEnergy)
			fmt.Println("---------------------")
		},
	})

	rootCmd.MarkPersistentFlagRequired("input")
	rootCmd.MarkPersistentFlagRequired("output")
}

func initFlags() {
	rootCmd.PersistentFlags().Float64VarP(&nucleationFreeEnergy, "a", "a", 0.0, "nucleation free energy in Kcal/mol")
	rootCmd.PersistentFlags().StringVarP(&superhelicityDomain, "N", "N", "0", "size of the superhelicity domain in nucleotides (use 'auto' for automatic sizing)")
	rootCmd.PersistentFlags().Float64VarP(&superhelicalDensity, "sigma", "s", 0.0, "superhelical density as a percentage (e.g., 0.07 for +7%)")
	rootCmd.PersistentFlags().IntVarP(&minRLoopLength, "minlength", "m", 0, "minimum length of an R-loop in nucleotides")
	rootCmd.PersistentFlags().BoolVarP(&reverse, "reverse", "r", false, "reverse the direction of the simulation")
	rootCmd.PersistentFlags().BoolVarP(&complement, "complement", "c", false, "use the complement strand for the simulation")
	rootCmd.PersistentFlags().BoolVarP(&unconstrained, "unconstrained", "u", false, "set the superhelicity modeling to unconstrained")
	rootCmd.PersistentFlags().BoolVarP(&cfg.Invert, "invert", "i", false, "invert the input sequence")
	rootCmd.PersistentFlags().BoolVarP(&cfg.Dump, "dump", "d", false, "dump all structures computed by the program to file")
	rootCmd.PersistentFlags().BoolVarP(&cfg.Circular, "circular", "C", false, "treat sequence as circular")
	rootCmd.PersistentFlags().BoolVarP(&cfg.Residuals, "residuals", "R", false, "calculate and output residual superhelicity for each structure")
	rootCmd.PersistentFlags().BoolVarP(&cfg.LocalAverageEnergy, "local-average-energy", "l", false, "use local average energy for the simulation")
	rootCmd.PersistentFlags().Float64VarP(&homopolymer, "homopolymer", "H", 0.0, "override base pairing energetics with constant value in Kcal/mol")
	rootCmd.PersistentFlags().StringVarP(&infilename, "input", "f", "", "input file name (required)")
	rootCmd.PersistentFlags().StringVarP(&outfilename, "output", "o", "", "output file name (required)")
}
