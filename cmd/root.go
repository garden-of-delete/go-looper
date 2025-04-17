package cmd

import (
	"fmt"
	"strconv"

	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

// Config holds the global configuration values for the application
type Config struct {
	NucleationFreeEnergy *float64
	SuperhelicityDomain  *int
	AutoDomainSize       bool
	SuperhelicalDensity  *float64
	MinRLoopLength       *int
	Reverse              *bool
	Complement           *bool
	Unconstrained        *bool
	Invert               bool
	Dump                 bool
	Circular             bool
	Residuals            bool
	LocalAverageEnergy   bool
	Homopolymer          *float64
	InfileName           string
	OutfileName          string
}

var (
	// config holds the global configuration
	config Config
	// temporary variables to hold flag values
	nucleationFreeEnergy float64
	superhelicityDomain  string
	superhelicalDensity  float64
	minRLoopLength       int
	reverse              bool
	complement           bool
	unconstrained        bool
	homopolymer          float64
	infilename           string
	outfilename          string
)

var rootCmd = &cobra.Command{
	Use:   "golooper [subcommand] [input_file] [output_path] [options]",
	Short: "golooper is a CLI application for running biophysical simulations on nucleic acid energetics.",
	Long: `Go Looper is a CLI application that provides various commands
for running biophysical simulations on nucleic acid energetics with an emphasis on R-loops (genomic DNA/RNA hybrids).`,
	PersistentPreRun: func(cmd *cobra.Command, args []string) {
		config.InfileName = infilename
		config.OutfileName = outfilename

		// Set the pointers in the config only if the flags are set
		cmd.Flags().Visit(func(f *pflag.Flag) {
			switch f.Name {
			case "a":
				config.NucleationFreeEnergy = &nucleationFreeEnergy
			case "N":
				if superhelicityDomain == "auto" {
					config.AutoDomainSize = true
				} else {
					value, err := strconv.Atoi(superhelicityDomain)
					if err == nil {
						config.SuperhelicityDomain = &value
					}
				}
			case "sigma":
				config.SuperhelicalDensity = &superhelicalDensity
			case "minlength":
				config.MinRLoopLength = &minRLoopLength
			case "reverse":
				config.Reverse = &reverse
			case "complement":
				config.Complement = &complement
			case "unconstrained":
				config.Unconstrained = &unconstrained
			case "homopolymer":
				config.Homopolymer = &homopolymer
			}
		})
	},
}

// GetConfig returns the current configuration
func GetConfig() Config {
	return config
}

// Execute adds all child commands to the root command and sets flags appropriately.
func Execute() error {
	return rootCmd.Execute()
}

func init() {
	rootCmd.PersistentFlags().Float64VarP(&nucleationFreeEnergy, "a", "a", 0.0, "nucleation free energy in Kcal/mol")
	rootCmd.PersistentFlags().StringVarP(&superhelicityDomain, "N", "N", "0", "size of the superhelicity domain in nucleotides (use 'auto' for automatic sizing)")
	rootCmd.PersistentFlags().Float64VarP(&superhelicalDensity, "sigma", "s", 0.0, "superhelical density as a percentage (e.g., 0.07 for +7%)")
	rootCmd.PersistentFlags().IntVarP(&minRLoopLength, "minlength", "m", 0, "minimum length of an R-loop in nucleotides")
	rootCmd.PersistentFlags().BoolVarP(&reverse, "reverse", "r", false, "reverse the direction of the simulation")
	rootCmd.PersistentFlags().BoolVarP(&complement, "complement", "c", false, "use the complement strand for the simulation")
	rootCmd.PersistentFlags().BoolVarP(&unconstrained, "unconstrained", "u", false, "set the superhelicity modeling to unconstrained")
	rootCmd.PersistentFlags().BoolVarP(&config.Invert, "invert", "i", false, "invert the input sequence")
	rootCmd.PersistentFlags().BoolVarP(&config.Dump, "dump", "d", false, "dump all structures computed by the program to file")
	rootCmd.PersistentFlags().BoolVarP(&config.Circular, "circular", "C", false, "treat sequence as circular")
	rootCmd.PersistentFlags().BoolVarP(&config.Residuals, "residuals", "R", false, "calculate and output residual superhelicity for each structure")
	rootCmd.PersistentFlags().BoolVarP(&config.LocalAverageEnergy, "local-average-energy", "l", false, "use local average energy for the simulation")
	rootCmd.PersistentFlags().Float64VarP(&homopolymer, "homopolymer", "H", 0.0, "override base pairing energetics with constant value in Kcal/mol")
	rootCmd.PersistentFlags().StringVarP(&infilename, "input", "f", "", "input file name (required)")
	rootCmd.PersistentFlags().StringVarP(&outfilename, "output", "o", "", "output file name (required)")

	rootCmd.AddCommand(&cobra.Command{
		Use:   "show-config",
		Short: "Display the current configuration values",
		Long: `Display the current configuration values that will be used for the simulation.
Values not explicitly set via command line flags will use the model's default values.`,
		Run: func(cmd *cobra.Command, args []string) {
			fmt.Println("Current Configuration:")
			fmt.Println("---------------------")
			if config.NucleationFreeEnergy == nil {
				fmt.Println("Nucleation Free Energy (--a): not set (will use model default)")
			} else {
				fmt.Printf("Nucleation Free Energy (--a): %.2f Kcal/mol\n", *config.NucleationFreeEnergy)
			}
			if config.AutoDomainSize {
				fmt.Println("Superhelicity Domain (--N): auto (will be determined automatically)")
			} else if config.SuperhelicityDomain == nil {
				fmt.Println("Superhelicity Domain (--N): not set (will use model default)")
			} else {
				fmt.Printf("Superhelicity Domain (--N): %d nucleotides\n", *config.SuperhelicityDomain)
			}
			if config.SuperhelicalDensity == nil {
				fmt.Println("Superhelical Density (--sigma): not set (will use model default)")
			} else {
				fmt.Printf("Superhelical Density (--sigma): %.3f (%.1f%%)\n", *config.SuperhelicalDensity, *config.SuperhelicalDensity*100)
			}
			if config.MinRLoopLength == nil {
				fmt.Println("Minimum R-loop Length (--minlength): not set (will use model default)")
			} else {
				fmt.Printf("Minimum R-loop Length (--minlength): %d nucleotides\n", *config.MinRLoopLength)
			}
			if config.Reverse == nil {
				fmt.Println("Reverse Direction (--reverse): not set (will use model default)")
			} else {
				fmt.Printf("Reverse Direction (--reverse): %v\n", *config.Reverse)
			}
			if config.Complement == nil {
				fmt.Println("Complement Strand (--complement): not set (will use model default)")
			} else {
				fmt.Printf("Complement Strand (--complement): %v\n", *config.Complement)
			}
			if config.Unconstrained == nil {
				fmt.Println("Unconstrained Model (--unconstrained): not set (will use model default)")
			} else {
				fmt.Printf("Unconstrained Model (--unconstrained): %v\n", *config.Unconstrained)
			}
			if config.Homopolymer == nil {
				fmt.Println("Homopolymer Energy (--homopolymer): not set (will use model default)")
			} else {
				fmt.Printf("Homopolymer Energy (--homopolymer): %.2f Kcal/mol\n", *config.Homopolymer)
			}
			fmt.Printf("Invert Output (--invert): %v\n", config.Invert)
			fmt.Printf("Dump Calculations (--dump): %v\n", config.Dump)
			fmt.Printf("Circular Sequence (--circular): %v\n", config.Circular)
			fmt.Printf("Calculate Residuals (--residuals): %v\n", config.Residuals)
			fmt.Printf("Local Average Energy (--local-average-energy): %v\n", config.LocalAverageEnergy)
			fmt.Println("---------------------")
		},
	})

	rootCmd.MarkPersistentFlagRequired("input")
	rootCmd.MarkPersistentFlagRequired("output")
}
