package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

var versionCmd = &cobra.Command{
	Use:   "version",
	Short: "Print the version number of Looper",
	Long:  `All software has versions. This is Looper's`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("golooper v0.5.0")
	},
}

// init registers the version subcommand on the root command.
func init() {
	rootCmd.AddCommand(versionCmd)
}
