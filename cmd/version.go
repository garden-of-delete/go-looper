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
		fmt.Println("Looper v0.1.0")
	},
}

func init() {
	rootCmd.AddCommand(versionCmd)
}
