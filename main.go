package main

import (
	"os"

	"golooper/cmd"
)

// main runs the golooper CLI and exits with status 1 on any command error.
func main() {
	if err := cmd.Execute(); err != nil {
		os.Exit(1)
	}
}
