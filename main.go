package main

import (
	"log"
	"os"

	"golooper/cmd"
)

func main() {
	if err := cmd.Execute(); err != nil {
		log.Fatal(err)
		os.Exit(1)
	}
}
