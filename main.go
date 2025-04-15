package main

import (
	"log"
	"os"

	"github.com/garden-of-delete/go-looper/cmd"
)

func main() {
	if err := cmd.Execute(); err != nil {
		log.Fatal(err)
		os.Exit(1)
	}
}
