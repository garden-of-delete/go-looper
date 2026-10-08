package sim

import (
	"fmt"
	"os"

	"golooper/config"
	"golooper/rlooper"
)

func sigmaGrid(lo, hi, step float64) []float64 {
	var sigmas []float64
	for i := 0; ; i++ {
		s := lo + float64(i)*step
		if s > hi+1e-4 {
			break
		}
		sigmas = append(sigmas, s)
	}
	return sigmas
}

func RunPerloop(cfg *config.Config, lo, hi, step float64, threads int) error {
	infile, err := os.Open(cfg.InfileName)
	if err != nil {
		return fmt.Errorf("error opening input file: %v", err)
	}
	defer infile.Close()

	genes, err := rlooper.ReadGenes(infile)
	if err != nil {
		return fmt.Errorf("error reading genes: %v", err)
	}
	if len(genes) == 0 {
		return fmt.Errorf("no genes found in input file")
	}
	g := genes[0]

	p := paramsFromConfig(cfg, threads)
	g.ApplyOrientation(p.complement, p.reverse, cfg.Invert)

	outFile, err := createFileWithDir(cfg.OutfileName + "_perloop.txt")
	if err != nil {
		return fmt.Errorf("error creating perloop output file: %v", err)
	}
	defer outFile.Close()

	for _, sigma := range sigmaGrid(lo, hi, step) {
		p.model.SetSigma(sigma)
		ens := rlooper.ComputeEnsemble(g.Sequence, p.model, p.minLength, p.threads)
		if _, err := fmt.Fprintf(outFile, "%s %s\n", formatValue(sigma), formatValue(1-ens.PGround)); err != nil {
			return fmt.Errorf("error writing perloop output: %v", err)
		}
	}

	return nil
}
