package sim

import (
	"fmt"
	"io"
	"math"
	"os"
	"slices"

	"golooper/config"
	"golooper/rlooper"
)

const defaultMinRLoopLength = 2
const defaultSensitivity = 1

func ModelFromConfig(cfg *config.Config) *rlooper.EnergyModel {
	m := rlooper.NewDefaultModel()
	if cfg.NucleationFreeEnergy != nil {
		m.Nucleation = *cfg.NucleationFreeEnergy
	}
	if cfg.SuperhelicalDensity != nil {
		m.SetSigma(*cfg.SuperhelicalDensity)
	}
	if cfg.SuperhelicityDomain != nil {
		m.SetN(float64(*cfg.SuperhelicityDomain))
	}
	if cfg.Homopolymer != nil {
		m.Homopolymer = true
		m.OverrideEnergy = *cfg.Homopolymer
	}
	if cfg.Unconstrained != nil {
		m.Unconstrained = *cfg.Unconstrained
	}
	return m
}

func minLengthFromConfig(cfg *config.Config) int {
	if cfg.MinRLoopLength != nil {
		return *cfg.MinRLoopLength
	}
	return defaultMinRLoopLength
}

func sensitivityFromConfig(cfg *config.Config) int {
	if cfg.Sensitivity != nil {
		return *cfg.Sensitivity
	}
	return defaultSensitivity
}

func topFromConfig(cfg *config.Config) int {
	if cfg.Top != nil {
		return *cfg.Top
	}
	return 0
}

func boolPtr(p *bool) bool {
	return p != nil && *p
}

type ensembleParams struct {
	model      *rlooper.EnergyModel
	minLength  int
	complement bool
	reverse    bool
	top        int
	threshold  float64
	threads    int
}

func paramsFromConfig(cfg *config.Config, threads int) ensembleParams {
	return ensembleParams{
		model:      ModelFromConfig(cfg),
		minLength:  minLengthFromConfig(cfg),
		complement: boolPtr(cfg.Complement),
		reverse:    boolPtr(cfg.Reverse),
		top:        topFromConfig(cfg),
		threshold:  math.Pow(10, float64(sensitivityFromConfig(cfg))),
		threads:    threads,
	}
}

func RunEnsembleSim(cfg *config.Config, threads int, stdout io.Writer) error {
	infile, err := os.Open(cfg.InfileName)
	if err != nil {
		return fmt.Errorf("error opening input file: %v", err)
	}
	defer infile.Close()

	genes, err := rlooper.ReadGenes(infile)
	if err != nil {
		return fmt.Errorf("error reading genes: %v", err)
	}

	outFiles, err := CreateOutputFiles(cfg)
	if err != nil {
		return fmt.Errorf("error creating output files: %v", err)
	}
	defer outFiles.Close()

	p := paramsFromConfig(cfg, threads)

	for _, g := range genes {
		if err := processGene(stdout, outFiles, g, cfg, p); err != nil {
			return err
		}
	}

	return nil
}

func processGene(stdout io.Writer, outFiles *FileOps, g *rlooper.Gene, cfg *config.Config, p ensembleParams) error {
	fmt.Fprintf(stdout, "processing gene: %s...", g.Name)

	if cfg.AutoDomainSize {
		p.model.SetN(float64(len(g.Sequence)))
	}
	g.ApplyOrientation(p.complement, p.reverse, cfg.Invert)

	ens := rlooper.ComputeEnsemble(g.Sequence, p.model, p.minLength, p.threads)
	fmt.Fprintln(stdout, ens.NumStructures)
	fmt.Fprintf(stdout, "P(ground state)= %s\n", formatValue(ens.PGround))

	sig := rlooper.ComputeSignals(g.Sequence, p.model, p.minLength, ens, cfg.LocalAverageEnergy, p.threads)
	if g.Loc.Strand == "-" {
		slices.Reverse(sig.BPProb)
		slices.Reverse(sig.MFE)
		if cfg.LocalAverageEnergy {
			slices.Reverse(sig.AvgG)
		}
	}

	if err := WriteWigSignal(outFiles.BasePairProbWig, g, sig.BPProb); err != nil {
		return err
	}
	if cfg.LocalAverageEnergy {
		if err := WriteWigSignal(outFiles.AverageEnergyWig, g, sig.AvgG); err != nil {
			return err
		}
	}
	if err := WriteWigSignal(outFiles.MinFreeEnergyWig, g, sig.MFE); err != nil {
		return err
	}

	if cfg.Bedfile {
		bpPeaks := CallPeaksAbsolute(sig.BPProb, p.threshold, g)
		if err := WriteBedPeaks(outFiles.BasePairProbBed, g, bpPeaks); err != nil {
			return err
		}
		mfePeaks := CallPeaksAbsolute(sig.MFE, p.threshold, g)
		if err := WriteBedPeaks(outFiles.MinFreeEnergyBed, g, mfePeaks); err != nil {
			return err
		}
	}

	fmt.Fprintln(stdout, "complete!")

	if cfg.Residuals {
		res := rlooper.ComputeResiduals(g.Sequence, p.model, p.minLength, ens)
		if err := PrintResiduals(stdout, res); err != nil {
			return err
		}
	}

	var recs []rlooper.StructureRecord
	if p.top > 0 || cfg.Dump {
		recs = rlooper.AllStructures(g.Sequence, p.model, p.minLength, ens)
	}

	if p.top > 0 {
		if err := PrintTop(stdout, g, p.model, recs, p.top, cfg.Residuals); err != nil {
			return err
		}
	}

	if cfg.Dump {
		if err := dumpGene(stdout, cfg, g, p.model, recs); err != nil {
			return err
		}
	}

	return nil
}

func dumpGene(stdout io.Writer, cfg *config.Config, g *rlooper.Gene, m *rlooper.EnergyModel, recs []rlooper.StructureRecord) error {
	dumpFile, err := createFileWithDir(cfg.OutfileName + g.Name + "_dump.txt")
	if err != nil {
		return fmt.Errorf("error creating dump file: %v", err)
	}
	if err := WriteDump(dumpFile, g, m.GroundStateEnergy(), recs); err != nil {
		dumpFile.Close()
		return err
	}
	if err := dumpFile.Close(); err != nil {
		return fmt.Errorf("error closing dump file: %v", err)
	}
	fmt.Fprintln(stdout, len(g.Sequence))
	return nil
}
