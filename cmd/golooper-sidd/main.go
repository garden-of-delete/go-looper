package main

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"golooper/rlooper"
	"golooper/sidd"
	"golooper/sim"

	"github.com/spf13/cobra"
)

type options struct {
	infile, outfile, tracks string
	domain, maxRuns         string
	clamp, circular         bool
	params                  sidd.Params
}

// parseDomain converts --N to the SIDD domain size, mapping "auto" to 0 and rejecting anything but a positive integer.
func parseDomain(s string) (float64, error) {
	if s == "auto" {
		return 0, nil
	}
	n, err := strconv.Atoi(s)
	if err != nil || n <= 0 {
		return 0, fmt.Errorf("--N must be a positive integer or \"auto\", got %q", s)
	}
	return float64(n), nil
}

// parseMaxRuns parses --max-runs, accepting only integers 1 through 4 (SIST's structural cap).
func parseMaxRuns(s string) (int, error) {
	n, err := strconv.Atoi(s)
	if err != nil || n < 1 || n > 4 {
		return 0, fmt.Errorf("--max-runs must be an integer in the range 1-4, got %q", s)
	}
	return n, nil
}

// resolveParams folds the flags into a sidd.Params, rejecting --clamp with --circular and out-of-range --N or --max-runs.
func (o *options) resolveParams() (sidd.Params, error) {
	p := o.params
	if o.clamp && o.circular {
		return p, errors.New("--clamp and --circular are mutually exclusive")
	}
	switch {
	case o.clamp:
		p.Geometry = sidd.GeomClamp
	case o.circular:
		p.Geometry = sidd.GeomCircular
	}
	n, err := parseDomain(o.domain)
	if err != nil {
		return p, err
	}
	p.N = n
	p.MaxRuns, err = parseMaxRuns(o.maxRuns)
	if err != nil {
		return p, err
	}
	return p, nil
}

// run reads the single input record, runs the SIDD model, writes the report to -o or stdout, and optionally writes the bpprob and gx wig tracks.
func (o *options) run(cmd *cobra.Command) error {
	p, err := o.resolveParams()
	if err != nil {
		return err
	}

	in, err := os.Open(o.infile)
	if err != nil {
		return err
	}
	defer in.Close()
	genes, err := readSequence(in, o.infile)
	if err != nil {
		return err
	}
	if len(genes) != 1 {
		return fmt.Errorf("input must contain exactly one FASTA record, found %d", len(genes))
	}
	gene := genes[0]

	model, err := sidd.NewModel(gene.Sequence, p)
	if err != nil {
		return err
	}
	result, err := sidd.Run(model)
	if err != nil {
		return err
	}

	if o.outfile == "" {
		err = sidd.WriteReport(cmd.OutOrStdout(), o.infile, model, result)
	} else {
		err = writeFile(o.outfile, func(f *os.File) error {
			return sidd.WriteReport(f, o.infile, model, result)
		})
	}
	if err != nil {
		return err
	}

	if o.tracks == "" {
		return nil
	}
	setTrackLocation(gene, o.infile)
	for suffix, signal := range map[string][]float64{"_sidd_bpprob.wig": result.P, "_sidd_gx.wig": result.G} {
		err := writeFile(o.tracks+suffix, func(f *os.File) error {
			return sim.WriteWigSignal(f, gene, signal)
		})
		if err != nil {
			return err
		}
	}
	return nil
}

// readSequence parses FASTA from r, accepting SIST-style headerless raw sequence by synthesizing a header from name.
func readSequence(r io.Reader, name string) ([]*rlooper.Gene, error) {
	data, err := io.ReadAll(r)
	if err != nil {
		return nil, err
	}
	if len(data) > 0 && !startsWithFasta(data) {
		data = append([]byte(">"+baseName(name)+"\n"), data...)
	}
	return rlooper.ReadGenes(bytes.NewReader(data))
}

// baseName returns the file's base name with whitespace runs replaced by underscores, so it is a valid FASTA/wig token.
func baseName(path string) string {
	return strings.Join(strings.Fields(filepath.Base(path)), "_")
}

// setTrackLocation fills the gene's chrom/start/end when the input had no UCSC range= header, so wig output is valid.
func setTrackLocation(g *rlooper.Gene, infile string) {
	if g.Loc.Chromosome != "" {
		return
	}
	g.Loc.Chromosome = g.Name
	if g.Loc.Chromosome == "" {
		g.Loc.Chromosome = baseName(infile)
	}
	g.Loc.Start = 1
	g.Loc.End = int64(len(g.Sequence))
}

// startsWithFasta reports whether the first non-blank line begins a FASTA record or comment ('>' or ';').
func startsWithFasta(data []byte) bool {
	for len(data) > 0 {
		line, rest, _ := bytes.Cut(data, []byte("\n"))
		if len(bytes.TrimSpace(line)) > 0 {
			return line[0] == '>' || line[0] == ';'
		}
		data = rest
	}
	return false
}

// writeFile creates name, runs write on it, and returns the first write or close error.
func writeFile(name string, write func(*os.File) error) error {
	f, err := os.Create(name)
	if err != nil {
		return err
	}
	if err := write(f); err != nil {
		f.Close()
		return err
	}
	return f.Close()
}

// newRootCmd builds the golooper-sidd command with its flags; --sigma is signed (SIST's -s is a positive magnitude) and --max-runs is a go-looper extension absent from SIST.
func newRootCmd() *cobra.Command {
	o := &options{params: sidd.DefaultParams()}
	cmd := &cobra.Command{
		Use:     "golooper-sidd",
		Short:   "golooper-sidd computes SIDD strand-separation profiles for a single DNA sequence.",
		Version: "golooper-sidd v0.5.0",
		Args:    cobra.NoArgs,
		RunE: func(cmd *cobra.Command, args []string) error {
			cmd.SilenceUsage = true
			return o.run(cmd)
		},
	}
	cmd.SetVersionTemplate("{{.Version}}\n")
	f := cmd.Flags()
	f.StringVarP(&o.infile, "infile", "f", "", "input FASTA file with exactly one record (required)")
	f.StringVarP(&o.outfile, "outfile", "o", "", "output report file (default stdout)")
	f.Float64Var(&o.params.Sigma, "sigma", o.params.Sigma, "superhelical density, signed; SIST's positive-magnitude convention maps to a negative value here")
	f.Float64Var(&o.params.Temperature, "temperature", o.params.Temperature, "temperature in kelvin")
	f.Float64Var(&o.params.Salt, "salt", o.params.Salt, "salt concentration in molar")
	f.Float64Var(&o.params.Threshold, "threshold", o.params.Threshold, "free energy threshold in kcal/mol")
	f.IntVar(&o.params.MinLen, "minlength", o.params.MinLen, "minimum melted-run length in bases")
	f.IntVar(&o.params.MaxLen, "maxlength", o.params.MaxLen, "maximum melted-run length in bases")
	f.StringVar(&o.domain, "N", "auto", "superhelical domain size: \"auto\" uses the sequence length (the padded length under --clamp); an explicit positive integer is used verbatim with no +50 adjustment")
	f.StringVar(&o.maxRuns, "max-runs", strconv.Itoa(o.params.MaxRuns), "maximum number of simultaneous denatured runs per state, 1-4 (SIST's structural cap is 4 and SIST has no such parameter; this is a go-looper extension); 1 reproduces the v0.4 single-run model byte-for-byte")
	f.BoolVar(&o.params.NearestNeighbor, "nearest-neighbor", false, "use nearest-neighbor energetics (requires --clamp or --circular)")
	f.BoolVar(&o.circular, "circular", false, "treat the sequence as a circular molecule")
	f.BoolVar(&o.clamp, "clamp", false, "append a 50-G clamp to the linear sequence (SIST's default geometry)")
	f.StringVar(&o.tracks, "tracks", "", "write <prefix>_sidd_bpprob.wig and <prefix>_sidd_gx.wig")
	_ = cmd.MarkFlagRequired("infile")
	return cmd
}

// main runs golooper-sidd and exits with status 1 on any command error.
func main() {
	if err := newRootCmd().Execute(); err != nil {
		os.Exit(1)
	}
}
