package sim

import (
	"fmt"
	"golooper/rlooper"
	"io"
	"strconv"
	"strings"
)

func formatValue(v float64) string {
	return strconv.FormatFloat(v, 'g', 6, 64)
}

func formatRow(vals ...float64) string {
	parts := make([]string, len(vals))
	for i, v := range vals {
		parts[i] = formatValue(v)
	}
	return strings.Join(parts, " ")
}

func WriteWigSignal(w io.Writer, g *rlooper.Gene, signal []float64) error {
	if _, err := fmt.Fprintf(w, "browser position %s:%d-%d\n#%s\nfixedStep chrom=%s start=%d step=1\n",
		g.Loc.Chromosome, g.Loc.Start, g.Loc.End, g.Name, g.Loc.Chromosome, g.Loc.Start); err != nil {
		return err
	}
	for _, v := range signal {
		if _, err := fmt.Fprintln(w, formatValue(v)); err != nil {
			return err
		}
	}
	return nil
}

func CallPeaksAbsolute(signal []float64, threshold float64, g *rlooper.Gene) []rlooper.Loci {
	var peaks []rlooper.Loci
	inPeak := false
	var peakStart int64

	for i, v := range signal {
		switch {
		case v > threshold && !inPeak:
			inPeak = true
			peakStart = g.Loc.Start + int64(i)
		case v <= threshold && inPeak:
			inPeak = false
			peaks = append(peaks, rlooper.Loci{
				Chromosome: g.Loc.Chromosome,
				Strand:     g.Loc.Strand,
				Start:      peakStart,
				End:        g.Loc.Start + int64(i),
			})
		}
	}

	return peaks
}

func WriteBedPeaks(w io.Writer, g *rlooper.Gene, peaks []rlooper.Loci) error {
	if _, err := fmt.Fprintf(w, "browser position %s:%d-%d\n#%s\n", g.Loc.Chromosome, g.Loc.Start, g.Loc.End, g.Name); err != nil {
		return err
	}

	name := "NEG"
	if g.Loc.Strand == "+" {
		name = "POS"
	}

	for i, p := range peaks {
		if _, err := fmt.Fprintf(w, "%s\t%d\t%d\t%s%d\t0\t%s\n", p.Chromosome, p.Start, p.End, name, i, p.Strand); err != nil {
			return err
		}
	}

	return nil
}

func WriteDump(w io.Writer, g *rlooper.Gene, groundEnergy float64, recs []rlooper.StructureRecord) error {
	if g.Loc.Strand != "+" && g.Loc.Strand != "-" {
		return fmt.Errorf("dump error: strand unspecified")
	}

	if _, err := fmt.Fprintln(w, "start_position stop_position energy probability"); err != nil {
		return err
	}
	if _, err := fmt.Fprintf(w, "0 0 %s\n", formatRow(groundEnergy, 0)); err != nil {
		return err
	}

	for _, r := range recs {
		start, stop := genomicSpan(g, r)
		if _, err := fmt.Fprintf(w, "%d %d %s\n", start, stop, formatRow(r.Energy, r.Prob)); err != nil {
			return err
		}
	}

	return nil
}

func genomicSpan(g *rlooper.Gene, r rlooper.StructureRecord) (start, stop int64) {
	if g.Loc.Strand == "-" {
		return g.Loc.End - int64(r.End), g.Loc.End - int64(r.Start)
	}
	return g.Loc.Start + int64(r.Start), g.Loc.Start + int64(r.End)
}

func PrintTop(w io.Writer, g *rlooper.Gene, m *rlooper.EnergyModel, recs []rlooper.StructureRecord, n int, withResiduals bool) error {
	if n > len(recs) {
		n = len(recs)
	}

	for i := 0; i < n; i++ {
		r := recs[i]

		start, stop := genomicSpan(g, r)

		var twist, linkingDifference, superhelicity float64
		if withResiduals {
			twist, linkingDifference = rlooper.StructureResiduals(m, r)
			superhelicity = linkingDifference / (m.N * m.A)
		}

		residuals := formatRow(r.Energy, r.Prob, twist, linkingDifference, superhelicity)
		if _, err := fmt.Fprintf(w, "%d %d %s\n", start, stop, residuals); err != nil {
			return err
		}
	}

	return nil
}

func PrintResiduals(w io.Writer, r rlooper.Residuals) error {
	lines := []struct {
		label string
		value float64
	}{
		{"ensemble_residual_twist: ", r.Twist},
		{"ensemble_residual_linking_difference: ", r.LinkingDifference},
		{"ensemble_residual_superhelicity: ", r.Superhelicity},
		{"ensemble_wrapping_absorption: ", r.WrappingAbsorption},
		{"ensemble_strand_separation_absorption: ", r.StrandSeparation},
	}

	for _, l := range lines {
		if _, err := fmt.Fprintf(w, "%s%s\n", l.label, formatValue(l.value)); err != nil {
			return err
		}
	}

	return nil
}
