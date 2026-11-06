package sim

import (
	"fmt"
	"golooper/rlooper"
	"io"
	"strconv"
	"strings"
)

// formatValue renders v with six significant digits in %g style, the number format pinned byte-for-byte by unit tests.
func formatValue(v float64) string {
	return strconv.FormatFloat(v, 'g', 6, 64)
}

// formatRow renders vals with formatValue, space-separated, as one row of the dump, top and residual outputs.
func formatRow(vals ...float64) string {
	parts := make([]string, len(vals))
	for i, v := range vals {
		parts[i] = formatValue(v)
	}
	return strings.Join(parts, " ")
}

// WriteWigSignal writes a UCSC WIG fixedStep track for signal, taking chrom and start from g.Loc, so a FASTA without a range= header yields an empty chrom and start 0 unless the caller populates g.Loc first; format pinned byte-for-byte by unit tests.
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

type Peak struct {
	Spans []rlooper.Loci
}

// CallPeaks returns the maximal runs of signal strictly above threshold as genomic peaks, merging a run touching the sequence end with one starting at index 0 when circular.
func CallPeaks(signal []float64, threshold float64, g *rlooper.Gene, circular bool) []Peak {
	if circular {
		return callPeaksCircular(signal, threshold, g)
	}
	return callPeaksLinear(signal, threshold, g)
}

// scanAboveThresholdRuns returns the closed above-threshold runs plus the start and presence of a run still open at the end of signal.
func scanAboveThresholdRuns(signal []float64, threshold float64, g *rlooper.Gene) (peaks []Peak, openStart int64, open bool) {
	inPeak := false
	var peakStart int64

	for i, v := range signal {
		switch {
		case v > threshold && !inPeak:
			inPeak = true
			peakStart = g.Loc.Start + int64(i)
		case v <= threshold && inPeak:
			inPeak = false
			peaks = append(peaks, Peak{Spans: []rlooper.Loci{{
				Chromosome: g.Loc.Chromosome,
				Strand:     g.Loc.Strand,
				Start:      peakStart,
				End:        g.Loc.Start + int64(i),
			}}})
		}
	}

	return peaks, peakStart, inPeak
}

// callPeaksLinear calls peaks on a linear sequence, so a run reaching the end of signal is dropped, never closed.
func callPeaksLinear(signal []float64, threshold float64, g *rlooper.Gene) []Peak {
	peaks, _, _ := scanAboveThresholdRuns(signal, threshold, g)
	return peaks
}

// callPeaksCircular closes a trailing open run at the sequence end, joining it to a run starting at the gene start as one two-span peak across the origin; a trailing run with no run at index 0 becomes its own single-span peak.
func callPeaksCircular(signal []float64, threshold float64, g *rlooper.Gene) []Peak {
	peaks, openStart, open := scanAboveThresholdRuns(signal, threshold, g)
	if !open {
		return peaks
	}

	endArc := rlooper.Loci{
		Chromosome: g.Loc.Chromosome,
		Strand:     g.Loc.Strand,
		Start:      openStart,
		End:        g.Loc.Start + int64(len(signal)),
	}

	if len(peaks) > 0 && peaks[0].Spans[0].Start == g.Loc.Start {
		peaks[0].Spans = []rlooper.Loci{endArc, peaks[0].Spans[0]}
		return peaks
	}

	return append(peaks, Peak{Spans: []rlooper.Loci{endArc}})
}

// WriteBedPeaks writes peaks as BED rows named POS or NEG by gene strand plus peak index, with spans of one wrapped peak sharing a name; format pinned byte-for-byte by unit tests.
func WriteBedPeaks(w io.Writer, g *rlooper.Gene, peaks []Peak) error {
	if _, err := fmt.Fprintf(w, "browser position %s:%d-%d\n#%s\n", g.Loc.Chromosome, g.Loc.Start, g.Loc.End, g.Name); err != nil {
		return err
	}

	name := "NEG"
	if g.Loc.Strand == "+" {
		name = "POS"
	}

	for i, peak := range peaks {
		for _, span := range peak.Spans {
			if _, err := fmt.Fprintf(w, "%s\t%d\t%d\t%s%d\t0\t%s\n", span.Chromosome, span.Start, span.End, name, i, span.Strand); err != nil {
				return err
			}
		}
	}

	return nil
}

// WriteDump writes the header, a ground-state row and one row per structure with genomic coordinates, erroring before any output if the gene strand is not + or -.
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

// genomicSpan maps a structure's sequence-relative coordinates to genomic ones, mirroring them about Loc.End on the minus strand; a wrapped record yields start greater than stop.
func genomicSpan(g *rlooper.Gene, r rlooper.StructureRecord) (start, stop int64) {
	if g.Loc.Strand == "-" {
		return g.Loc.End - int64(r.End), g.Loc.End - int64(r.Start)
	}
	return g.Loc.Start + int64(r.Start), g.Loc.Start + int64(r.End)
}

// PrintTop writes the first n of recs as genomic start, stop, energy and probability rows, adding twist, linking difference and superhelicity when withResiduals is set.
func PrintTop(w io.Writer, g *rlooper.Gene, m *rlooper.EnergyModel, recs []rlooper.StructureRecord, n int, withResiduals bool) error {
	if n > len(recs) {
		n = len(recs)
	}

	for i := 0; i < n; i++ {
		r := recs[i]

		start, stop := genomicSpan(g, r)

		var twist, linkingDifference, superhelicity float64
		if withResiduals {
			twist, linkingDifference = rlooper.StructureResiduals(m, r.Length(len(g.Sequence)))
			superhelicity = linkingDifference / (m.N * m.A)
		}

		residuals := formatRow(r.Energy, r.Prob, twist, linkingDifference, superhelicity)
		if _, err := fmt.Fprintf(w, "%d %d %s\n", start, stop, residuals); err != nil {
			return err
		}
	}

	return nil
}

// PrintResiduals writes the five labelled ensemble residual lines to w in the fixed order and labels pinned by the goldens.
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
