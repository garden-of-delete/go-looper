package sim

import (
	"bufio"
	"bytes"
	"golooper/rlooper"
	"io"
	"math"
	"os"
	"reflect"
	"slices"
	"strconv"
	"strings"
	"testing"
)

func loadOrientedFixture(t *testing.T) *rlooper.Gene {
	t.Helper()
	f, err := os.Open("../testdata/fixtures/synthetic200.fa")
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	genes, err := rlooper.ReadGenes(f)
	if err != nil {
		t.Fatal(err)
	}
	g := genes[0]
	g.ApplyOrientation(false, false, false)
	return g
}

func checkTol(t *testing.T, got, want float64) {
	t.Helper()
	tol := math.Max(1e-9, 1e-5*math.Abs(want))
	if math.Abs(got-want) > tol {
		t.Fatalf("got %g want %g (tol %g)", got, want, tol)
	}
}

type dumpRow struct {
	energy, prob float64
}

func parseDump(t *testing.T, r io.Reader) map[[2]int64]dumpRow {
	t.Helper()
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64*1024), 1024*1024)
	rows := make(map[[2]int64]dumpRow)
	first := true
	for sc.Scan() {
		if first {
			first = false
			continue
		}
		fields := strings.Fields(sc.Text())
		if len(fields) != 4 {
			t.Fatalf("malformed dump line: %q", sc.Text())
		}
		start, err := strconv.ParseInt(fields[0], 10, 64)
		if err != nil {
			t.Fatal(err)
		}
		stop, err := strconv.ParseInt(fields[1], 10, 64)
		if err != nil {
			t.Fatal(err)
		}
		if start == 0 && stop == 0 {
			continue
		}
		energy, err := strconv.ParseFloat(fields[2], 64)
		if err != nil {
			t.Fatal(err)
		}
		prob, err := strconv.ParseFloat(fields[3], 64)
		if err != nil {
			t.Fatal(err)
		}
		rows[[2]int64{start, stop}] = dumpRow{energy: energy, prob: prob}
	}
	if err := sc.Err(); err != nil {
		t.Fatal(err)
	}
	return rows
}

func parseBedPeaks(t *testing.T, path string) []rlooper.Loci {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	var peaks []rlooper.Loci
	sc := bufio.NewScanner(f)
	lineNum := 0
	for sc.Scan() {
		lineNum++
		if lineNum <= 3 {
			continue
		}
		fields := strings.Split(sc.Text(), "\t")
		if len(fields) != 6 {
			t.Fatalf("malformed bed line: %q", sc.Text())
		}
		start, err := strconv.ParseInt(fields[1], 10, 64)
		if err != nil {
			t.Fatal(err)
		}
		end, err := strconv.ParseInt(fields[2], 10, 64)
		if err != nil {
			t.Fatal(err)
		}
		peaks = append(peaks, rlooper.Loci{Chromosome: fields[0], Strand: fields[5], Start: start, End: end})
	}
	if err := sc.Err(); err != nil {
		t.Fatal(err)
	}
	return peaks
}

func peakWindows(peaks []rlooper.Loci) [][2]int64 {
	out := make([][2]int64, len(peaks))
	for i, p := range peaks {
		out[i] = [2]int64{p.Start, p.End}
	}
	return out
}

type topRow struct {
	start, stop                                           int64
	energy, prob, twist, linkingDifference, superhelicity float64
}

func parseTopRows(t *testing.T, r io.Reader) []topRow {
	t.Helper()
	var rows []topRow
	sc := bufio.NewScanner(r)
	for sc.Scan() {
		fields := strings.Fields(sc.Text())
		if len(fields) != 7 {
			continue
		}
		vals := make([]float64, 7)
		ok := true
		for i, fld := range fields {
			v, err := strconv.ParseFloat(fld, 64)
			if err != nil {
				ok = false
				break
			}
			vals[i] = v
		}
		if !ok {
			continue
		}
		rows = append(rows, topRow{
			start: int64(vals[0]), stop: int64(vals[1]),
			energy: vals[2], prob: vals[3],
			twist: vals[4], linkingDifference: vals[5], superhelicity: vals[6],
		})
	}
	if err := sc.Err(); err != nil {
		t.Fatal(err)
	}
	return rows
}

func parseResiduals(t *testing.T, r io.Reader) map[string]float64 {
	t.Helper()
	keys := []string{
		"ensemble_residual_twist",
		"ensemble_residual_linking_difference",
		"ensemble_residual_superhelicity",
		"ensemble_wrapping_absorption",
		"ensemble_strand_separation_absorption",
	}
	values := make(map[string]float64)
	sc := bufio.NewScanner(r)
	for sc.Scan() {
		key, rest, ok := strings.Cut(sc.Text(), ": ")
		if !ok || !slices.Contains(keys, key) {
			continue
		}
		v, err := strconv.ParseFloat(strings.TrimSpace(rest), 64)
		if err != nil {
			t.Fatal(err)
		}
		values[key] = v
	}
	if err := sc.Err(); err != nil {
		t.Fatal(err)
	}
	for _, k := range keys {
		if _, ok := values[k]; !ok {
			t.Fatalf("missing %q in residuals output", k)
		}
	}
	return values
}

func TestCallPeaksAbsoluteUnit(t *testing.T) {
	g := &rlooper.Gene{Loc: rlooper.Loci{Chromosome: "chrT", Strand: "+", Start: 100}}
	signal := []float64{0, 0.5, 0.5, 0, 0.2, 0, 0.9}
	got := CallPeaksAbsolute(signal, 0.1, g)
	want := []rlooper.Loci{
		{Chromosome: "chrT", Strand: "+", Start: 101, End: 103},
		{Chromosome: "chrT", Strand: "+", Start: 104, End: 105},
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("peaks = %+v, want %+v", got, want)
	}
}

func TestWriteWigSignalUnit(t *testing.T) {
	g := &rlooper.Gene{Name: "X", Loc: rlooper.Loci{Chromosome: "chrT", Start: 1000, End: 1003}}
	var buf bytes.Buffer
	if err := WriteWigSignal(&buf, g, []float64{0, 0.5, 1e-7, 1}); err != nil {
		t.Fatal(err)
	}
	want := "browser position chrT:1000-1003\n#X\nfixedStep chrom=chrT start=1000 step=1\n0\n0.5\n1e-07\n1\n"
	if buf.String() != want {
		t.Fatalf("got %q, want %q", buf.String(), want)
	}
}

func TestWriteBedPeaksNegativeStrand(t *testing.T) {
	g := &rlooper.Gene{Name: "X", Loc: rlooper.Loci{Chromosome: "chrT", Strand: "-", Start: 100, End: 110}}
	peaks := []rlooper.Loci{{Chromosome: "chrT", Strand: "-", Start: 101, End: 103}}
	var buf bytes.Buffer
	if err := WriteBedPeaks(&buf, g, peaks); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(buf.String(), "NEG0") {
		t.Fatalf("output %q missing NEG0 label", buf.String())
	}
}

func TestWriteDumpNegativeStrandFlip(t *testing.T) {
	g := &rlooper.Gene{Loc: rlooper.Loci{Strand: "-", Start: 10, End: 19}}
	recs := []rlooper.StructureRecord{{Start: 2, End: 5, Energy: 1, Prob: 0.5}}
	var buf bytes.Buffer
	if err := WriteDump(&buf, g, 0, recs); err != nil {
		t.Fatal(err)
	}
	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	if len(lines) != 3 {
		t.Fatalf("expected 3 lines, got %d: %q", len(lines), buf.String())
	}
	if !strings.HasPrefix(lines[2], "14 17 ") {
		t.Fatalf("row = %q, want prefix %q", lines[2], "14 17 ")
	}
}

func TestWriteDumpInvalidStrand(t *testing.T) {
	g := &rlooper.Gene{Loc: rlooper.Loci{Strand: ""}}
	var buf bytes.Buffer
	if err := WriteDump(&buf, g, 0, nil); err == nil {
		t.Fatal("expected error for unspecified strand")
	}
	if buf.Len() != 0 {
		t.Fatalf("expected no output on error, got %q", buf.String())
	}
}

func TestPrintTopNegativeStrandFlip(t *testing.T) {
	g := &rlooper.Gene{Loc: rlooper.Loci{Strand: "-", Start: 10, End: 19}, Sequence: make([]byte, 50)}
	recs := []rlooper.StructureRecord{{Start: 2, End: 5, Energy: 1, Prob: 0.5}}
	var buf bytes.Buffer
	if err := PrintTop(&buf, g, rlooper.NewDefaultModel(), recs, 1, false); err != nil {
		t.Fatal(err)
	}
	if !strings.HasPrefix(buf.String(), "14 17 ") {
		t.Fatalf("row = %q, want prefix %q", buf.String(), "14 17 ")
	}
}

func TestLargeCoordinatesPrintAsIntegers(t *testing.T) {
	g := &rlooper.Gene{
		Name:     "BIG",
		Sequence: make([]byte, 200),
		Loc:      rlooper.Loci{Chromosome: "chr1", Strand: "+", Start: 123456789, End: 123456988},
	}
	recs := []rlooper.StructureRecord{{Start: 0, End: 199, Energy: 1, Prob: 0.5}}

	var wigBuf bytes.Buffer
	if err := WriteWigSignal(&wigBuf, g, []float64{0.5}); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(wigBuf.String(), "start=123456789") {
		t.Fatalf("wig output missing exact integer coordinate: %q", wigBuf.String())
	}

	peaks := []rlooper.Loci{{Chromosome: "chr1", Strand: "+", Start: 123456789, End: 123456988}}
	var bedBuf bytes.Buffer
	if err := WriteBedPeaks(&bedBuf, g, peaks); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(bedBuf.String(), "123456789\t123456988") {
		t.Fatalf("bed output missing exact integer coordinates: %q", bedBuf.String())
	}

	var dumpBuf bytes.Buffer
	if err := WriteDump(&dumpBuf, g, 0, recs); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(dumpBuf.String(), "123456789 123456988") {
		t.Fatalf("dump output missing exact integer coordinates: %q", dumpBuf.String())
	}

	m := rlooper.NewDefaultModel()
	var topBuf bytes.Buffer
	if err := PrintTop(&topBuf, g, m, recs, 1, false); err != nil {
		t.Fatal(err)
	}
	if !strings.HasPrefix(topBuf.String(), "123456789 123456988 ") {
		t.Fatalf("top output missing exact integer coordinates: %q", topBuf.String())
	}

	for name, out := range map[string]string{"wig": wigBuf.String(), "bed": bedBuf.String(), "dump": dumpBuf.String(), "top": topBuf.String()} {
		if strings.Contains(out, "e+") {
			t.Fatalf("%s output used scientific notation for a coordinate: %q", name, out)
		}
	}
}

func TestPrintResidualsUnit(t *testing.T) {
	r := rlooper.Residuals{
		Twist:              -0.200898,
		LinkingDifference:  -2.51641,
		Superhelicity:      -0.0174471,
		WrappingAbsorption: -1.83643,
		StrandSeparation:   -5.74331,
	}
	var buf bytes.Buffer
	if err := PrintResiduals(&buf, r); err != nil {
		t.Fatal(err)
	}
	want := "ensemble_residual_twist: -0.200898\n" +
		"ensemble_residual_linking_difference: -2.51641\n" +
		"ensemble_residual_superhelicity: -0.0174471\n" +
		"ensemble_wrapping_absorption: -1.83643\n" +
		"ensemble_strand_separation_absorption: -5.74331\n"
	if buf.String() != want {
		t.Fatalf("got %q, want %q", buf.String(), want)
	}
}

func TestWriteDumpAgainstGoldenDefaults(t *testing.T) {
	g := loadOrientedFixture(t)
	m := rlooper.NewDefaultModel()
	ens := rlooper.ComputeEnsemble(g.Sequence, m, 2, 1)
	recs := rlooper.AllStructures(g.Sequence, m, 2, ens)

	var buf bytes.Buffer
	if err := WriteDump(&buf, g, m.GroundStateEnergy(), recs); err != nil {
		t.Fatal(err)
	}

	lines := strings.SplitN(buf.String(), "\n", 3)
	if lines[0] != "start_position stop_position energy probability" {
		t.Fatalf("header = %q", lines[0])
	}
	if lines[1] != "0 0 36.018 0" {
		t.Fatalf("ground row = %q, want %q", lines[1], "0 0 36.018 0")
	}

	got := parseDump(t, strings.NewReader(buf.String()))

	goldenFile, err := os.Open("../testdata/golden/defaults/outSYNTH200_dump.txt")
	if err != nil {
		t.Fatal(err)
	}
	defer goldenFile.Close()
	want := parseDump(t, goldenFile)

	if len(got) != len(want) {
		t.Fatalf("row count got %d, want %d", len(got), len(want))
	}
	for key, w := range want {
		g2, ok := got[key]
		if !ok {
			t.Fatalf("missing row for %v", key)
		}
		checkTol(t, g2.energy, w.energy)
		checkTol(t, g2.prob, w.prob)
	}
}

func TestCallPeaksAbsoluteAgainstGoldenDefaults(t *testing.T) {
	g := loadOrientedFixture(t)
	m := rlooper.NewDefaultModel()
	ens := rlooper.ComputeEnsemble(g.Sequence, m, 2, 1)
	sig := rlooper.ComputeSignals(g.Sequence, m, 2, ens, true, 1)

	bpPeaks := CallPeaksAbsolute(sig.BPProb, 0.1, g)
	wantBP := parseBedPeaks(t, "../testdata/golden/defaults/out_bpprob.bed")
	if !reflect.DeepEqual(peakWindows(bpPeaks), peakWindows(wantBP)) {
		t.Fatalf("bpprob peaks = %v, want %v", peakWindows(bpPeaks), peakWindows(wantBP))
	}

	mfePeaks := CallPeaksAbsolute(sig.MFE, 0.1, g)
	wantMFE := parseBedPeaks(t, "../testdata/golden/defaults/out_mfe.bed")
	if !reflect.DeepEqual(peakWindows(mfePeaks), peakWindows(wantMFE)) {
		t.Fatalf("mfe peaks = %v, want %v", peakWindows(mfePeaks), peakWindows(wantMFE))
	}
}

func TestPrintTopAgainstGoldenDefaults(t *testing.T) {
	g := loadOrientedFixture(t)
	m := rlooper.NewDefaultModel()
	ens := rlooper.ComputeEnsemble(g.Sequence, m, 2, 1)
	recs := rlooper.AllStructures(g.Sequence, m, 2, ens)

	recsByCoord := make(map[[2]int]float64, len(recs))
	for _, r := range recs {
		recsByCoord[[2]int{r.Start, r.End}] = r.Energy
	}

	var buf bytes.Buffer
	if err := PrintTop(&buf, g, m, recs, 5, true); err != nil {
		t.Fatal(err)
	}
	got := parseTopRows(t, strings.NewReader(buf.String()))

	goldenFile, err := os.Open("../testdata/golden/defaults/stdout.txt")
	if err != nil {
		t.Fatal(err)
	}
	defer goldenFile.Close()
	want := parseTopRows(t, goldenFile)

	if len(got) != 5 || len(want) != 5 {
		t.Fatalf("expected 5 top rows, got %d, want %d", len(got), len(want))
	}

	for i := range want {
		w := want[i]
		g2 := got[i]
		checkTol(t, g2.energy, w.energy)
		checkTol(t, g2.prob, w.prob)
		checkTol(t, g2.twist, w.twist)
		checkTol(t, g2.linkingDifference, w.linkingDifference)
		checkTol(t, g2.superhelicity, w.superhelicity)

		if g2.start == w.start && g2.stop == w.stop {
			continue
		}
		relStart := int(w.start - g.Loc.Start)
		relStop := int(w.stop - g.Loc.Start)
		energy, ok := recsByCoord[[2]int{relStart, relStop}]
		if !ok {
			t.Fatalf("row %d: Go coords (%d,%d) vs golden coords (%d,%d) disagree, and golden window is not a Go structure at all", i, g2.start, g2.stop, w.start, w.stop)
		}
		checkTol(t, energy, g2.energy)
		t.Logf("row %d: tie between Go (%d,%d) and golden (%d,%d) at energy %g", i, g2.start, g2.stop, w.start, w.stop, g2.energy)
	}
}

func TestPrintResidualsAgainstGoldenDefaults(t *testing.T) {
	g := loadOrientedFixture(t)
	m := rlooper.NewDefaultModel()
	ens := rlooper.ComputeEnsemble(g.Sequence, m, 2, 1)
	res := rlooper.ComputeResiduals(g.Sequence, m, 2, ens)

	var buf bytes.Buffer
	if err := PrintResiduals(&buf, res); err != nil {
		t.Fatal(err)
	}
	got := parseResiduals(t, strings.NewReader(buf.String()))

	goldenFile, err := os.Open("../testdata/golden/defaults/stdout.txt")
	if err != nil {
		t.Fatal(err)
	}
	defer goldenFile.Close()
	want := parseResiduals(t, goldenFile)

	for k, w := range want {
		checkTol(t, got[k], w)
	}
}
