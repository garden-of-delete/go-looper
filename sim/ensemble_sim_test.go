package sim

import (
	"bufio"
	"bytes"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"testing"

	"golooper/config"
	"golooper/rlooper"
)

func readWigValues(t *testing.T, path string) []float64 {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	var values []float64
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "track") || strings.HasPrefix(line, "browser") ||
			strings.HasPrefix(line, "#") || strings.HasPrefix(line, "fixedStep") {
			continue
		}
		v, err := strconv.ParseFloat(line, 64)
		if err != nil {
			t.Fatal(err)
		}
		values = append(values, v)
	}
	if err := sc.Err(); err != nil {
		t.Fatal(err)
	}
	return values
}

func compareWig(t *testing.T, got []float64, path string) {
	t.Helper()
	want := readWigValues(t, path)
	if len(got) != len(want) {
		t.Fatalf("%s: length %d vs %d", path, len(got), len(want))
	}
	for i := range got {
		checkTol(t, got[i], want[i])
	}
}

func onesSpan(t *testing.T, signal []float64) (start, end int) {
	t.Helper()
	start = -1
	for i, v := range signal {
		if v == 1 {
			start = i
			break
		}
	}
	if start < 0 {
		t.Fatal("signal has no MFE span")
	}
	end = start
	for end < len(signal) && signal[end] == 1 {
		end++
	}
	return start, end
}

func compareMFEWig(t *testing.T, ens ensembleTieInfo, got []float64, path string) {
	t.Helper()
	gotStart, gotEnd := onesSpan(t, got)
	wantStart, wantEnd := onesSpan(t, readWigValues(t, path))
	if gotStart == wantStart && gotEnd == wantEnd {
		return
	}
	energies := make(map[[2]int64]float64, len(ens.recs))
	for _, r := range ens.recs {
		start, stop := genomicSpan(ens.gene, r)
		energies[[2]int64{start, stop}] = r.Energy
	}
	wantAbsStart := ens.gene.Loc.Start + int64(wantStart)
	wantAbsEnd := ens.gene.Loc.Start + int64(wantEnd)
	wantEnergy, ok := energies[[2]int64{wantAbsStart, wantAbsEnd}]
	if !ok {
		t.Fatalf("%s: golden MFE window (%d,%d) is not a Go structure at all, Go chose (%d,%d)", path, wantStart, wantEnd, gotStart, gotEnd)
	}
	checkTol(t, ens.mfeEnergy, wantEnergy)
	t.Logf("%s: tied MFE, C++ chose window (%d,%d), Go chose (%d,%d), energy %g", path, wantStart, wantEnd, gotStart, gotEnd, ens.mfeEnergy)
}

func parsePGroundState(t *testing.T, text string) float64 {
	t.Helper()
	for _, line := range strings.Split(text, "\n") {
		rest, ok := strings.CutPrefix(line, "P(ground state)= ")
		if !ok {
			continue
		}
		v, err := strconv.ParseFloat(strings.TrimSpace(rest), 64)
		if err != nil {
			t.Fatal(err)
		}
		return v
	}
	t.Fatal("P(ground state) line not found in stdout")
	return 0
}

func runEnsembleCase(t *testing.T, cfg *config.Config) (dir string, stdout string) {
	t.Helper()
	var buf bytes.Buffer
	if err := RunEnsembleSim(cfg, runtime.NumCPU(), &buf); err != nil {
		t.Fatalf("RunEnsembleSim: %v", err)
	}
	return filepath.Dir(cfg.OutfileName), buf.String()
}

func TestRunEnsembleSimAgainstGoldenDefaults(t *testing.T) {
	dir := t.TempDir()
	sensitivity := -1
	top := 5
	cfg := &config.Config{
		InfileName:         "../testdata/fixtures/synthetic200.fa",
		OutfileName:        filepath.Join(dir, "out"),
		LocalAverageEnergy: true,
		Dump:               true,
		Residuals:          true,
		Bedfile:            true,
		Sensitivity:        &sensitivity,
		Top:                &top,
	}

	_, stdout := runEnsembleCase(t, cfg)
	goldenDir := "../testdata/golden/defaults"

	if !strings.Contains(stdout, "processing gene: SYNTH200...19900") {
		t.Fatalf("stdout missing expected processing line: %q", stdout)
	}

	goldenStdout, err := os.ReadFile(filepath.Join(goldenDir, "stdout.txt"))
	if err != nil {
		t.Fatal(err)
	}
	checkTol(t, parsePGroundState(t, stdout), parsePGroundState(t, string(goldenStdout)))

	gotResiduals := parseResiduals(t, strings.NewReader(stdout))
	wantResiduals := parseResiduals(t, bytes.NewReader(goldenStdout))
	for k, w := range wantResiduals {
		checkTol(t, gotResiduals[k], w)
	}

	gotTop := parseTopRows(t, strings.NewReader(stdout))
	wantTop := parseTopRows(t, bytes.NewReader(goldenStdout))
	if len(gotTop) != len(wantTop) {
		t.Fatalf("top row count got %d, want %d", len(gotTop), len(wantTop))
	}

	ens := loadEnsembleForTie(t, cfg)
	recsByCoord := make(map[[2]int]float64, len(ens.recs))
	for _, r := range ens.recs {
		recsByCoord[[2]int{r.Start, r.End}] = r.Energy
	}
	for i := range wantTop {
		w := wantTop[i]
		g2 := gotTop[i]
		checkTol(t, g2.energy, w.energy)
		checkTol(t, g2.prob, w.prob)
		checkTol(t, g2.twist, w.twist)
		checkTol(t, g2.linkingDifference, w.linkingDifference)
		checkTol(t, g2.superhelicity, w.superhelicity)
		if g2.start == w.start && g2.stop == w.stop {
			continue
		}
		relStart := int(w.start - ens.gene.Loc.Start)
		relStop := int(w.stop - ens.gene.Loc.Start)
		energy, ok := recsByCoord[[2]int{relStart, relStop}]
		if !ok {
			t.Fatalf("top row %d: Go coords (%d,%d) vs golden coords (%d,%d) disagree, and golden window is not a Go structure at all", i, g2.start, g2.stop, w.start, w.stop)
		}
		checkTol(t, energy, g2.energy)
		t.Logf("top row %d: tie between Go (%d,%d) and golden (%d,%d) at energy %g", i, g2.start, g2.stop, w.start, w.stop, g2.energy)
	}

	compareWig(t, readWigValues(t, filepath.Join(dir, "out_bpprob.wig")), filepath.Join(goldenDir, "out_bpprob.wig"))
	compareWig(t, readWigValues(t, filepath.Join(dir, "out_avgG.wig")), filepath.Join(goldenDir, "out_avgG.wig"))
	compareMFEWig(t, ens, readWigValues(t, filepath.Join(dir, "out_mfe.wig")), filepath.Join(goldenDir, "out_mfe.wig"))

	wantBPPeaks := parseBedPeaks(t, filepath.Join(goldenDir, "out_bpprob.bed"))
	gotBPPeaks := parseBedPeaks(t, filepath.Join(dir, "out_bpprob.bed"))
	if !equalWindows(peakWindows(gotBPPeaks), peakWindows(wantBPPeaks)) {
		t.Fatalf("bpprob peaks = %v, want %v", peakWindows(gotBPPeaks), peakWindows(wantBPPeaks))
	}
	wantMFEPeaks := parseBedPeaks(t, filepath.Join(goldenDir, "out_mfe.bed"))
	gotMFEPeaks := parseBedPeaks(t, filepath.Join(dir, "out_mfe.bed"))
	if !equalWindows(peakWindows(gotMFEPeaks), peakWindows(wantMFEPeaks)) {
		t.Fatalf("mfe peaks = %v, want %v", peakWindows(gotMFEPeaks), peakWindows(wantMFEPeaks))
	}

	gotDump, err := os.Open(filepath.Join(dir, "outSYNTH200_dump.txt"))
	if err != nil {
		t.Fatal(err)
	}
	defer gotDump.Close()
	gotDumpRows := parseDump(t, gotDump)
	goldenDump, err := os.Open(filepath.Join(goldenDir, "outSYNTH200_dump.txt"))
	if err != nil {
		t.Fatal(err)
	}
	defer goldenDump.Close()
	wantDumpRows := parseDump(t, goldenDump)
	if len(gotDumpRows) != len(wantDumpRows) {
		t.Fatalf("dump row count got %d, want %d", len(gotDumpRows), len(wantDumpRows))
	}
	for key, w := range wantDumpRows {
		g2, ok := gotDumpRows[key]
		if !ok {
			t.Fatalf("missing dump row for %v", key)
		}
		checkTol(t, g2.energy, w.energy)
		checkTol(t, g2.prob, w.prob)
	}
}

func TestRunEnsembleSimAgainstGoldenSigma010(t *testing.T) {
	dir := t.TempDir()
	sigma := -0.10
	cfg := &config.Config{
		InfileName:          "../testdata/fixtures/synthetic200.fa",
		OutfileName:         filepath.Join(dir, "out"),
		SuperhelicalDensity: &sigma,
	}

	_, stdout := runEnsembleCase(t, cfg)
	goldenDir := "../testdata/golden/sigma010"

	goldenStdout, err := os.ReadFile(filepath.Join(goldenDir, "stdout.txt"))
	if err != nil {
		t.Fatal(err)
	}
	checkTol(t, parsePGroundState(t, stdout), parsePGroundState(t, string(goldenStdout)))

	ens := loadEnsembleForTie(t, cfg)
	compareWig(t, readWigValues(t, filepath.Join(dir, "out_bpprob.wig")), filepath.Join(goldenDir, "out_bpprob.wig"))
	compareMFEWig(t, ens, readWigValues(t, filepath.Join(dir, "out_mfe.wig")), filepath.Join(goldenDir, "out_mfe.wig"))
}

func TestRunEnsembleSimAgainstGoldenPfc53(t *testing.T) {
	dir := t.TempDir()
	cfg := &config.Config{
		InfileName:  "../res/pfc53_full.fa",
		OutfileName: filepath.Join(dir, "out"),
	}

	runEnsembleCase(t, cfg)
	goldenDir := "../testdata/golden/pfc53"

	ens := loadEnsembleForTie(t, cfg)
	compareWig(t, readWigValues(t, filepath.Join(dir, "out_bpprob.wig")), filepath.Join(goldenDir, "out_bpprob.wig"))
	compareMFEWig(t, ens, readWigValues(t, filepath.Join(dir, "out_mfe.wig")), filepath.Join(goldenDir, "out_mfe.wig"))
}

type ensembleTieInfo struct {
	gene      *rlooper.Gene
	recs      []rlooper.StructureRecord
	mfeEnergy float64
}

func loadEnsembleForTie(t *testing.T, cfg *config.Config) ensembleTieInfo {
	t.Helper()
	f, err := os.Open(cfg.InfileName)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	genes, err := rlooper.ReadGenes(f)
	if err != nil {
		t.Fatal(err)
	}
	g := genes[0]
	m := ModelFromConfig(cfg)
	if cfg.AutoDomainSize {
		m.SetN(float64(len(g.Sequence)))
	}
	complement := cfg.Complement != nil && *cfg.Complement
	reverse := cfg.Reverse != nil && *cfg.Reverse
	g.ApplyOrientation(complement, reverse, cfg.Invert)
	minLength := minLengthFromConfig(cfg)
	ens := rlooper.ComputeEnsemble(g.Sequence, m, minLength, 1)
	recs := rlooper.AllStructures(g.Sequence, m, minLength, ens)
	return ensembleTieInfo{gene: g, recs: recs, mfeEnergy: ens.MFE.Energy}
}

func equalWindows(a, b [][2]int64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}
