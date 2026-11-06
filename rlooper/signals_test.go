package rlooper

import (
	"bufio"
	"math"
	"os"
	"slices"
	"strconv"
	"strings"
	"testing"
)

// loadOrientedFixture loads the synthetic200 fixture gene and applies default orientation.
func loadOrientedFixture(t *testing.T) *Gene {
	t.Helper()
	genes := readFixture(t, "../testdata/fixtures/synthetic200.fa")
	g := genes[0]
	g.ApplyOrientation(false, false, false)
	return g
}

// readWigValues parses the data lines of a wiggle file, skipping track, browser, comment and fixedStep header lines.
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

// compareWig asserts that got matches the wiggle file at path element-wise within checkTol.
func compareWig(t *testing.T, got []float64, path string) {
	t.Helper()
	want := readWigValues(t, path)
	if len(got) != len(want) {
		t.Fatalf("length %d vs %d", len(got), len(want))
	}
	for i := range got {
		checkTol(t, got[i], want[i])
	}
}

// compareMFEWig asserts the Go MFE track marks the ensemble MFE and equals the C++ golden span, accepting a different golden window only if its energy matches within the golden tolerance (1e-9 absolute / 1e-5 relative).
func compareMFEWig(t *testing.T, seq []byte, m *EnergyModel, ens EnsembleResult, got []float64, path string) {
	t.Helper()
	gotStart, gotEnd := onesSpan(t, got)
	if gotStart != ens.MFE.Start || gotEnd != ens.MFE.End {
		t.Fatalf("Go MFE signal marks [%d,%d), want [%d,%d)", gotStart, gotEnd, ens.MFE.Start, ens.MFE.End)
	}
	start, end := onesSpan(t, readWigValues(t, path))
	if start == ens.MFE.Start && end == ens.MFE.End {
		return
	}
	checkTol(t, windowEnergy(m, prefixBP(seq, m, false), start, end-start+1), ens.MFE.Energy)
	t.Logf("%s: tied MFE, C++ chose window (%d,%d), Go chose (%d,%d)", path, start, end, ens.MFE.Start, ens.MFE.End)
}

// onesSpan returns the half-open span of the contiguous run of 1s in signal and fails if any other element is nonzero.
func onesSpan(t *testing.T, signal []float64) (start, end int) {
	t.Helper()
	start = slices.Index(signal, 1)
	if start < 0 {
		t.Fatal("signal has no MFE span")
	}
	end = start
	for end < len(signal) && signal[end] == 1 {
		end++
	}
	for i, v := range signal {
		if (i < start || i >= end) && v != 0 {
			t.Fatalf("signal[%d] = %g outside contiguous MFE span [%d,%d)", i, v, start, end)
		}
	}
	return start, end
}

// TestSignalsAgainstGoldenWigDefaults pins BPProb, AvgG and MFE tracks against the C++ golden wigs for default parameters.
func TestSignalsAgainstGoldenWigDefaults(t *testing.T) {
	g := loadOrientedFixture(t)
	m := NewDefaultModel()
	w := Windows{MinLength: 2, Workers: 1}
	ens := ComputeEnsemble(g.Sequence, m, w)
	sig := ComputeSignals(g.Sequence, m, w, ens, true)
	compareWig(t, sig.BPProb, "../testdata/golden/defaults/out_bpprob.wig")
	compareWig(t, sig.AvgG, "../testdata/golden/defaults/out_avgG.wig")
	compareMFEWig(t, g.Sequence, m, ens, sig.MFE, "../testdata/golden/defaults/out_mfe.wig")
}

// TestSignalsAgainstGoldenWigSigma010 pins BPProb and MFE tracks against the C++ golden wigs at sigma -0.10.
func TestSignalsAgainstGoldenWigSigma010(t *testing.T) {
	g := loadOrientedFixture(t)
	m := NewDefaultModel()
	m.SetSigma(-0.10)
	w := Windows{MinLength: 2, Workers: 1}
	ens := ComputeEnsemble(g.Sequence, m, w)
	sig := ComputeSignals(g.Sequence, m, w, ens, false)
	compareWig(t, sig.BPProb, "../testdata/golden/sigma010/out_bpprob.wig")
	compareMFEWig(t, g.Sequence, m, ens, sig.MFE, "../testdata/golden/sigma010/out_mfe.wig")
}

// TestSignalsAgainstGoldenWigNauto pins BPProb and MFE tracks against the C++ golden wigs at N=200.
func TestSignalsAgainstGoldenWigNauto(t *testing.T) {
	g := loadOrientedFixture(t)
	m := NewDefaultModel()
	m.SetN(200)
	w := Windows{MinLength: 2, Workers: 1}
	ens := ComputeEnsemble(g.Sequence, m, w)
	sig := ComputeSignals(g.Sequence, m, w, ens, false)
	compareWig(t, sig.BPProb, "../testdata/golden/nauto/out_bpprob.wig")
	compareMFEWig(t, g.Sequence, m, ens, sig.MFE, "../testdata/golden/nauto/out_mfe.wig")
}

// TestSignalsAgainstGoldenWigHomopolymer05 pins BPProb and MFE tracks against the C++ golden wigs in homopolymer mode with energy 0.5.
func TestSignalsAgainstGoldenWigHomopolymer05(t *testing.T) {
	g := loadOrientedFixture(t)
	m := NewDefaultModel()
	m.Homopolymer, m.OverrideEnergy = true, 0.5
	w := Windows{MinLength: 2, Workers: 1}
	ens := ComputeEnsemble(g.Sequence, m, w)
	sig := ComputeSignals(g.Sequence, m, w, ens, false)
	compareWig(t, sig.BPProb, "../testdata/golden/homopolymer05/out_bpprob.wig")
	compareMFEWig(t, g.Sequence, m, ens, sig.MFE, "../testdata/golden/homopolymer05/out_mfe.wig")
}

// TestSignalsAgainstNaiveDirectSum pins BPProb and AvgG against a brute-force per-window sum on GATTACA, independent of the running-sum implementation.
func TestSignalsAgainstNaiveDirectSum(t *testing.T) {
	m := NewDefaultModel()
	seq := []byte("GATTACA")
	w := Windows{MinLength: 2, Workers: 1}
	ens := ComputeEnsemble(seq, m, w)
	sig := ComputeSignals(seq, m, w, ens, true)

	p := prefixBP(seq, m, false)
	rt := m.RT()
	wantBP := make([]float64, len(seq))
	wantNum := make([]float64, len(seq))
	wantDen := make([]float64, len(seq))
	for s := 0; s < len(seq); s++ {
		for e := s + 1; e < len(seq); e++ {
			g := windowEnergy(m, p, s, e-s+1)
			prob := math.Exp(-g/rt - ens.LogZ)
			for i := s; i < e; i++ {
				wantBP[i] += prob
				wantNum[i] += prob * g
				wantDen[i] += prob
			}
		}
	}
	for i := range seq {
		if math.Abs(sig.BPProb[i]-wantBP[i]) > 1e-12 {
			t.Fatalf("BPProb[%d] = %g, want %g", i, sig.BPProb[i], wantBP[i])
		}
		wantAvgG := 0.0
		if wantDen[i] > 0 {
			wantAvgG = wantNum[i] / wantDen[i]
		}
		if math.Abs(sig.AvgG[i]-wantAvgG) > 1e-12 {
			t.Fatalf("AvgG[%d] = %g, want %g", i, sig.AvgG[i], wantAvgG)
		}
	}
}

// loadPfc53Fixture loads and orients the first record of the pfc53 sequence.
func loadPfc53Fixture(t *testing.T) *Gene {
	t.Helper()
	genes := readFixture(t, "../res/pfc53_full.fa")
	g := genes[0]
	g.ApplyOrientation(false, false, false)
	return g
}

// TestSignalsPfc53NoNegativeBPProb pins that BPProb is never negative and the final base has exactly zero BPProb and AvgG.
func TestSignalsPfc53NoNegativeBPProb(t *testing.T) {
	g := loadPfc53Fixture(t)
	m := NewDefaultModel()
	w := Windows{MinLength: 2, Workers: 1}
	ens := ComputeEnsemble(g.Sequence, m, w)
	sig := ComputeSignals(g.Sequence, m, w, ens, true)

	n := len(sig.BPProb)
	if sig.BPProb[n-1] != 0 {
		t.Fatalf("BPProb[%d] (final base) = %g, want exactly 0", n-1, sig.BPProb[n-1])
	}
	for i, v := range sig.BPProb {
		if v < 0 {
			t.Fatalf("BPProb[%d] = %g, want >= 0", i, v)
		}
	}
	if sig.AvgG[n-1] != 0 {
		t.Fatalf("AvgG[%d] (final base, uncovered) = %g, want exactly 0", n-1, sig.AvgG[n-1])
	}
}

// TestSignalsZeroStructures pins that all signals are zero when MinLength admits no windows.
func TestSignalsZeroStructures(t *testing.T) {
	m := NewDefaultModel()
	seq := []byte("ACGT")
	w := Windows{MinLength: 10, Workers: 1}
	ens := ComputeEnsemble(seq, m, w)
	sig := ComputeSignals(seq, m, w, ens, true)
	for i, v := range sig.BPProb {
		if v != 0 {
			t.Fatalf("BPProb[%d] = %g, want 0", i, v)
		}
	}
	for i, v := range sig.AvgG {
		if v != 0 {
			t.Fatalf("AvgG[%d] = %g, want 0", i, v)
		}
	}
	for i, v := range sig.MFE {
		if v != 0 {
			t.Fatalf("MFE[%d] = %g, want 0", i, v)
		}
	}
}

// TestComputeResidualsAgainstGolden pins all five ensemble residuals against the C++ golden stdout.
func TestComputeResidualsAgainstGolden(t *testing.T) {
	g := loadOrientedFixture(t)
	m := NewDefaultModel()
	w := Windows{MinLength: 2, Workers: 1}
	ens := ComputeEnsemble(g.Sequence, m, w)
	res := ComputeResiduals(g.Sequence, m, w, ens)
	want := parseResidualsGolden(t, "../testdata/golden/defaults/stdout.txt")
	checkTol(t, res.Twist, want["ensemble_residual_twist"])
	checkTol(t, res.LinkingDifference, want["ensemble_residual_linking_difference"])
	checkTol(t, res.Superhelicity, want["ensemble_residual_superhelicity"])
	checkTol(t, res.WrappingAbsorption, want["ensemble_wrapping_absorption"])
	checkTol(t, res.StrandSeparation, want["ensemble_strand_separation_absorption"])
}

// TestComputeResidualsZeroStructures pins that with no windows the residuals reduce to the ground-state terms and zero wrapping and strand separation.
func TestComputeResidualsZeroStructures(t *testing.T) {
	m := NewDefaultModel()
	seq := []byte("ACGT")
	w := Windows{MinLength: 10, Workers: 1}
	ens := ComputeEnsemble(seq, m, w)
	res := ComputeResiduals(seq, m, w, ens)
	wantTwist := m.K * m.C * m.Alpha / (2 * math.Pi)
	if res.Twist != wantTwist {
		t.Fatalf("Twist = %g, want %g", res.Twist, wantTwist)
	}
	if res.LinkingDifference != m.Alpha {
		t.Fatalf("LinkingDifference = %g, want %g", res.LinkingDifference, m.Alpha)
	}
	if rel(res.Superhelicity, m.Sigma) > 1e-12 {
		t.Fatalf("Superhelicity = %g, want %g", res.Superhelicity, m.Sigma)
	}
	if res.WrappingAbsorption != 0 {
		t.Fatalf("WrappingAbsorption = %g, want 0", res.WrappingAbsorption)
	}
	if res.StrandSeparation != 0 {
		t.Fatalf("StrandSeparation = %g, want 0", res.StrandSeparation)
	}
}

// parseResidualsGolden extracts the five ensemble residual values from a C++ golden stdout file and fails if any is missing.
func parseResidualsGolden(t *testing.T, path string) map[string]float64 {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()

	want := make(map[string]float64)
	keys := []string{
		"ensemble_residual_twist",
		"ensemble_residual_linking_difference",
		"ensemble_residual_superhelicity",
		"ensemble_wrapping_absorption",
		"ensemble_strand_separation_absorption",
	}
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		key, rest, ok := strings.Cut(sc.Text(), ": ")
		if !ok || !slices.Contains(keys, key) {
			continue
		}
		v, err := strconv.ParseFloat(strings.TrimSpace(rest), 64)
		if err != nil {
			t.Fatal(err)
		}
		want[key] = v
	}
	if err := sc.Err(); err != nil {
		t.Fatal(err)
	}
	for _, k := range keys {
		if _, ok := want[k]; !ok {
			t.Fatalf("missing %q in %s", k, path)
		}
	}
	return want
}
