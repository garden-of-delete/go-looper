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

func loadOrientedFixture(t *testing.T) *Gene {
	t.Helper()
	genes := readFixture(t, "../testdata/fixtures/synthetic200.fa")
	g := genes[0]
	g.ApplyOrientation(false, false, false)
	return g
}

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
		t.Fatalf("length %d vs %d", len(got), len(want))
	}
	for i := range got {
		checkTol(t, got[i], want[i])
	}
}

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
	checkTol(t, windowEnergy(m, PrefixBPEnergies(seq, m), start, end), ens.MFE.Energy)
	t.Logf("%s: tied MFE, C++ chose window (%d,%d), Go chose (%d,%d)", path, start, end, ens.MFE.Start, ens.MFE.End)
}

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

func TestSignalsAgainstGoldenWigDefaults(t *testing.T) {
	g := loadOrientedFixture(t)
	m := NewDefaultModel()
	ens := ComputeEnsemble(g.Sequence, m, 2, 1)
	sig := ComputeSignals(g.Sequence, m, 2, ens, true, 1)
	compareWig(t, sig.BPProb, "../testdata/golden/defaults/out_bpprob.wig")
	compareWig(t, sig.AvgG, "../testdata/golden/defaults/out_avgG.wig")
	compareMFEWig(t, g.Sequence, m, ens, sig.MFE, "../testdata/golden/defaults/out_mfe.wig")
}

func TestSignalsAgainstGoldenWigSigma010(t *testing.T) {
	g := loadOrientedFixture(t)
	m := NewDefaultModel()
	m.SetSigma(-0.10)
	ens := ComputeEnsemble(g.Sequence, m, 2, 1)
	sig := ComputeSignals(g.Sequence, m, 2, ens, false, 1)
	compareWig(t, sig.BPProb, "../testdata/golden/sigma010/out_bpprob.wig")
	compareMFEWig(t, g.Sequence, m, ens, sig.MFE, "../testdata/golden/sigma010/out_mfe.wig")
}

func TestSignalsAgainstGoldenWigNauto(t *testing.T) {
	g := loadOrientedFixture(t)
	m := NewDefaultModel()
	m.SetN(200)
	ens := ComputeEnsemble(g.Sequence, m, 2, 1)
	sig := ComputeSignals(g.Sequence, m, 2, ens, false, 1)
	compareWig(t, sig.BPProb, "../testdata/golden/nauto/out_bpprob.wig")
	compareMFEWig(t, g.Sequence, m, ens, sig.MFE, "../testdata/golden/nauto/out_mfe.wig")
}

func TestSignalsAgainstGoldenWigHomopolymer05(t *testing.T) {
	g := loadOrientedFixture(t)
	m := NewDefaultModel()
	m.Homopolymer, m.OverrideEnergy = true, 0.5
	ens := ComputeEnsemble(g.Sequence, m, 2, 1)
	sig := ComputeSignals(g.Sequence, m, 2, ens, false, 1)
	compareWig(t, sig.BPProb, "../testdata/golden/homopolymer05/out_bpprob.wig")
	compareMFEWig(t, g.Sequence, m, ens, sig.MFE, "../testdata/golden/homopolymer05/out_mfe.wig")
}

func TestSignalsAgainstNaiveDirectSum(t *testing.T) {
	m := NewDefaultModel()
	seq := []byte("GATTACA")
	ens := ComputeEnsemble(seq, m, 2, 1)
	sig := ComputeSignals(seq, m, 2, ens, true, 1)

	p := PrefixBPEnergies(seq, m)
	rt := m.RT()
	wantBP := make([]float64, len(seq))
	wantNum := make([]float64, len(seq))
	wantDen := make([]float64, len(seq))
	for s := 0; s < len(seq); s++ {
		for e := s + 1; e < len(seq); e++ {
			g := windowEnergy(m, p, s, e)
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

func loadPfc53Fixture(t *testing.T) *Gene {
	t.Helper()
	genes := readFixture(t, "../res/pfc53_full.fa")
	g := genes[0]
	g.ApplyOrientation(false, false, false)
	return g
}

func TestSignalsPfc53NoNegativeBPProb(t *testing.T) {
	g := loadPfc53Fixture(t)
	m := NewDefaultModel()
	ens := ComputeEnsemble(g.Sequence, m, 2, 1)
	sig := ComputeSignals(g.Sequence, m, 2, ens, true, 1)

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

func TestSignalsZeroStructures(t *testing.T) {
	m := NewDefaultModel()
	seq := []byte("ACGT")
	ens := ComputeEnsemble(seq, m, 10, 1)
	sig := ComputeSignals(seq, m, 10, ens, true, 1)
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

func TestComputeResidualsAgainstGolden(t *testing.T) {
	g := loadOrientedFixture(t)
	m := NewDefaultModel()
	ens := ComputeEnsemble(g.Sequence, m, 2, 1)
	res := ComputeResiduals(g.Sequence, m, 2, ens, 1)
	want := parseResidualsGolden(t, "../testdata/golden/defaults/stdout.txt")
	checkTol(t, res.Twist, want["ensemble_residual_twist"])
	checkTol(t, res.LinkingDifference, want["ensemble_residual_linking_difference"])
	checkTol(t, res.Superhelicity, want["ensemble_residual_superhelicity"])
	checkTol(t, res.WrappingAbsorption, want["ensemble_wrapping_absorption"])
	checkTol(t, res.StrandSeparation, want["ensemble_strand_separation_absorption"])
}

func TestComputeResidualsZeroStructures(t *testing.T) {
	m := NewDefaultModel()
	seq := []byte("ACGT")
	ens := ComputeEnsemble(seq, m, 10, 1)
	res := ComputeResiduals(seq, m, 10, ens, 1)
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
