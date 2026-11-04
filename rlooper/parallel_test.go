package rlooper

import (
	"bytes"
	"fmt"
	"math"
	"os"
	"runtime"
	"slices"
	"testing"
)

func relOrAbs(got, want float64) float64 {
	if want == 0 {
		return math.Abs(got - want)
	}
	return rel(got, want)
}

func withinSignalTol(got, want float64) bool {
	return math.Abs(got-want) <= math.Max(1e-12, 1e-12*math.Abs(want))
}

func assertEnsembleEquivalent(t *testing.T, p, s EnsembleResult) {
	t.Helper()
	if rel(p.LogZ, s.LogZ) > 1e-12 {
		t.Fatalf("LogZ %g vs %g", p.LogZ, s.LogZ)
	}
	if p.NumStructures != s.NumStructures {
		t.Fatalf("NumStructures %d vs %d", p.NumStructures, s.NumStructures)
	}
	if p.MFE.Start != s.MFE.Start || p.MFE.End != s.MFE.End || p.MFE.Energy != s.MFE.Energy {
		t.Fatalf("MFE %+v vs %+v", p.MFE, s.MFE)
	}
	if relOrAbs(p.MFE.Prob, s.MFE.Prob) > 1e-12 {
		t.Fatalf("MFE.Prob %g vs %g", p.MFE.Prob, s.MFE.Prob)
	}
}

func assertSignalsEquivalent(t *testing.T, p, s Signals) {
	t.Helper()
	for i := range s.BPProb {
		if !withinSignalTol(p.BPProb[i], s.BPProb[i]) {
			t.Fatalf("BPProb[%d] %g vs %g", i, p.BPProb[i], s.BPProb[i])
		}
	}
	for i := range s.AvgG {
		if !withinSignalTol(p.AvgG[i], s.AvgG[i]) {
			t.Fatalf("AvgG[%d] %g vs %g", i, p.AvgG[i], s.AvgG[i])
		}
	}
}

func assertResidualsEquivalent(t *testing.T, p, s Residuals) {
	t.Helper()
	if relOrAbs(p.Twist, s.Twist) > 1e-12 {
		t.Fatalf("Twist %g vs %g", p.Twist, s.Twist)
	}
	if relOrAbs(p.LinkingDifference, s.LinkingDifference) > 1e-12 {
		t.Fatalf("LinkingDifference %g vs %g", p.LinkingDifference, s.LinkingDifference)
	}
	if relOrAbs(p.Superhelicity, s.Superhelicity) > 1e-12 {
		t.Fatalf("Superhelicity %g vs %g", p.Superhelicity, s.Superhelicity)
	}
	if relOrAbs(p.WrappingAbsorption, s.WrappingAbsorption) > 1e-12 {
		t.Fatalf("WrappingAbsorption %g vs %g", p.WrappingAbsorption, s.WrappingAbsorption)
	}
	if relOrAbs(p.StrandSeparation, s.StrandSeparation) > 1e-12 {
		t.Fatalf("StrandSeparation %g vs %g", p.StrandSeparation, s.StrandSeparation)
	}
}

func TestParallelMatchesSerial(t *testing.T) {
	g := loadOrientedFixture(t)
	m := NewDefaultModel()
	s := ComputeEnsemble(g.Sequence, m, Windows{MinLength: 2, Workers: 1})
	p := ComputeEnsemble(g.Sequence, m, Windows{MinLength: 2, Workers: 8})
	assertEnsembleEquivalent(t, p, s)

	ss := ComputeSignals(g.Sequence, m, Windows{MinLength: 2, Workers: 1}, s, true)
	ps := ComputeSignals(g.Sequence, m, Windows{MinLength: 2, Workers: 8}, s, true)
	assertSignalsEquivalent(t, ps, ss)
}

func TestCircularParallelMatchesSerial(t *testing.T) {
	g := loadOrientedFixture(t)
	m := NewDefaultModel()
	s := ComputeEnsemble(g.Sequence, m, Windows{MinLength: 2, Workers: 1, Circular: true})
	p := ComputeEnsemble(g.Sequence, m, Windows{MinLength: 2, Workers: 8, Circular: true})
	assertEnsembleEquivalent(t, p, s)

	ss := ComputeSignals(g.Sequence, m, Windows{MinLength: 2, Workers: 1, Circular: true}, s, true)
	ps := ComputeSignals(g.Sequence, m, Windows{MinLength: 2, Workers: 8, Circular: true}, s, true)
	assertSignalsEquivalent(t, ps, ss)
}

func TestParallelMoreWorkersThanRows(t *testing.T) {
	m := NewDefaultModel()
	ens := ComputeEnsemble([]byte("GATTACA"), m, Windows{MinLength: 2, Workers: 64})
	if ens.NumStructures != 21 {
		t.Fatalf("count=%d", ens.NumStructures)
	}
}

func TestParallelEmptyChunks(t *testing.T) {
	m := NewDefaultModel()
	seq := []byte("AAAAAAAAAA")
	s := ComputeEnsemble(seq, m, Windows{MinLength: 8, Workers: 1})
	if s.NumStructures != 6 {
		t.Fatalf("setup: NumStructures = %d, want 6", s.NumStructures)
	}
	p := ComputeEnsemble(seq, m, Windows{MinLength: 8, Workers: 8})
	assertEnsembleEquivalent(t, p, s)

	ss := ComputeSignals(seq, m, Windows{MinLength: 8, Workers: 1}, s, true)
	ps := ComputeSignals(seq, m, Windows{MinLength: 8, Workers: 8}, s, true)
	assertSignalsEquivalent(t, ps, ss)
}

func TestCircularParallelEmptyChunks(t *testing.T) {
	m := NewDefaultModel()
	seq := []byte("AAAAAAAAAA")
	s := ComputeEnsemble(seq, m, Windows{MinLength: 8, Workers: 1, Circular: true})
	if s.NumStructures != 20 {
		t.Fatalf("setup: NumStructures = %d, want 20", s.NumStructures)
	}
	p := ComputeEnsemble(seq, m, Windows{MinLength: 8, Workers: 8, Circular: true})
	assertEnsembleEquivalent(t, p, s)

	ss := ComputeSignals(seq, m, Windows{MinLength: 8, Workers: 1, Circular: true}, s, true)
	ps := ComputeSignals(seq, m, Windows{MinLength: 8, Workers: 8, Circular: true}, s, true)
	assertSignalsEquivalent(t, ps, ss)
}

func TestParallelDeterministic(t *testing.T) {
	g := loadOrientedFixture(t)
	m := NewDefaultModel()

	e1 := ComputeEnsemble(g.Sequence, m, Windows{MinLength: 2, Workers: 4})
	s1 := ComputeSignals(g.Sequence, m, Windows{MinLength: 2, Workers: 4}, e1, true)

	e2 := ComputeEnsemble(g.Sequence, m, Windows{MinLength: 2, Workers: 4})
	s2 := ComputeSignals(g.Sequence, m, Windows{MinLength: 2, Workers: 4}, e2, true)

	if e1.LogZ != e2.LogZ {
		t.Fatalf("LogZ not deterministic: %v vs %v", e1.LogZ, e2.LogZ)
	}
	if e1.MFE != e2.MFE {
		t.Fatalf("MFE not deterministic: %+v vs %+v", e1.MFE, e2.MFE)
	}
	if !slices.Equal(s1.BPProb, s2.BPProb) {
		t.Fatalf("BPProb not deterministic")
	}
	if !slices.Equal(s1.AvgG, s2.AvgG) {
		t.Fatalf("AvgG not deterministic")
	}
	if !slices.Equal(s1.MFE, s2.MFE) {
		t.Fatalf("MFE signal not deterministic")
	}
}

func TestCircularParallelDeterministic(t *testing.T) {
	g := loadOrientedFixture(t)
	m := NewDefaultModel()

	e1 := ComputeEnsemble(g.Sequence, m, Windows{MinLength: 2, Workers: 4, Circular: true})
	s1 := ComputeSignals(g.Sequence, m, Windows{MinLength: 2, Workers: 4, Circular: true}, e1, true)

	e2 := ComputeEnsemble(g.Sequence, m, Windows{MinLength: 2, Workers: 4, Circular: true})
	s2 := ComputeSignals(g.Sequence, m, Windows{MinLength: 2, Workers: 4, Circular: true}, e2, true)

	if e1.LogZ != e2.LogZ {
		t.Fatalf("LogZ not deterministic: %v vs %v", e1.LogZ, e2.LogZ)
	}
	if e1.MFE != e2.MFE {
		t.Fatalf("MFE not deterministic: %+v vs %+v", e1.MFE, e2.MFE)
	}
	if !slices.Equal(s1.BPProb, s2.BPProb) {
		t.Fatalf("BPProb not deterministic")
	}
	if !slices.Equal(s1.AvgG, s2.AvgG) {
		t.Fatalf("AvgG not deterministic")
	}
	if !slices.Equal(s1.MFE, s2.MFE) {
		t.Fatalf("MFE signal not deterministic")
	}
}

func TestParallelMFETieBreakMatchesSerial(t *testing.T) {
	m := NewDefaultModel()
	m.Homopolymer, m.OverrideEnergy = true, 0.5
	seq := bytes.Repeat([]byte("A"), 150)

	serial := ComputeEnsemble(seq, m, Windows{MinLength: 2, Workers: 1})
	parallel := ComputeEnsemble(seq, m, Windows{MinLength: 2, Workers: 8})

	assertEnsembleEquivalent(t, parallel, serial)
}

func TestParallelResidualsMatchesSerial(t *testing.T) {
	g := loadOrientedFixture(t)
	m := NewDefaultModel()
	w := Windows{MinLength: 2, Workers: 1}
	ens := ComputeEnsemble(g.Sequence, m, w)
	s := ComputeResiduals(g.Sequence, m, w, ens)
	p := ComputeResiduals(g.Sequence, m, Windows{MinLength: 2, Workers: 8}, ens)
	assertResidualsEquivalent(t, p, s)
}

func TestCircularResidualsMatchesSerial(t *testing.T) {
	g := loadOrientedFixture(t)
	m := NewDefaultModel()
	w := Windows{MinLength: 2, Workers: 1, Circular: true}
	ens := ComputeEnsemble(g.Sequence, m, w)
	s := ComputeResiduals(g.Sequence, m, w, ens)
	p := ComputeResiduals(g.Sequence, m, Windows{MinLength: 2, Workers: 8, Circular: true}, ens)
	assertResidualsEquivalent(t, p, s)

	linW := Windows{MinLength: 2, Workers: 1}
	lin := ComputeResiduals(g.Sequence, m, linW, ComputeEnsemble(g.Sequence, m, linW))
	differs := relOrAbs(s.Twist, lin.Twist) > 1e-9 ||
		relOrAbs(s.LinkingDifference, lin.LinkingDifference) > 1e-9 ||
		relOrAbs(s.Superhelicity, lin.Superhelicity) > 1e-9 ||
		relOrAbs(s.WrappingAbsorption, lin.WrappingAbsorption) > 1e-9 ||
		relOrAbs(s.StrandSeparation, lin.StrandSeparation) > 1e-9
	if !differs {
		t.Fatalf("circular residuals identical to linear residuals; circular bounds may not have engaged: %+v vs %+v", s, lin)
	}
}

func TestParallelResidualsDeterministic(t *testing.T) {
	g := loadOrientedFixture(t)
	m := NewDefaultModel()
	w := Windows{MinLength: 2, Workers: 1}
	ens := ComputeEnsemble(g.Sequence, m, w)

	r1 := ComputeResiduals(g.Sequence, m, Windows{MinLength: 2, Workers: 4}, ens)
	r2 := ComputeResiduals(g.Sequence, m, Windows{MinLength: 2, Workers: 4}, ens)

	if r1 != r2 {
		t.Fatalf("Residuals not deterministic: %+v vs %+v", r1, r2)
	}
}

func TestParallelResidualsOversubscription(t *testing.T) {
	m := NewDefaultModel()
	seq := []byte("GATTACA")
	w := Windows{MinLength: 2, Workers: 1}
	ens := ComputeEnsemble(seq, m, w)
	s := ComputeResiduals(seq, m, w, ens)
	p := ComputeResiduals(seq, m, Windows{MinLength: 2, Workers: 64}, ens)
	assertResidualsEquivalent(t, p, s)
}

func loadPfc53ForBench(b *testing.B) *Gene {
	b.Helper()
	f, err := os.Open("../res/pfc53_full.fa")
	if err != nil {
		b.Fatal(err)
	}
	defer f.Close()
	genes, err := ReadGenes(f)
	if err != nil {
		b.Fatal(err)
	}
	g := genes[0]
	g.ApplyOrientation(false, false, false)
	return g
}

func BenchmarkEnsemblePfc53(b *testing.B) {
	g := loadPfc53ForBench(b)
	m := NewDefaultModel()
	for _, workers := range []int{1, runtime.NumCPU()} {
		b.Run(fmt.Sprintf("workers=%d", workers), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ComputeEnsemble(g.Sequence, m, Windows{MinLength: 2, Workers: workers})
			}
		})
	}
}

func BenchmarkSignalsPfc53(b *testing.B) {
	g := loadPfc53ForBench(b)
	m := NewDefaultModel()
	ens := ComputeEnsemble(g.Sequence, m, Windows{MinLength: 2, Workers: 1})
	for _, workers := range []int{1, runtime.NumCPU()} {
		b.Run(fmt.Sprintf("workers=%d", workers), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				ComputeSignals(g.Sequence, m, Windows{MinLength: 2, Workers: workers}, ens, true)
			}
		})
	}
}
