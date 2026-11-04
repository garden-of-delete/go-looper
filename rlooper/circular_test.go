package rlooper

import (
	"bytes"
	"fmt"
	"math"
	"testing"
)

func TestCircularEnsembleMatchesIndependentEnumeration(t *testing.T) {
	seq := []byte("ACGTAC")
	n := len(seq)
	for _, minLength := range []int{2, 3, 5} {
		t.Run(fmt.Sprintf("minLength=%d", minLength), func(t *testing.T) {
			m := NewDefaultModel()
			w := Windows{MinLength: minLength, Workers: 1, Circular: true}
			ens := ComputeEnsemble(seq, m, w)
			if ens.NumStructures != n*(n-minLength) {
				t.Fatalf("count=%d want %d", ens.NumStructures, n*(n-minLength))
			}

			var z float64
			naive := make(map[[2]int]float64)
			for s := 0; s < n; s++ {
				for mm := minLength; mm <= n-1; mm++ {
					e := m.SupercoilingEnergy(mm)
					for j := 0; j < mm-1; j++ {
						e += m.StepEnergy(seq[(s+j)%n], seq[(s+j+1)%n])
					}
					naive[[2]int{s, (s + mm - 1) % n}] = e
					z += math.Exp(-e / m.RT())
				}
			}
			z += math.Exp(-m.GroundStateEnergy() / m.RT())
			logZ := math.Log(z)
			if rel(ens.LogZ, logZ) > 1e-12 {
				t.Fatalf("logZ %g vs %g", ens.LogZ, logZ)
			}

			fourPiSqC := 4 * math.Pi * math.Pi * m.C
			var naiveRes Residuals
			for s := 0; s < n; s++ {
				for mm := minLength; mm <= n-1; mm++ {
					e := naive[[2]int{s, (s + mm - 1) % n}]
					prob := math.Exp(-e/m.RT() - logZ)
					l := float64(mm - 1)
					denom := fourPiSqC + m.K*l
					rld := (fourPiSqC / denom) * (m.Alpha + l*m.A)
					rtw := (2 * math.Pi * m.K / denom) * (m.Alpha + l*m.A)
					naiveRes.Twist += rtw * prob
					naiveRes.LinkingDifference += rld * prob
					naiveRes.WrappingAbsorption += (m.Alpha - rld + l*m.A) * prob
					naiveRes.StrandSeparation -= l * m.A * prob
				}
			}
			pGround := math.Exp(-m.GroundStateEnergy()/m.RT() - logZ)
			naiveRes.Twist += pGround * m.K * m.C * m.Alpha / (2 * math.Pi)
			naiveRes.LinkingDifference += pGround * m.Alpha
			naiveRes.Superhelicity = naiveRes.LinkingDifference / (m.N * m.A)

			recs := AllStructures(seq, m, w, ens)
			if len(recs) != len(naive) {
				t.Fatalf("recs %d vs naive %d", len(recs), len(naive))
			}
			var sum float64
			for _, r := range recs {
				want, ok := naive[[2]int{r.Start, r.End}]
				if !ok {
					t.Fatalf("unexpected window %d-%d", r.Start, r.End)
				}
				if rel(r.Energy, want) > 1e-12 {
					t.Fatalf("window %d-%d energy %g want %g", r.Start, r.End, r.Energy, want)
				}
				sum += r.Prob
			}
			if math.Abs(sum+ens.PGround-1) > 1e-9 {
				t.Fatalf("probability sum %g", sum+ens.PGround)
			}

			res := ComputeResiduals(seq, m, w, ens)
			if rel(res.Twist, naiveRes.Twist) > 1e-12 ||
				rel(res.LinkingDifference, naiveRes.LinkingDifference) > 1e-12 ||
				rel(res.Superhelicity, naiveRes.Superhelicity) > 1e-12 ||
				rel(res.WrappingAbsorption, naiveRes.WrappingAbsorption) > 1e-12 ||
				rel(res.StrandSeparation, naiveRes.StrandSeparation) > 1e-12 {
				t.Fatalf("residuals %+v vs naive %+v", res, naiveRes)
			}
		})
	}
}

func TestCircularLinearSubsetBitParity(t *testing.T) {
	m := NewDefaultModel()
	g := loadOrientedFixture(t)
	lin := AllStructures(g.Sequence, m, Windows{MinLength: 2, Workers: 1}, ComputeEnsemble(g.Sequence, m, Windows{MinLength: 2, Workers: 1}))
	circ := AllStructures(g.Sequence, m, Windows{MinLength: 2, Workers: 1, Circular: true}, ComputeEnsemble(g.Sequence, m, Windows{MinLength: 2, Workers: 1, Circular: true}))
	linByKey := make(map[[2]int]float64, len(lin))
	for _, r := range lin {
		linByKey[[2]int{r.Start, r.End}] = r.Energy
	}
	checked := 0
	for _, r := range circ {
		if r.End < r.Start {
			continue
		}
		if want, ok := linByKey[[2]int{r.Start, r.End}]; ok {
			if r.Energy != want {
				t.Fatalf("window %d-%d: circular %v != linear %v", r.Start, r.End, r.Energy, want)
			}
			checked++
		}
	}
	if checked != len(lin)-1 {
		t.Fatalf("checked %d, want %d", checked, len(lin)-1)
	}
}

func TestCircularAnalyticHomopolymerZ(t *testing.T) {
	m := NewDefaultModel()
	m.Homopolymer, m.OverrideEnergy = true, 0.5
	n := 150
	seq := bytes.Repeat([]byte("A"), n)
	w := Windows{MinLength: 2, Workers: 1, Circular: true}
	ens := ComputeEnsemble(seq, m, w)
	var z float64
	for mm := 2; mm <= n-1; mm++ {
		e := m.SupercoilingEnergy(mm) + float64(mm-1)*0.5
		z += float64(n) * math.Exp(-e/m.RT())
	}
	z += math.Exp(-m.GroundStateEnergy() / m.RT())
	if rel(ens.LogZ, math.Log(z)) > 1e-12 {
		t.Fatalf("logZ %g vs analytic %g", ens.LogZ, math.Log(z))
	}
}

func TestCircularNoWindows(t *testing.T) {
	m := NewDefaultModel()
	ens := ComputeEnsemble([]byte("AC"), m, Windows{MinLength: 2, Workers: 1, Circular: true})
	if ens.NumStructures != 0 || math.Abs(ens.PGround-1) > 1e-12 || math.IsNaN(ens.LogZ) {
		t.Fatalf("%+v", ens)
	}
	ens = ComputeEnsemble([]byte("ACGTAC"), m, Windows{MinLength: 6, Workers: 1, Circular: true})
	if ens.NumStructures != 0 || math.Abs(ens.PGround-1) > 1e-12 {
		t.Fatalf("%+v", ens)
	}
}

func TestCircularTopStructuresMatchesAllStructures(t *testing.T) {
	m := NewDefaultModel()
	seq := []byte("ACGTAC")
	w := Windows{MinLength: 2, Workers: 1, Circular: true}
	ens := ComputeEnsemble(seq, m, w)
	all := AllStructures(seq, m, w, ens)
	const k = 3
	top := TopStructures(seq, m, w, ens, k)
	if len(top) != k {
		t.Fatalf("len(top) = %d, want %d", len(top), k)
	}
	wrapped := false
	for i := range top {
		if top[i] != all[i] {
			t.Fatalf("TopStructures[%d] = %+v, want AllStructures[%d] = %+v", i, top[i], i, all[i])
		}
		if top[i].End < top[i].Start {
			wrapped = true
		}
	}
	if !wrapped {
		t.Fatalf("top-%d set contains no wrapped record (End < Start); choose a different k/n", k)
	}
}

func TestCircularSignalsMatchIndependentEnumeration(t *testing.T) {
	seq := []byte("ACGTACGT")
	n := len(seq)
	for _, minLength := range []int{2, 3, 5} {
		t.Run(fmt.Sprintf("minLength=%d", minLength), func(t *testing.T) {
			m := NewDefaultModel()
			w := Windows{MinLength: minLength, Workers: 1, Circular: true}
			ens := ComputeEnsemble(seq, m, w)
			sig := ComputeSignals(seq, m, w, ens, true)
			bp := make([]float64, n)
			num := make([]float64, n)
			den := make([]float64, n)
			for s := 0; s < n; s++ {
				for mm := minLength; mm <= n-1; mm++ {
					g := m.SupercoilingEnergy(mm)
					for j := 0; j < mm-1; j++ {
						g += m.StepEnergy(seq[(s+j)%n], seq[(s+j+1)%n])
					}
					prob := math.Exp(-g/m.RT() - ens.LogZ)
					for j := 0; j <= mm-2; j++ {
						bp[(s+j)%n] += prob
						num[(s+j)%n] += prob * g
						den[(s+j)%n] += prob
					}
				}
			}
			for i := 0; i < n; i++ {
				if rel(sig.BPProb[i], bp[i]) > 1e-12 {
					t.Fatalf("bpprob[%d] %g vs %g", i, sig.BPProb[i], bp[i])
				}
				wantAvg := 0.0
				if den[i] > 0 {
					wantAvg = num[i] / den[i]
				}
				if rel(sig.AvgG[i], wantAvg) > 1e-12 {
					t.Fatalf("avgG[%d] %g vs %g", i, sig.AvgG[i], wantAvg)
				}
			}
		})
	}
}

func TestCircularRotationInvariance(t *testing.T) {
	m := NewDefaultModel()
	g := loadOrientedFixture(t)
	n := len(g.Sequence)
	k := 37
	rot := make([]byte, n)
	for i := range rot {
		rot[i] = g.Sequence[(i+k)%n]
	}
	w := Windows{MinLength: 2, Workers: 1, Circular: true}
	e1 := ComputeEnsemble(g.Sequence, m, w)
	e2 := ComputeEnsemble(rot, m, w)
	if rel(e1.LogZ, e2.LogZ) > 1e-12 || rel(e1.PGround, e2.PGround) > 1e-12 {
		t.Fatalf("Z not rotation invariant: %g vs %g", e1.LogZ, e2.LogZ)
	}
	r1 := AllStructures(g.Sequence, m, w, e1)
	if r1[0].Energy == r1[1].Energy {
		t.Fatal("fixture has tied minimum; rotation MFE mapping requires a unique one")
	}
	if e2.MFE.Start != (e1.MFE.Start-k+n)%n || e2.MFE.End != (e1.MFE.End-k+n)%n || rel(e2.MFE.Energy, e1.MFE.Energy) > 1e-12 {
		t.Fatalf("MFE did not rotate: orig %+v rot %+v", e1.MFE, e2.MFE)
	}
	s1 := ComputeSignals(g.Sequence, m, w, e1, true)
	s2 := ComputeSignals(rot, m, w, e2, true)
	for i := 0; i < n; i++ {
		j := (i + k) % n
		if rel(s1.BPProb[j], s2.BPProb[i]) > 1e-12 {
			t.Fatalf("bpprob not rotated at %d", i)
		}
		if rel(s1.AvgG[j], s2.AvgG[i]) > 1e-12 {
			t.Fatalf("avgG not rotated at %d", i)
		}
		if s1.MFE[j] != s2.MFE[i] {
			t.Fatalf("MFE track not rotated at %d: orig[%d]=%g rot[%d]=%g", i, j, s1.MFE[j], i, s2.MFE[i])
		}
	}
}

func TestCircularSignalsHomopolymerUniform(t *testing.T) {
	m := NewDefaultModel()
	m.Homopolymer, m.OverrideEnergy = true, 0.5
	seq := bytes.Repeat([]byte("A"), 120)
	w := Windows{MinLength: 2, Workers: 1, Circular: true}
	ens := ComputeEnsemble(seq, m, w)
	sig := ComputeSignals(seq, m, w, ens, false)
	for i, v := range sig.BPProb {
		if rel(v, sig.BPProb[0]) > 1e-12 {
			t.Fatalf("bpprob[%d]=%g != bpprob[0]=%g", i, v, sig.BPProb[0])
		}
	}
}

func TestCircularSignalsWrappedMFE(t *testing.T) {
	m := NewDefaultModel()
	g := loadOrientedFixture(t)
	n := len(g.Sequence)
	w := Windows{MinLength: 2, Workers: 1, Circular: true}
	base := ComputeEnsemble(g.Sequence, m, w)
	k := (base.MFE.Start + base.MFE.Length(n)/2) % n
	seq := make([]byte, n)
	for i := range seq {
		seq[i] = g.Sequence[(i+k)%n]
	}
	ens := ComputeEnsemble(seq, m, w)
	if ens.MFE.End >= ens.MFE.Start {
		t.Fatalf("rotation by %d was supposed to wrap the MFE, got %+v", k, ens.MFE)
	}
	sig := ComputeSignals(seq, m, w, ens, false)
	count := 0
	for j := 0; j <= ens.MFE.Length(n)-2; j++ {
		if sig.MFE[(ens.MFE.Start+j)%n] != 1.0 {
			t.Fatalf("MFE track missing credit at offset %d", j)
		}
		count++
	}
	for i, v := range sig.MFE {
		if v != 0 && v != 1 {
			t.Fatalf("MFE[%d]=%g", i, v)
		}
	}
	total := 0
	for _, v := range sig.MFE {
		if v == 1 {
			total++
		}
	}
	if total != count {
		t.Fatalf("MFE track has %d ones, expected %d", total, count)
	}
}
