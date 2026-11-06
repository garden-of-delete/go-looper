package sidd

import (
	"math"
	"testing"
)

// naiveProfile recomputes P, G and melted bases directly from the naive window energies; it is the independent oracle for the profile code.
func naiveProfile(m *Model) (p, g []float64, melted float64) {
	_, all := naiveEnergies(m)
	ground := m.GroundEnergy()
	groundW := math.Exp(-ground/m.RT + Shift)
	zb, zg := groundW, ground*groundW
	sumW := make([]float64, m.EngLen)
	sumEW := make([]float64, m.EngLen)
	for key, E := range all {
		w := math.Exp(-E/m.RT + Shift)
		zb += w
		zg += E * w
		for i := key[0]; i < key[0]+key[1]; i++ {
			sumW[i%m.EngLen] += w
			sumEW[i%m.EngLen] += E * w
		}
	}
	p = make([]float64, m.SeqLen)
	g = make([]float64, m.EngLen)
	maxG := math.Inf(-1)
	for x := range g {
		g[x] = sumEW[x]/sumW[x] - zg/zb
		maxG = math.Max(maxG, g[x])
	}
	for x := range g {
		g[x] -= maxG - GMaxAdj
	}
	for x := range p {
		p[x] = sumW[x] / zb
	}
	for _, w := range sumW {
		melted += w / zb
	}
	return p, g[:m.SeqLen], melted
}

// minMax returns the smallest and largest element of v.
func minMax(v []float64) (lo, hi float64) {
	lo, hi = v[0], v[0]
	for _, x := range v {
		lo, hi = math.Min(lo, x), math.Max(hi, x)
	}
	return lo, hi
}

// TestNaiveProfile pins the single-run P, G and MeltedBases to the naive oracle across geometries.
func TestNaiveProfile(t *testing.T) {
	seq := randomSeq(300, 1)
	cases := []struct {
		name string
		geom Geometry
	}{
		{"linear", GeomLinear},
		{"clamp", GeomClamp},
		{"circular", GeomCircular},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			m := newTestModel(t, seq, func(p *Params) {
				p.Threshold = 1e9
				p.Geometry = tc.geom
			})
			res, err := Run(m)
			if err != nil {
				t.Fatal(err)
			}
			wantP, wantG, wantMelted := naiveProfile(m)
			if d := relDiff(res.MeltedBases, wantMelted); d > 1e-12 {
				t.Errorf("MeltedBases = %v, want %v (rel %g)", res.MeltedBases, wantMelted, d)
			}
			if len(res.P) != m.SeqLen || len(res.G) != m.SeqLen {
				t.Fatalf("len(P) = %d, len(G) = %d, want %d", len(res.P), len(res.G), m.SeqLen)
			}
			for x := range wantP {
				if d := relDiff(res.P[x], wantP[x]); d > 1e-12 {
					t.Fatalf("P[%d] = %v, want %v (rel %g)", x, res.P[x], wantP[x], d)
				}
				if d := math.Abs(res.G[x] - wantG[x]); d > 1e-9 {
					t.Fatalf("G[%d] = %v, want %v", x, res.G[x], wantG[x])
				}
			}
		})
	}
}

// TestRotationInvariance pins that rotating a circular sequence rotates P by the same amount.
func TestRotationInvariance(t *testing.T) {
	const shift = 37
	seq := randomSeq(240, 3)
	rotated := append(append([]byte(nil), seq[shift:]...), seq[:shift]...)
	run := func(s []byte) *Result {
		res, err := Run(newTestModel(t, s, func(p *Params) {
			p.Geometry = GeomCircular
			p.MaxLen = 200
		}))
		if err != nil {
			t.Fatal(err)
		}
		return res
	}
	base, rot := run(seq), run(rotated)
	for i := range rot.P {
		if d := relDiff(rot.P[i], base.P[(i+shift)%len(seq)]); d > 1e-9 {
			t.Fatalf("P[%d] rel diff %g after rotation", i, d)
		}
	}
}

// TestHomopolymerUniform pins a flat P and a maximum G of GMaxAdj on a circular homopolymer.
func TestHomopolymerUniform(t *testing.T) {
	m := newTestModel(t, uniformSeq('A', 300), func(p *Params) { p.Geometry = GeomCircular })
	res, err := Run(m)
	if err != nil {
		t.Fatal(err)
	}
	minP, maxP := minMax(res.P)
	if d := (maxP - minP) / maxP; d > 1e-12 {
		t.Errorf("P spread rel %g", d)
	}
	_, maxG := minMax(res.G)
	if math.Abs(maxG-GMaxAdj) > 1e-9 {
		t.Errorf("max G = %v, want %v", maxG, GMaxAdj)
	}
}

// TestClampDropsPads pins that clamp output covers only the real sequence while the pad bases still affect ZB.
func TestClampDropsPads(t *testing.T) {
	seq := randomSeq(300, 4)
	run := func(g Geometry) *Result {
		res, err := Run(newTestModel(t, seq, func(p *Params) { p.Geometry = g }))
		if err != nil {
			t.Fatal(err)
		}
		return res
	}
	clamp, linear := run(GeomClamp), run(GeomLinear)
	if len(clamp.P) != len(seq) || len(clamp.G) != len(seq) {
		t.Fatalf("len(P) = %d, len(G) = %d, want %d", len(clamp.P), len(clamp.G), len(seq))
	}
	if clamp.ZB == linear.ZB {
		t.Error("clamp ZB equals linear ZB; pad windows contribute nothing")
	}
}

// TestUncoveredFill pins that bases no retained state covers get the maximum G.
func TestUncoveredFill(t *testing.T) {
	seq := append(append(uniformSeq('G', 100), uniformSeq('A', 20)...), uniformSeq('G', 100)...)
	m := newTestModel(t, seq, func(p *Params) {
		p.Geometry = GeomCircular
		p.Sigma = -0.001
		p.MaxLen = 200
	})
	res, err := Run(m)
	if err != nil {
		t.Fatal(err)
	}
	_, maxG := minMax(res.G)
	if math.Abs(maxG-GMaxAdj) > 1e-9 {
		t.Errorf("max G = %v, want %v", maxG, GMaxAdj)
	}
	uncovered, covered := 0, 0
	for x, p := range res.P {
		if p != 0 {
			covered++
			continue
		}
		uncovered++
		if res.G[x] != maxG {
			t.Errorf("uncovered G[%d] = %v, want max %v", x, res.G[x], maxG)
		}
	}
	if uncovered == 0 || covered == 0 {
		t.Fatalf("covered = %d, uncovered = %d; test needs both", covered, uncovered)
	}
}
