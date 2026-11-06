package rlooper

import (
	"math"
	"math/rand"
	"testing"
)

// TestDefaults pins the derived K, Alpha, RT and ground-state energy of the default model.
func TestDefaults(t *testing.T) {
	m := NewDefaultModel()
	near(t, m.K, 0.9029123033, 1e-9)
	near(t, m.Alpha, -10.0961538, 1e-6)
	near(t, m.RT(), 0.6156220250, 1e-9)
	near(t, m.GroundStateEnergy(), 36.017974, 1e-4)
}

// TestSetterCoupling pins that SetN and SetSigma recompute the derived K and Alpha.
func TestSetterCoupling(t *testing.T) {
	m := NewDefaultModel()
	m.SetN(200)
	near(t, m.K, 6.771842275, 1e-9)
	near(t, m.Alpha, -1.346153846153846, 1e-12)
	m.SetSigma(-0.10)
	near(t, m.Alpha, -1.923076923076923, 1e-12)
}

// TestSupercoilingEnergy pins the zero minimum at m=105, reference values at 2 and 500, and that Unconstrained drops the term.
func TestSupercoilingEnergy(t *testing.T) {
	m := NewDefaultModel()
	if v := m.SupercoilingEnergy(105); math.Abs(v) > 1e-9 {
		t.Fatalf("minimum at m=105 not zero: %g", v)
	}
	near(t, m.SupercoilingEnergy(2), 43.18, 0.01)
	near(t, m.SupercoilingEnergy(500), 88.57, 0.01)
	m.Unconstrained = true
	if m.SupercoilingEnergy(50) != 0 {
		t.Fatal("unconstrained must drop the term")
	}
}

// TestStepEnergyGattaca pins the GATTACA step-energy sum and the asymmetry of the CG and GC steps.
func TestStepEnergyGattaca(t *testing.T) {
	m := NewDefaultModel()
	seq := []byte("GATTACA")
	var sum float64
	for i := 0; i+1 < len(seq); i++ {
		sum += m.StepEnergy(seq[i], seq[i+1])
	}
	near(t, sum, 1.53, 1e-12)
	near(t, m.StepEnergy('C', 'G'), -0.16, 1e-12)
	near(t, m.StepEnergy('G', 'C'), 0.97, 1e-12)
}

// TestPrefixMatchesNaive pins prefixBP differences against naive per-window step sums on a random sequence.
func TestPrefixMatchesNaive(t *testing.T) {
	m := NewDefaultModel()
	r := rand.New(rand.NewSource(42))
	seq := make([]byte, 100)
	for i := range seq {
		seq[i] = "ACGT"[r.Intn(4)]
	}
	p := prefixBP(seq, m, false)
	for s := 0; s < len(seq); s++ {
		for e := s + 1; e < len(seq); e++ {
			var naive float64
			for i := s; i < e; i++ {
				naive += m.StepEnergy(seq[i], seq[i+1])
			}
			if math.Abs((p[e]-p[s])-naive) > 1e-9 {
				t.Fatalf("s=%d e=%d", s, e)
			}
		}
	}
}

// near fails unless got is within the absolute tolerance tol of want.
func near(t *testing.T, got, want, tol float64) {
	t.Helper()
	if math.Abs(got-want) > tol {
		t.Fatalf("got %g want %g", got, want)
	}
}
