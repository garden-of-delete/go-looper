package rlooper

import (
	"math"
	"math/rand"
	"testing"
)

func TestDefaults(t *testing.T) {
	m := NewDefaultModel()
	near(t, m.K, 0.9029123033, 1e-9)
	near(t, m.Alpha, -10.0961538, 1e-6)
	near(t, m.RT(), 0.6156220250, 1e-9)
	near(t, m.GroundStateEnergy(), 36.017974, 1e-4)
}

func TestSetterCoupling(t *testing.T) {
	m := NewDefaultModel()
	m.SetN(200)
	near(t, m.K, 2200*GasConstant*310/200, 1e-12)
	near(t, m.Alpha, 200*-0.07/10.4, 1e-12)
	m.SetSigma(-0.10)
	near(t, m.Alpha, 200*-0.10/10.4, 1e-12)
}

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

func TestPrefixMatchesNaive(t *testing.T) {
	m := NewDefaultModel()
	r := rand.New(rand.NewSource(42))
	seq := make([]byte, 100)
	for i := range seq {
		seq[i] = "ACGT"[r.Intn(4)]
	}
	p := PrefixBPEnergies(seq, m)
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

func near(t *testing.T, got, want, tol float64) {
	t.Helper()
	if math.Abs(got-want) > tol {
		t.Fatalf("got %g want %g", got, want)
	}
}
