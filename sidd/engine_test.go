package sidd

import (
	"math/rand"
	"strings"
	"testing"
)

// randomSeq returns n pseudo-random bases from a fixed seed so tests are reproducible.
func randomSeq(n int, seed int64) []byte {
	r := rand.New(rand.NewSource(seed))
	seq := make([]byte, n)
	for i := range seq {
		seq[i] = "ACGT"[r.Intn(4)]
	}
	return seq
}

// uniformSeq returns n copies of base.
func uniformSeq(base byte, n int) []byte {
	return []byte(strings.Repeat(string(base), n))
}

// newTestModel builds a MaxRuns-1 model from default parameters after applying mutate, failing the test on error.
func newTestModel(t *testing.T, seq []byte, mutate func(*Params)) *Model {
	t.Helper()
	p := DefaultParams()
	p.MaxRuns = 1
	if mutate != nil {
		mutate(&p)
	}
	m, err := NewModel(seq, p)
	if err != nil {
		t.Fatal(err)
	}
	return m
}

// TestNaiveAgreement pins the single-run engine to the naive oracle (states, ZB, ZG, MinE within 1e-12) across all geometries with full retention.
func TestNaiveAgreement(t *testing.T) {
	seq := randomSeq(300, 1)
	cases := []struct {
		name   string
		mutate func(*Params)
	}{
		{"linear", func(p *Params) { p.Geometry = GeomLinear }},
		{"clamp", func(p *Params) { p.Geometry = GeomClamp }},
		{"circular", func(p *Params) { p.Geometry = GeomCircular }},
		{"clamp nearest-neighbor", func(p *Params) { p.Geometry = GeomClamp; p.NearestNeighbor = true }},
		{"circular nearest-neighbor", func(p *Params) { p.Geometry = GeomCircular; p.NearestNeighbor = true }},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			m := newTestModel(t, seq, func(p *Params) {
				p.Threshold = 1e9
				tc.mutate(p)
			})
			res, err := Run(m)
			if err != nil {
				t.Fatal(err)
			}
			want := naiveRun(m)
			if res.OneRunStates != want.states {
				t.Errorf("states = %d, want %d", res.OneRunStates, want.states)
			}
			if d := relDiff(res.ZB, want.zb); d > 1e-12 {
				t.Errorf("ZB = %v, want %v (rel %g)", res.ZB, want.zb, d)
			}
			if d := relDiff(res.ZG, want.zg); d > 1e-12 {
				t.Errorf("ZG = %v, want %v (rel %g)", res.ZG, want.zg, d)
			}
			if d := relDiff(res.MinE, want.minE); d > 1e-12 {
				t.Errorf("MinE = %v, want %v", res.MinE, want.minE)
			}
		})
	}
}

// TestThetaRetention pins forEachRetained to exactly the windows with energy below minE+Threshold.
func TestThetaRetention(t *testing.T) {
	m := newTestModel(t, uniformSeq('A', 60), func(p *Params) {
		p.Geometry = GeomCircular
		p.MaxLen = 30
		p.Threshold = 12
	})
	res, err := Run(m)
	if err != nil {
		t.Fatal(err)
	}
	minE, _ := naiveEnergies(m)
	limit := minE + 12

	retained := map[[2]int]bool{}
	forEachRetained(m, candidates(m), res.MinE, func(s, n int, E, w float64) {
		retained[[2]int{s, n}] = true
		if E >= limit {
			t.Errorf("retained (%d,%d) E=%v >= %v", s, n, E, limit)
		}
	})

	_, all := naiveEnergies(m)
	for key, E := range all {
		if E < limit && !retained[key] {
			t.Errorf("missing (%d,%d) E=%v < %v", key[0], key[1], E, limit)
		}
		if E >= limit && retained[key] {
			t.Errorf("unexpected (%d,%d) E=%v", key[0], key[1], E)
		}
	}
	if res.OneRunStates != len(retained) {
		t.Errorf("OneRunStates = %d, retained %d", res.OneRunStates, len(retained))
	}
}

// TestScanOrderTies pins that tied candidate energies keep ascending start order, matching SIST's scan order.
func TestScanOrderTies(t *testing.T) {
	m := newTestModel(t, uniformSeq('A', 60), func(p *Params) {
		p.Geometry = GeomCircular
		p.MaxLen = 10
		p.Threshold = 1e9
	})
	lists := candidates(m)
	for n := m.MinLen; n <= m.MaxLen; n++ {
		if len(lists[n]) != m.EngLen {
			t.Fatalf("length %d: %d candidates, want %d", n, len(lists[n]), m.EngLen)
		}
		for i, c := range lists[n] {
			if c.start != i {
				t.Fatalf("length %d: position %d holds start %d", n, i, c.start)
			}
		}
	}

	next := map[int]int{}
	forEachRetained(m, lists, 0, func(s, n int, E, w float64) {
		if s != next[n] {
			t.Fatalf("length %d: visited start %d, want %d", n, s, next[n])
		}
		next[n]++
	})
}

// TestLinearNeverWraps pins that linear geometry generates no window that runs past the sequence end.
func TestLinearNeverWraps(t *testing.T) {
	m := newTestModel(t, randomSeq(40, 2), func(p *Params) {
		p.Geometry = GeomLinear
		p.MaxLen = 20
	})
	lists := candidates(m)
	for n := m.MinLen; n <= m.MaxLen; n++ {
		if want := m.SeqLen - n + 1; len(lists[n]) != want {
			t.Errorf("length %d: %d candidates, want %d", n, len(lists[n]), want)
		}
		for _, c := range lists[n] {
			if c.start+n > m.SeqLen {
				t.Errorf("window (%d,%d) wraps", c.start, n)
			}
		}
	}
}

// TestOverflowGuard pins the above-melting-point overflow error.
func TestOverflowGuard(t *testing.T) {
	m := newTestModel(t, uniformSeq('A', 2000), func(p *Params) {
		p.Geometry = GeomClamp
		p.Temperature = 400
	})
	_, err := Run(m)
	if err == nil || !strings.Contains(err.Error(), "melting point") {
		t.Fatalf("error = %v, want overflow error", err)
	}
}

// TestUnderflowGuard pins the underflow error when every Boltzmann weight underflows.
func TestUnderflowGuard(t *testing.T) {
	m := newTestModel(t, uniformSeq('G', 20000), func(p *Params) {
		p.Geometry = GeomCircular
		p.Sigma = -0.2
	})
	_, err := Run(m)
	if err == nil || !strings.Contains(err.Error(), "underflows") {
		t.Fatalf("error = %v, want underflow error", err)
	}
}

// TestPartialUnderflowGuard pins the underflow error when the best weight is only subnormal rather than flushed to zero, so the near-minimum weights are the only ones with value left, at reduced precision.
func TestPartialUnderflowGuard(t *testing.T) {
	m := newTestModel(t, randomSeq(10000, 3), func(p *Params) {
		p.Geometry = GeomClamp
		p.Sigma = -0.119
	})
	_, err := Run(m)
	if err == nil || !strings.Contains(err.Error(), "underflows") {
		t.Fatalf("error = %v, want underflow error", err)
	}
}

// TestSigmaZero pins that sigma 0 yields MinE 0 and a transition probability strictly between 0 and 1.
func TestSigmaZero(t *testing.T) {
	m := newTestModel(t, uniformSeq('A', 600), func(p *Params) {
		p.Geometry = GeomLinear
		p.Sigma = 0
	})
	res, err := Run(m)
	if err != nil {
		t.Fatal(err)
	}
	if res.TransitionProb <= 0 || res.TransitionProb >= 1 {
		t.Errorf("TransitionProb = %v", res.TransitionProb)
	}
	if res.MinE != 0 {
		t.Errorf("MinE = %v, want 0", res.MinE)
	}
}
