package sidd

import (
	"math"
	"strings"
	"testing"
)

func covers(start, length, engLen, cell int) bool {
	return (cell-start+engLen)%engLen < length
}

func refSeparated(p1, n1, p2, n2, engLen int) bool {
	for c := 0; c < engLen; c++ {
		if covers(p1, n1, engLen, c) && covers(p2, n2, engLen, c) {
			return false
		}
	}
	arc := func(from, length, to int) bool {
		start := (from + length) % engLen
		gap := (to - start + engLen) % engLen
		crossesOrigin := gap >= 1 && start+gap >= engLen
		switch {
		case crossesOrigin:
			return gap >= 1
		case from+length >= engLen:
			return gap >= 3
		default:
			return gap >= 2
		}
	}
	return arc(p1, n1, p2) && arc(p2, n2, p1)
}

func TestSeparated(t *testing.T) {
	const L = 20
	cases := []struct {
		name           string
		p1, n1, p2, n2 int
		circular       bool
		want           bool
	}{
		{"off-origin gap 1", 2, 3, 6, 3, true, false},
		{"off-origin gap 2", 2, 3, 7, 3, true, true},
		{"off-origin gap 1 swapped", 6, 3, 2, 3, true, false},
		{"off-origin gap 2 swapped", 7, 3, 2, 3, true, true},
		{"origin-crossing run ends at last index gap 2", 17, 3, 2, 3, true, false},
		{"origin-crossing run ends at last index gap 3", 17, 3, 3, 3, true, true},
		{"origin-crossing wrapped run gap 2", 18, 4, 4, 3, true, false},
		{"origin-crossing wrapped run gap 3", 18, 4, 5, 3, true, true},
		{"origin-crossing swapped last-index gap 2", 2, 3, 17, 3, true, false},
		{"origin-crossing swapped last-index gap 3", 3, 3, 17, 3, true, true},
		{"origin-crossing neither reaches last index gap 1", 14, 5, 0, 3, true, true},
		{"origin-crossing neither reaches last index gap 2", 14, 4, 0, 3, true, true},
		{"both reach last index", 15, 5, 17, 4, true, false},
		{"adjacent runs", 2, 3, 5, 3, true, false},
		{"overlapping runs", 2, 5, 4, 3, true, false},
		{"linear gap 1", 2, 3, 6, 3, false, false},
		{"linear gap 2", 2, 3, 7, 3, false, true},
		{"linear run ending at last index gap 2", 15, 5, 2, 3, false, true},
		{"linear wraparound pair", 18, 4, 4, 3, false, false},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := separated(tc.p1, tc.n1, tc.p2, tc.n2, L, tc.circular); got != tc.want {
				t.Errorf("separated = %v, want %v", got, tc.want)
			}
			if got := separated(tc.p2, tc.n2, tc.p1, tc.n1, L, tc.circular); got != tc.want {
				t.Errorf("separated (swapped arguments) = %v, want %v", got, tc.want)
			}
		})
	}
}

func TestSeparatedMatchesReference(t *testing.T) {
	const L = 14
	for p1 := 0; p1 < L; p1++ {
		for n1 := 1; n1 <= 5; n1++ {
			for p2 := 0; p2 < L; p2++ {
				for n2 := 1; n2 <= 5; n2++ {
					if got, want := separated(p1, n1, p2, n2, L, true), refSeparated(p1, n1, p2, n2, L); got != want {
						t.Fatalf("separated(%d,%d,%d,%d) = %v, reference %v", p1, n1, p2, n2, got, want)
					}
				}
			}
		}
	}
}

func naiveMultiCounts(m *Model, minE float64) (counts [5]int, sumB, sumRuns float64) {
	type win struct {
		start, length int
		e             float64
	}
	var wins []win
	for _, w := range naiveWindows(m) {
		wins = append(wins, win{w[0], w[1], Nucleation + m.WindowEnergy(w[0], w[1])})
	}
	limit := minE + m.Threshold
	capFor := func(r int) int {
		if r == 2 {
			return m.MaxLen
		}
		return 200
	}
	circular := m.Geometry != GeomLinear
	ok := func(chosen []int, cand int) bool {
		for _, q := range chosen {
			a, b := wins[q], wins[cand]
			if circular {
				if !refSeparated(a.start, a.length, b.start, b.length, m.EngLen) {
					return false
				}
			} else if !(a.start+a.length < b.start-1 || b.start+b.length < a.start-1) {
				return false
			}
		}
		return true
	}
	var rec func(chosen []int, first, r int, totalLen int, sumE float64)
	rec = func(chosen []int, first, r int, totalLen int, sumE float64) {
		if len(chosen) == r {
			E := sumE + m.GRes(totalLen)
			if E < limit {
				counts[r]++
				w := math.Exp(-E/m.RT + Shift)
				sumB += w
				sumRuns += float64(r) * w
			}
			return
		}
		for i := first; i < len(wins); i++ {
			w := wins[i]
			if totalLen+w.length > capFor(r) {
				continue
			}
			if r == 2 && len(chosen) == 0 && w.length > m.MaxLen/2 {
				continue
			}
			if !ok(chosen, i) {
				continue
			}
			rec(append(chosen, i), i+1, r, totalLen+w.length, sumE+w.e)
		}
	}
	for r := 2; r <= m.MaxRuns; r++ {
		rec(nil, 0, r, 0, 0)
	}
	return counts, sumB, sumRuns
}

func TestPairDedupTinyCase(t *testing.T) {
	m := newTestModel(t, uniformSeq('A', 12), func(p *Params) {
		p.Geometry = GeomCircular
		p.MaxLen = 3
		p.MaxRuns = 2
		p.Threshold = 1e9
	})
	res, err := Run(m)
	if err != nil {
		t.Fatal(err)
	}
	want := 0
	for n1 := 1; n1 <= m.MaxLen/2; n1++ {
		for n2 := n1; n1+n2 <= m.MaxLen; n2++ {
			for p1 := 0; p1 < 12; p1++ {
				for p2 := 0; p2 < 12; p2++ {
					if n1 == n2 && p2 <= p1 {
						continue
					}
					if refSeparated(p1, n1, p2, n2, 12) {
						want++
					}
				}
			}
		}
	}
	if want == 0 {
		t.Fatal("closed-form count is zero; test setup is vacuous")
	}
	if res.TwoRunStates != want {
		t.Errorf("TwoRunStates = %d, want %d", res.TwoRunStates, want)
	}
}

func TestMultiRunMatchesNaive(t *testing.T) {
	cases := []struct {
		name      string
		seq       []byte
		geom      Geometry
		maxLen    int
		threshold float64
		fourRuns  bool
	}{
		{"circular", randomSeq(15, 7), GeomCircular, 4, 70, true},
		{"linear", randomSeq(16, 8), GeomLinear, 4, 70, true},
		{"clamp", randomSeq(8, 9), GeomClamp, 4, 30, false},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			m := newTestModel(t, tc.seq, func(p *Params) {
				p.Geometry = tc.geom
				p.MaxLen = tc.maxLen
				p.MaxRuns = 4
				p.Threshold = tc.threshold
			})
			res, err := Run(m)
			if err != nil {
				t.Fatal(err)
			}
			counts, sumB, sumRuns := naiveMultiCounts(m, res.MinE)
			got := [5]int{0, 0, res.TwoRunStates, res.ThreeRunStates, res.FourRunStates}
			for r := 2; r <= 4; r++ {
				if got[r] != counts[r] {
					t.Errorf("%d-run states = %d, want %d", r, got[r], counts[r])
				}
			}
			if tc.fourRuns && counts[4] == 0 {
				t.Fatal("no four-run states; test setup is vacuous")
			}
			oneB := res.TransitionProb*res.ZB - sumB
			wantAvg := (oneB + sumRuns) / res.ZB
			if relDiff(res.AvgRuns, wantAvg) > 1e-9 {
				t.Errorf("AvgRuns = %v, want %v", res.AvgRuns, wantAvg)
			}
		})
	}
}

func TestLengthBounds(t *testing.T) {
	const maxLen = 250
	for r := 2; r <= 4; r++ {
		total := classCap(maxLen, r)
		var got, want int
		var walk func(level, prev, sum int)
		walk = func(level, prev, sum int) {
			if level == r {
				got++
				return
			}
			for n := max(prev, 1); n <= min(maxLen, lengthBound(total, sum, r-level)); n++ {
				walk(level+1, n, sum+n)
			}
		}
		walk(0, 1, 0)
		var naive func(level, prev, sum int)
		naive = func(level, prev, sum int) {
			if level == r {
				want++
				return
			}
			for n := prev; n <= maxLen; n++ {
				if sum+n*(r-level) > total {
					break
				}
				naive(level+1, n, sum+n)
			}
		}
		naive(0, 1, 0)
		if got != want {
			t.Errorf("%d-run length tuples = %d, want %d", r, got, want)
		}
		if r == 2 && lengthBound(total, 0, 2) != maxLen/2 {
			t.Errorf("2R first-length bound = %d, want %d", lengthBound(total, 0, 2), maxLen/2)
		}
		if r > 2 && classCap(maxLen, r) != 200 {
			t.Errorf("%d-run cap = %d, want 200", r, classCap(maxLen, r))
		}
	}
}

func islandSeq(island, gap, count int) []byte {
	unit := strings.Repeat("A", island) + strings.Repeat("G", gap)
	return []byte(strings.Repeat(unit, count))
}

func TestCapsOnVisitedStates(t *testing.T) {
	cases := []struct {
		name                         string
		island, gap, count           int
		maxLen                       int
		sigma, theta                 float64
		wantClasses                  [5]bool
		wantMaxTotal2, wantMaxTotal3 int
	}{
		{name: "all classes visited", island: 3, gap: 12, count: 5, maxLen: 6, sigma: -0.4, theta: 14,
			wantClasses: [5]bool{false, true, true, true, true}},
		{name: "200 and MaxLen caps bind", island: 80, gap: 30, count: 3, maxLen: 250, sigma: -1.0, theta: 1,
			wantClasses: [5]bool{false, false, true, true, false}, wantMaxTotal2: 250, wantMaxTotal3: 200},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			m := newTestModel(t, islandSeq(tc.island, tc.gap, tc.count), func(p *Params) {
				p.Geometry = GeomCircular
				p.MaxLen = tc.maxLen
				p.MaxRuns = 4
				p.Sigma = tc.sigma
				p.Threshold = tc.theta
			})
			lists := candidates(m)
			minE := globalMinimum(m, lists)
			var seen, maxTotal [5]int
			forEachRetainedState(m, lists, minE, func(runs []runSpan, E, w float64) {
				total := 0
				for _, r := range runs {
					total += r.length
				}
				seen[len(runs)]++
				maxTotal[len(runs)] = max(maxTotal[len(runs)], total)
				switch {
				case len(runs) == 2 && (min(runs[0].length, runs[1].length) > m.MaxLen/2 || total > m.MaxLen):
					t.Fatalf("2R cap violated: %v", runs)
				case len(runs) > 2 && total > 200:
					t.Fatalf("%dR cap violated: %v", len(runs), runs)
				}
			})
			t.Logf("visited by class (1R..4R): %v, max totals: %v", seen[1:], maxTotal[1:])
			for r := 1; r <= 4; r++ {
				if tc.wantClasses[r] && seen[r] == 0 {
					t.Errorf("no %dR states visited; test setup is vacuous", r)
				}
			}
			if tc.wantMaxTotal2 != 0 && maxTotal[2] != tc.wantMaxTotal2 {
				t.Errorf("max 2R total = %d, want the MaxLen cap %d to bind", maxTotal[2], tc.wantMaxTotal2)
			}
			if tc.wantMaxTotal3 != 0 && maxTotal[3] != tc.wantMaxTotal3 {
				t.Errorf("max 3R total = %d, want the 200 cap to bind", maxTotal[3])
			}
		})
	}
}

func TestTwoPhaseWindow(t *testing.T) {
	island := strings.Repeat("A", 15)
	gc := strings.Repeat("G", 40)
	seq := []byte(island + gc + island + gc)
	m := newTestModel(t, seq, func(p *Params) {
		p.Geometry = GeomCircular
		p.MaxLen = 30
		p.MaxRuns = 2
		p.Sigma = -0.4
		p.Threshold = 3
	})

	oneMin, all := naiveEnergies(m)
	twoMin := math.Inf(1)
	windows := naiveWindows(m)
	for _, a := range windows {
		for _, b := range windows {
			if a[1] > b[1] || a[1] > m.MaxLen/2 || a[1]+b[1] > m.MaxLen {
				continue
			}
			E := Nucleation + m.WindowEnergy(a[0], a[1]) + Nucleation + m.WindowEnergy(b[0], b[1]) + m.GRes(a[1]+b[1])
			if E < twoMin && refSeparated(a[0], a[1], b[0], b[1], m.EngLen) {
				twoMin = E
			}
		}
	}
	if !(twoMin < oneMin) {
		t.Fatalf("setup: two-run minimum %v is not below the one-run minimum %v", twoMin, oneMin)
	}

	res, err := Run(m)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(res.MinE-twoMin) > 1e-12 {
		t.Errorf("MinE = %v, want two-run minimum %v", res.MinE, twoMin)
	}
	want := 0
	for _, E := range all {
		if E < twoMin+m.Threshold {
			want++
		}
	}
	if want == 0 {
		t.Fatal("setup: no one-run states within the window of the two-run minimum")
	}
	if res.OneRunStates != want {
		t.Errorf("OneRunStates = %d, want %d", res.OneRunStates, want)
	}
}

func TestMultiRunGuards(t *testing.T) {
	m := newTestModel(t, uniformSeq('A', 2000), func(p *Params) {
		p.Geometry = GeomClamp
		p.Temperature = 400
		p.MaxRuns = 4
	})
	if _, err := Run(m); err == nil || !strings.Contains(err.Error(), "melting point") {
		t.Fatalf("error = %v, want overflow error", err)
	}

	m = newTestModel(t, uniformSeq('G', 20000), func(p *Params) {
		p.Geometry = GeomCircular
		p.Sigma = -0.2
		p.MaxRuns = 4
	})
	if _, err := Run(m); err == nil || !strings.Contains(err.Error(), "underflows") {
		t.Fatalf("error = %v, want underflow error", err)
	}
}

func TestMaxRunsValidation(t *testing.T) {
	for _, runs := range []int{0, 5, -1} {
		p := DefaultParams()
		p.MaxRuns = runs
		_, err := NewModel(uniformSeq('A', 500), p)
		if err == nil || !strings.Contains(err.Error(), "maximum runs") {
			t.Errorf("MaxRuns %d: error = %v, want maximum runs rejection", runs, err)
		}
	}
	for _, runs := range []int{1, 4} {
		p := DefaultParams()
		p.MaxRuns = runs
		if _, err := NewModel(uniformSeq('A', 500), p); err != nil {
			t.Errorf("MaxRuns %d: unexpected error %v", runs, err)
		}
	}
}

func TestAvgRuns1RMode(t *testing.T) {
	m := newTestModel(t, randomSeq(300, 5), func(p *Params) {
		p.Geometry = GeomCircular
		p.MaxRuns = 1
	})
	res, err := Run(m)
	if err != nil {
		t.Fatal(err)
	}
	if res.TwoRunStates != 0 || res.ThreeRunStates != 0 || res.FourRunStates != 0 {
		t.Errorf("multi-run counts = %d/%d/%d, want zeros", res.TwoRunStates, res.ThreeRunStates, res.FourRunStates)
	}
	if res.AvgRuns != res.TransitionProb {
		t.Errorf("AvgRuns = %v, TransitionProb = %v", res.AvgRuns, res.TransitionProb)
	}
}
