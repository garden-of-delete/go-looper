package sidd

import (
	"fmt"
	"math"
	"slices"
	"testing"
	"time"
)

const (
	anchorTotalCap  = 200
	anchorBracket   = 1e-12
	anchorTolerance = 1e-12
)

type anchorSum struct{ sum, carry float64 }

func (a *anchorSum) add(x float64) {
	t := a.sum + x
	if math.Abs(a.sum) >= math.Abs(x) {
		a.carry += (a.sum - t) + x
	} else {
		a.carry += (x - t) + a.sum
	}
	a.sum = t
}

func (a *anchorSum) value() float64 { return a.sum + a.carry }

func anchorClose(got, want float64) bool {
	return got == want || math.Abs(got-want) <= anchorTolerance*math.Max(math.Abs(got), math.Abs(want))
}

type anchorRun struct{ start, length int }

type naiveMultiResult struct {
	ZB, ZG, AvgRuns, MinE float64
	ZBLow, ZBHigh         float64
	Counts                [5]int
	DefinitelyIn          [5]int
	PossiblyIn            [5]int
}

type anchorSpace struct {
	m        *Model
	maxRuns  int
	circular bool
	windows  [][]float64
}

func newAnchorSpace(m *Model, maxRuns int) *anchorSpace {
	s := &anchorSpace{m: m, maxRuns: maxRuns, circular: m.Geometry != GeomLinear}
	s.windows = make([][]float64, m.MaxLen+1)
	for n := m.MinLen; n <= m.MaxLen; n++ {
		s.windows[n] = make([]float64, s.starts(n))
		for st := range s.windows[n] {
			s.windows[n][st] = m.WindowEnergy(st, n)
		}
	}
	return s
}

func (s *anchorSpace) starts(n int) int {
	if s.circular {
		return s.m.EngLen
	}
	return s.m.SeqLen - n + 1
}

func (s *anchorSpace) apart(a, b anchorRun) bool {
	if !s.circular {
		if a.start > b.start {
			a, b = b, a
		}
		return b.start-(a.start+a.length) >= 2
	}
	size := s.m.EngLen
	if (b.start-a.start+size)%size < a.length || (a.start-b.start+size)%size < b.length {
		return false
	}
	return s.arcOK(a, b) && s.arcOK(b, a)
}

func (s *anchorSpace) arcOK(from, to anchorRun) bool {
	size := s.m.EngLen
	next := from.start + from.length
	gap := ((to.start-next)%size + size) % size
	switch {
	case next-1 >= size-1:
		return gap >= 3
	case next+gap >= size:
		return gap >= 1
	default:
		return gap >= 2
	}
}

func (s *anchorSpace) capsHold(runs []anchorRun) bool {
	total := 0
	for _, r := range runs {
		total += r.length
	}
	switch len(runs) {
	case 1:
		return true
	case 2:
		smaller := min(runs[0].length, runs[1].length)
		return smaller <= s.m.MaxLen/2 && total <= s.m.MaxLen
	default:
		return total <= anchorTotalCap
	}
}

func (s *anchorSpace) each(visit func(runs []anchorRun, E float64)) {
	var extend func(runs []anchorRun, sum float64, total int)
	extend = func(runs []anchorRun, sum float64, total int) {
		first := 0
		if len(runs) > 0 {
			first = runs[len(runs)-1].start + 1
		}
		for n := s.m.MinLen; n <= s.m.MaxLen; n++ {
			for st := first; st < s.starts(n); st++ {
				r := anchorRun{st, n}
				if slices.ContainsFunc(runs, func(p anchorRun) bool { return !s.apart(p, r) }) {
					continue
				}
				next := append(runs, r)
				nextSum := sum + (Nucleation + s.windows[n][st])
				nextTotal := total + n
				if s.capsHold(next) {
					visit(next, nextSum+s.m.GRes(nextTotal))
				}
				if len(next) < s.maxRuns {
					extend(next, nextSum, nextTotal)
				}
			}
		}
	}
	extend(make([]anchorRun, 0, s.maxRuns), 0, 0)
}

func naiveMultiRun(m *Model, maxRuns int, theta float64) naiveMultiResult {
	space := newAnchorSpace(m, maxRuns)
	ground := m.GroundEnergy()
	minE := ground
	space.each(func(_ []anchorRun, E float64) { minE = math.Min(minE, E) })

	edge := minE + theta
	lo := edge - anchorBracket*math.Abs(edge)
	hi := edge + anchorBracket*math.Abs(edge)
	groundW := math.Exp(-ground/m.RT + Shift)
	res := naiveMultiResult{MinE: minE}
	var sumB, sumG, sumRuns, sumLow, sumHigh anchorSum
	space.each(func(runs []anchorRun, E float64) {
		k := len(runs)
		w := math.Exp(-E/m.RT + Shift)
		if E < lo {
			res.DefinitelyIn[k]++
			sumLow.add(w)
		}
		if E < hi {
			res.PossiblyIn[k]++
			sumHigh.add(w)
		}
		if E < edge {
			res.Counts[k]++
			sumB.add(w)
			sumG.add(E * w)
			sumRuns.add(float64(k) * w)
		}
	})
	sumB.add(groundW)
	sumLow.add(groundW)
	sumHigh.add(groundW)
	sumG.add(ground * groundW)
	res.ZB, res.ZBLow, res.ZBHigh, res.ZG = sumB.value(), sumLow.value(), sumHigh.value(), sumG.value()
	res.AvgRuns = sumRuns.value() / res.ZB
	return res
}

type anchorLeg struct {
	name  string
	geom  Geometry
	nn    bool
	sigma float64
	sizes []anchorSize
}

type anchorSize struct {
	runs, seqLen, maxLen int
}

func anchorModel(t *testing.T, seq []byte, leg anchorLeg, size anchorSize, theta float64) *Model {
	t.Helper()
	p := DefaultParams()
	p.Geometry = leg.geom
	p.NearestNeighbor = leg.nn
	p.Sigma = leg.sigma
	p.MinLen = 1
	p.MaxLen = size.maxLen
	p.MaxRuns = size.runs
	p.Threshold = theta
	m, err := NewModel(seq, p)
	if err != nil {
		t.Fatal(err)
	}
	return m
}

func assertAnchorAgreement(t *testing.T, m *Model, theta float64) naiveMultiResult {
	t.Helper()
	want := naiveMultiRun(m, m.MaxRuns, theta)
	got, err := Run(m)
	if err != nil {
		t.Fatal(err)
	}
	if !anchorClose(got.MinE, want.MinE) {
		t.Errorf("MinE = %v, want %v", got.MinE, want.MinE)
	}
	if want.DefinitelyIn == want.PossiblyIn {
		for _, c := range []struct {
			name      string
			got, want float64
		}{
			{"ZB", got.ZB, want.ZB},
			{"ZG", got.ZG, want.ZG},
			{"AvgRuns", got.AvgRuns, want.AvgRuns},
		} {
			if !anchorClose(c.got, c.want) {
				t.Errorf("%s = %v, want %v (rel %.3g)", c.name, c.got, c.want, relDiff(c.got, c.want))
			}
		}
	} else if got.ZB < want.ZBLow*(1-anchorTolerance) || got.ZB > want.ZBHigh*(1+anchorTolerance) {
		t.Errorf("ZB = %v, want within [%v, %v]", got.ZB, want.ZBLow, want.ZBHigh)
	}
	engine := [5]int{0, got.OneRunStates, got.TwoRunStates, got.ThreeRunStates, got.FourRunStates}
	for k := 1; k <= 4; k++ {
		switch {
		case k > m.MaxRuns:
			if engine[k] != 0 {
				t.Errorf("%d-run count = %d with MaxRuns %d, want 0", k, engine[k], m.MaxRuns)
			}
		case theta >= 1e9:
			if engine[k] != want.Counts[k] {
				t.Errorf("%d-run count = %d, want %d", k, engine[k], want.Counts[k])
			}
			if want.Counts[k] == 0 {
				t.Errorf("%d-run class is empty at full retention; test is vacuous", k)
			}
		default:
			if engine[k] < want.DefinitelyIn[k] || engine[k] > want.PossiblyIn[k] {
				t.Errorf("%d-run count = %d, want within [%d, %d] (naive %d)", k, engine[k], want.DefinitelyIn[k], want.PossiblyIn[k], want.Counts[k])
			}
		}
	}
	return want
}

func logAnchor(t *testing.T, want naiveMultiResult, elapsed time.Duration) {
	t.Helper()
	nontrivial := 0
	for k := 1; k <= 4; k++ {
		if want.DefinitelyIn[k] < want.PossiblyIn[k] {
			nontrivial++
		}
	}
	t.Logf("elapsed %v counts %v definitely %v possibly %v nontrivial-brackets %d", elapsed.Round(time.Millisecond), want.Counts[1:], want.DefinitelyIn[1:], want.PossiblyIn[1:], nontrivial)
}

func partialMultiRun(want naiveMultiResult, theta float64) bool {
	if theta >= 1e9 {
		return false
	}
	for k := 2; k <= 4; k++ {
		if want.Counts[k] > 0 {
			return true
		}
	}
	return false
}

func TestNaiveMultiRunAgreement(t *testing.T) {
	partial := 0
	clampSizes := []anchorSize{{2, 90, 12}, {3, 36, 5}, {4, 10, 2}}
	openSizes := []anchorSize{{2, 90, 12}, {3, 50, 6}, {4, 30, 4}}
	legs := []anchorLeg{
		{"clamp-copolymeric", GeomClamp, false, -0.06, clampSizes},
		{"circular-copolymeric", GeomCircular, false, -0.2, openSizes},
		{"clamp-nearest-neighbor", GeomClamp, true, -0.2, clampSizes},
		{"circular-nearest-neighbor", GeomCircular, true, -0.06, openSizes},
		{"linear-copolymeric", GeomLinear, false, -0.2, openSizes},
	}
	for _, leg := range legs {
		for _, size := range leg.sizes {
			for _, theta := range []float64{1e9, 6, 25, 55} {
				t.Run(fmt.Sprintf("%s/runs%d/theta%g", leg.name, size.runs, theta), func(t *testing.T) {
					m := anchorModel(t, randomSeq(size.seqLen, 20260411), leg, size, theta)
					start := time.Now()
					want := assertAnchorAgreement(t, m, theta)
					logAnchor(t, want, time.Since(start))
					if partialMultiRun(want, theta) {
						partial++
					}
				})
			}
		}
	}
	if partial == 0 {
		t.Error("no finite-theta subtest retained a multi-run state; partial-retention coverage is vacuous")
	}
}

func TestNaiveShortCircular(t *testing.T) {
	partial := 0
	leg := anchorLeg{name: "circular", geom: GeomCircular, sigma: -0.06}
	sizes := []anchorSize{{2, 60, 20}, {3, 60, 8}, {4, 40, 4}}
	for _, size := range sizes {
		for _, theta := range []float64{1e9, 3, 25, 55} {
			t.Run(fmt.Sprintf("runs%d/theta%g", size.runs, theta), func(t *testing.T) {
				m := anchorModel(t, uniformSeq('A', size.seqLen), leg, size, theta)
				start := time.Now()
				want := assertAnchorAgreement(t, m, theta)
				logAnchor(t, want, time.Since(start))
				if partialMultiRun(want, theta) {
					partial++
				}
			})
		}
	}
	if partial == 0 {
		t.Error("no finite-theta subtest retained a multi-run state; partial-retention coverage is vacuous")
	}
}

func TestNaiveTiedWindowEdge(t *testing.T) {
	leg := anchorLeg{name: "circular", geom: GeomCircular, sigma: -0.06}
	size := anchorSize{2, 60, 20}
	m := anchorModel(t, uniformSeq('A', size.seqLen), leg, size, 0)
	tied := 2*(Nucleation+m.WindowEnergy(0, 5)) + m.GRes(10)
	theta := tied - m.GroundEnergy()
	m = anchorModel(t, uniformSeq('A', size.seqLen), leg, size, theta)
	want := assertAnchorAgreement(t, m, theta)
	if want.DefinitelyIn[2] >= want.PossiblyIn[2] {
		t.Fatalf("two-run bracket [%d, %d] is trivial; the tied edge was not hit", want.DefinitelyIn[2], want.PossiblyIn[2])
	}
	logAnchor(t, want, 0)
}

func TestNaiveMultiRunProfile(t *testing.T) {
	clampSizes := []anchorSize{{3, 36, 5}, {4, 10, 2}}
	openSizes := []anchorSize{{3, 50, 6}, {4, 30, 4}}
	legs := []anchorLeg{
		{"clamp-copolymeric", GeomClamp, false, -0.06, clampSizes},
		{"circular-copolymeric", GeomCircular, false, -0.2, openSizes},
		{"linear-copolymeric", GeomLinear, false, -0.2, openSizes},
		{"clamp-nearest-neighbor", GeomClamp, true, -0.2, clampSizes},
		{"circular-nearest-neighbor", GeomCircular, true, -0.06, openSizes},
	}
	for _, leg := range legs {
		for _, size := range leg.sizes {
			classStates := 0
			for _, theta := range []float64{25, 40, 1e9} {
				t.Run(fmt.Sprintf("%s/runs%d/theta%g", leg.name, size.runs, theta), func(t *testing.T) {
					classStates += assertProfileAnchor(t, leg, size, theta)
				})
			}
			if classStates == 0 {
				t.Errorf("%s/runs%d: no %d-run state retained at any theta; the class's profile is unanchored", leg.name, size.runs, size.runs)
			}
		}
	}
}

func assertProfileAnchor(t *testing.T, leg anchorLeg, size anchorSize, theta float64) int {
	t.Helper()
	m := anchorModel(t, randomSeq(size.seqLen, 20260411), leg, size, theta)
	res, err := Run(m)
	if err != nil {
		t.Fatal(err)
	}

	space := newAnchorSpace(m, m.MaxRuns)
	ground := m.GroundEnergy()
	minE := ground
	space.each(func(_ []anchorRun, E float64) { minE = math.Min(minE, E) })
	edge := minE + m.Threshold
	groundW := math.Exp(-ground/m.RT + Shift)
	var zbSum, zgSum anchorSum
	multi, classStates := 0, 0
	sumB := make([]anchorSum, m.EngLen)
	sumG := make([]anchorSum, m.EngLen)
	space.each(func(runs []anchorRun, E float64) {
		if E >= edge {
			return
		}
		w := math.Exp(-E/m.RT + Shift)
		if len(runs) >= 2 {
			multi++
		}
		if len(runs) == size.runs {
			classStates++
		}
		zbSum.add(w)
		zgSum.add(E * w)
		for _, r := range runs {
			for i := r.start; i < r.start+r.length; i++ {
				sumB[i%m.EngLen].add(w)
				sumG[i%m.EngLen].add(E * w)
			}
		}
	})

	if multi == 0 {
		t.Fatal("no retained state with two or more runs; profile test is vacuous")
	}
	zbSum.add(groundW)
	zgSum.add(ground * groundW)
	zb, zg := zbSum.value(), zgSum.value()
	gx := make([]float64, m.EngLen)
	maxG := math.Inf(-1)
	for x := range gx {
		if sumB[x].value() > 0 {
			gx[x] = sumG[x].value()/sumB[x].value() - zg/zb
			maxG = math.Max(maxG, gx[x])
		}
	}
	if math.IsInf(maxG, 0) || math.IsNaN(maxG) {
		t.Fatalf("maxG = %v; profile test is vacuous", maxG)
	}
	for x := range gx {
		if sumB[x].value() == 0 {
			gx[x] = maxG
		}
		gx[x] += GMaxAdj - maxG
	}
	t.Logf("retained multi-run states: %d", multi)
	if len(res.P) != m.SeqLen || len(res.G) != m.SeqLen {
		t.Fatalf("profile lengths P %d, G %d, want %d", len(res.P), len(res.G), m.SeqLen)
	}
	var melted anchorSum
	for x := range sumB {
		melted.add(sumB[x].value() / zb)
	}
	if !anchorClose(res.MeltedBases, melted.value()) {
		t.Errorf("MeltedBases = %v, want %v", res.MeltedBases, melted.value())
	}
	for x := 0; x < m.SeqLen; x++ {
		if !anchorClose(res.P[x], sumB[x].value()/zb) {
			t.Errorf("P[%d] = %v, want %v", x, res.P[x], sumB[x].value()/zb)
		}
		if !(math.Abs(res.G[x]-gx[x]) <= anchorTolerance*math.Max(1, math.Abs(gx[x]))) {
			t.Errorf("G[%d] = %v, want %v", x, res.G[x], gx[x])
		}
	}
	return classStates
}
