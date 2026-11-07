package sidd

import (
	"fmt"
	"math"
)

const (
	maxRunClasses  = 4
	multiRunLenCap = 200
)

type runSpan struct{ start, length int }

// separated reports whether two runs are far enough apart, mirroring SIST's overlap(), including its deliberately asymmetric origin-crossing rule (at least 2 intact base pairs off-origin, at least 3 when a run reaches the last index, else at least 1).
func separated(p1, n1, p2, n2 int, engLen int, circular bool) bool {
	if !circular {
		if p1+n1 > engLen || p2+n2 > engLen {
			return false
		}
		return p1+n1 < p2-1 || p2+n2 < p1-1
	}
	reach1 := p1+n1 > engLen-1
	reach2 := p2+n2 > engLen-1
	switch {
	case !reach1 && !reach2:
		return p1+n1 < p2-1 || p2+n2 < p1-1
	case reach2 && !reach1:
		return p1+n1 < p2-1 && n2-(engLen-1-p2) < p1-1
	case reach1 && !reach2:
		return p2+n2 < p1-1 && n1-(engLen-1-p1) < p2-1
	default:
		return false
	}
}

// classCap returns the cap on total melted length for a runs-run class: MaxLen for two runs, 200 for three and four.
func classCap(maxLen, runs int) int {
	if runs == 2 {
		return maxLen
	}
	return multiRunLenCap
}

// lengthBound returns the largest length the next run may take so that the remaining runs, each at least as long, still fit within total.
func lengthBound(total, sum, remaining int) int {
	return (total - sum) / remaining
}

type classScan struct {
	m        *Model
	lists    [][]candidate
	runs     int
	total    int
	minRun   float64
	minGres  float64
	limit    *float64
	visit    func(runs []runSpan, E float64)
	circular bool
	lens     [maxRunClasses]int
	spans    [maxRunClasses]runSpan
	picked   [maxRunClasses]int
	heads    [maxRunClasses]float64
	gres     float64
}

// newClassScan prepares a branch-and-bound scan of the runs-run class, caching the minimum single-run energy and minimum GRes used by its pruning bounds.
func newClassScan(m *Model, lists [][]candidate, runs int, limit *float64, visit func(runs []runSpan, E float64)) *classScan {
	s := &classScan{
		m:        m,
		lists:    lists,
		runs:     runs,
		total:    classCap(m.MaxLen, runs),
		minRun:   math.Inf(1),
		minGres:  math.Inf(1),
		limit:    limit,
		visit:    visit,
		circular: m.Geometry != GeomLinear,
	}
	for n := m.MinLen; n <= m.MaxLen; n++ {
		if len(lists[n]) > 0 {
			s.minRun = math.Min(s.minRun, lists[n][0].e)
		}
	}
	for t := runs * m.MinLen; t <= min(s.total, runs*m.MaxLen); t++ {
		s.minGres = math.Min(s.minGres, m.GRes(t))
	}
	return s
}

// run executes the scan, doing nothing when the class cannot fit under its length cap.
func (s *classScan) run() {
	if s.runs*s.m.MinLen > s.total || math.IsInf(s.minRun, 1) {
		return
	}
	s.lengths(0, s.m.MinLen, 0, 0)
}

// lengths enumerates non-decreasing run-length tuples, pruning any prefix whose energy lower bound reaches *limit.
func (s *classScan) lengths(level, prev, sum int, headSum float64) {
	last := level == s.runs-1
	hi := min(s.m.MaxLen, lengthBound(s.total, sum, s.runs-level))
	for n := prev; n <= hi; n++ {
		list := s.lists[n]
		if len(list) == 0 {
			continue
		}
		hs := headSum + list[0].e
		var lb float64
		if last {
			lb = hs + s.m.GRes(sum+n)
		} else {
			lb = hs
			for k := level + 1; k < s.runs; k++ {
				lb += s.minRun
			}
			lb += s.minGres
		}
		if lb >= *s.limit {
			continue
		}
		s.lens[level] = n
		s.heads[level] = list[0].e
		if last {
			s.gres = s.m.GRes(sum + n)
			s.positions(0, 0)
		} else {
			s.lengths(level+1, n, sum+n, hs)
		}
	}
}

// positions enumerates start positions for the runs of the current length tuple, starting past the previous pick when lengths are equal so each unordered set of runs is visited once, where SIST visits such a set under both orderings. The runs slice handed to visit is a reused window into s and is overwritten on the next state; callers must not retain it.
func (s *classScan) positions(level int, eSum float64) {
	n := s.lens[level]
	list := s.lists[n]
	first := 0
	if level > 0 && s.lens[level-1] == n {
		first = s.picked[level-1] + 1
	}
	last := level == s.runs-1
	for i := first; i < len(list); i++ {
		c := list[i]
		part := eSum + c.e
		bound := part
		for k := level + 1; k < s.runs; k++ {
			bound += s.heads[k]
		}
		bound += s.gres
		if bound >= *s.limit {
			break
		}
		if !s.clear(level, c.start, n) {
			continue
		}
		s.spans[level] = runSpan{c.start, n}
		s.picked[level] = i
		if last {
			s.visit(s.spans[:s.runs], bound)
		} else {
			s.positions(level+1, part)
		}
	}
}

// clear reports whether the run at level is separated from every run already placed.
func (s *classScan) clear(level, start, length int) bool {
	for q := 0; q < level; q++ {
		if !separated(s.spans[q].start, s.spans[q].length, start, length, s.m.EngLen, s.circular) {
			return false
		}
	}
	return true
}

// globalMinimum returns the true minimum energy over all states of up to MaxRuns runs, found before any retention window is applied (phase one of the two-phase theta window, fixing SIST's moving window).
func globalMinimum(m *Model, lists [][]candidate) float64 {
	minE := minimumEnergy(m, lists)
	for runs := 2; runs <= m.MaxRuns; runs++ {
		newClassScan(m, lists, runs, &minE, func(_ []runSpan, E float64) {
			minE = E
		}).run()
	}
	return minE
}

// forEachRetainedState visits every state within Threshold of the frozen global minimum minE (phase two of the theta window), single-run states first and then each multi-run class. The runs slice passed to visit is a reused window into the scan state and is overwritten on the next state; a visitor must copy it to retain it.
func forEachRetainedState(m *Model, lists [][]candidate, minE float64, visit func(runs []runSpan, E, w float64)) {
	var single [1]runSpan
	forEachRetained(m, lists, minE, func(start, n int, E, w float64) {
		single[0] = runSpan{start, n}
		visit(single[:], E, w)
	})
	maxE := minE + m.Threshold
	for runs := 2; runs <= m.MaxRuns; runs++ {
		newClassScan(m, lists, runs, &maxE, func(spans []runSpan, E float64) {
			visit(spans, E, math.Exp(-E/m.RT+Shift))
		}).run()
	}
}

type compensatedSum struct{ sum, carry float64 }

// add adds x using Neumaier compensated summation.
func (c *compensatedSum) add(x float64) {
	t := c.sum + x
	if math.Abs(c.sum) >= math.Abs(x) {
		c.carry += (c.sum - t) + x
	} else {
		c.carry += (x - t) + c.sum
	}
	c.sum = t
}

// value returns the compensated total, falling back to the plain sum if the carry is not finite.
func (c *compensatedSum) value() float64 {
	if !isFinite(c.carry) {
		return c.sum
	}
	return c.sum + c.carry
}

type compensatedCell struct{ sumG, sumB compensatedSum }

type compensatedMatrix [][]compensatedCell

// newCompensatedMatrix allocates the per-length, per-start accumulator grid for the multi-run profile.
func newCompensatedMatrix(m *Model) compensatedMatrix {
	cells := make(compensatedMatrix, m.MaxLen+1)
	for n := m.MinLen; n <= m.MaxLen; n++ {
		cells[n] = make([]compensatedCell, m.EngLen)
	}
	return cells
}

// add credits a state's weighted energy and weight to the cell of each of its runs.
func (c compensatedMatrix) add(runs []runSpan, E, w float64) {
	ew := float64(E * w)
	for _, r := range runs {
		cl := &c[r.length][r.start]
		cl.sumG.add(ew)
		cl.sumB.add(w)
	}
}

// collapse reduces the compensated cells to plain cells so the shared profile accumulation can consume them.
func (c compensatedMatrix) collapse() cellMatrix {
	out := make(cellMatrix, len(c))
	for n, row := range c {
		if row == nil {
			continue
		}
		out[n] = make([]cell, len(row))
		for s := range row {
			out[n][s] = cell{sumG: row[s].sumG.value(), sumB: row[s].sumB.value()}
		}
	}
	return out
}

// runMulti computes the ensemble for MaxRuns above 1, using Neumaier-compensated sums because plain summation lost up to 6e-11 relative; runSingle keeps plain sums for byte parity.
func runMulti(m *Model) (*Result, error) {
	lists := candidates(m)
	minE := globalMinimum(m, lists)
	if err := checkRepresentable(m, minE); err != nil {
		return nil, err
	}
	res := &Result{MinE: minE, MaxE: minE + m.Threshold}

	var sumB, sumG, sumRuns compensatedSum
	cells := newCompensatedMatrix(m)
	forEachRetainedState(m, lists, minE, func(runs []runSpan, E, w float64) {
		switch len(runs) {
		case 1:
			res.OneRunStates++
		case 2:
			res.TwoRunStates++
		case 3:
			res.ThreeRunStates++
		default:
			res.FourRunStates++
		}
		sumB.add(w)
		sumG.add(float64(E * w))
		sumRuns.add(float64(float64(len(runs)) * w))
		cells.add(runs, E, w)
	})

	res.GroundW = math.Exp(-m.Alpha*m.Alpha*m.K/2/m.RT + Shift)
	retainedB := sumB.value()
	sumB.add(res.GroundW)
	sumG.add(float64(m.Alpha * m.Alpha * m.K / 2 * res.GroundW))
	res.ZB = sumB.value()
	res.ZG = sumG.value()

	if res.ZB == 0 || !isFinite(res.ZB) || !isFinite(res.ZG) {
		return nil, fmt.Errorf("sidd: partition function (ZB %v, ZG %v) is not representable at these parameters; shorten the sequence or use a less negative sigma", res.ZB, res.ZG)
	}
	res.TransitionProb = retainedB / res.ZB
	res.AvgRuns = sumRuns.value() / res.ZB
	cells.collapse().accumulate(m).finish(m, res.ZB, res.ZG, res)
	return res, nil
}
