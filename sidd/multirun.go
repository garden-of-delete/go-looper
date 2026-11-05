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

func classCap(maxLen, runs int) int {
	if runs == 2 {
		return maxLen
	}
	return multiRunLenCap
}

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

func (s *classScan) run() {
	if s.runs*s.m.MinLen > s.total || math.IsInf(s.minRun, 1) {
		return
	}
	s.lengths(0, s.m.MinLen, 0, 0)
}

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

func (s *classScan) clear(level, start, length int) bool {
	for q := 0; q < level; q++ {
		if !separated(s.spans[q].start, s.spans[q].length, start, length, s.m.EngLen, s.circular) {
			return false
		}
	}
	return true
}

func globalMinimum(m *Model, lists [][]candidate) float64 {
	minE := minimumEnergy(m, lists)
	for runs := 2; runs <= m.MaxRuns; runs++ {
		newClassScan(m, lists, runs, &minE, func(_ []runSpan, E float64) {
			minE = E
		}).run()
	}
	return minE
}

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

func (c *compensatedSum) add(x float64) {
	t := c.sum + x
	if math.Abs(c.sum) >= math.Abs(x) {
		c.carry += (c.sum - t) + x
	} else {
		c.carry += (x - t) + c.sum
	}
	c.sum = t
}

func (c *compensatedSum) value() float64 {
	if !isFinite(c.carry) {
		return c.sum
	}
	return c.sum + c.carry
}

func runMulti(m *Model) (*Result, error) {
	lists := candidates(m)
	minE := globalMinimum(m, lists)
	if err := checkRepresentable(m, minE); err != nil {
		return nil, err
	}
	res := &Result{MinE: minE, MaxE: minE + m.Threshold}

	var sumB, sumG, sumRuns compensatedSum
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
	return res, nil
}
