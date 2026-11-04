package sidd

import (
	"cmp"
	"fmt"
	"math"
	"slices"
)

type Result struct {
	MinE, MaxE     float64
	OneRunStates   int
	ZB, ZG         float64
	GroundW        float64
	AvgRuns        float64
	TransitionProb float64
	MeltedBases    float64
	P, G           []float64
}

type candidate struct {
	start int
	e     float64
}

func candidates(m *Model) [][]candidate {
	lists := make([][]candidate, m.MaxLen+1)
	for n := m.MinLen; n <= m.MaxLen; n++ {
		starts := m.EngLen
		if m.Geometry == GeomLinear {
			starts = m.SeqLen - n + 1
		}
		if starts < 0 {
			starts = 0
		}
		list := make([]candidate, starts)
		for s := range list {
			list[s] = candidate{start: s, e: Nucleation + m.WindowEnergy(s, n)}
		}
		slices.SortStableFunc(list, func(a, b candidate) int { return cmp.Compare(a.e, b.e) })
		lists[n] = list
	}
	return lists
}

func minimumEnergy(m *Model, lists [][]candidate) float64 {
	minE := m.GroundEnergy()
	for n := m.MinLen; n <= m.MaxLen; n++ {
		if len(lists[n]) == 0 {
			continue
		}
		minE = math.Min(minE, lists[n][0].e+m.GRes(n))
	}
	return minE
}

func forEachRetained(m *Model, lists [][]candidate, minE float64, visit func(s, n int, E, w float64)) {
	maxE := minE + m.Threshold
	for n := m.MinLen; n <= m.MaxLen; n++ {
		gres := m.GRes(n)
		for _, c := range lists[n] {
			E := c.e + gres
			if E >= maxE {
				break
			}
			visit(c.start, n, E, math.Exp(-E/m.RT+Shift))
		}
	}
}

func Run(m *Model) (*Result, error) {
	lists := candidates(m)
	minE := minimumEnergy(m, lists)
	res := &Result{MinE: minE, MaxE: minE + m.Threshold}

	var sumB, sumG float64
	forEachRetained(m, lists, minE, func(s, n int, E, w float64) {
		res.OneRunStates++
		sumB += w
		sumG += float64(E * w)
	})

	res.GroundW = math.Exp(-m.Alpha*m.Alpha*m.K/2/m.RT + Shift)
	res.ZB = sumB + res.GroundW
	res.ZG = sumG + float64(m.Alpha*m.Alpha*m.K/2*res.GroundW)

	if res.ZB == 0 || math.IsNaN(res.ZB) || math.IsInf(res.ZB, 0) {
		return nil, fmt.Errorf("sidd: partition function %v is not representable; adjust the shift or parameters", res.ZB)
	}
	res.AvgRuns = sumB / res.ZB
	res.TransitionProb = (res.ZB - res.GroundW) / res.ZB
	return res, nil
}
