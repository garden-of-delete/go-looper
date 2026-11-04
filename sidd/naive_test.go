package sidd

import "math"

type naiveResult struct {
	zb, zg, minE float64
	states       int
}

func naiveWindows(m *Model) [][2]int {
	var windows [][2]int
	for n := m.MinLen; n <= m.MaxLen; n++ {
		last := m.EngLen - 1
		if m.Geometry == GeomLinear {
			last = m.SeqLen - n
		}
		for s := 0; s <= last; s++ {
			windows = append(windows, [2]int{s, n})
		}
	}
	return windows
}

func naiveEnergies(m *Model) (float64, map[[2]int]float64) {
	minE := m.GroundEnergy()
	all := map[[2]int]float64{}
	for _, w := range naiveWindows(m) {
		E := Nucleation + m.WindowEnergy(w[0], w[1]) + m.GRes(w[1])
		all[w] = E
		minE = math.Min(minE, E)
	}
	return minE, all
}

func naiveRun(m *Model) naiveResult {
	minE, all := naiveEnergies(m)
	limit := minE + m.Threshold
	ground := m.GroundEnergy()
	groundW := math.Exp(-ground/m.RT + Shift)
	res := naiveResult{minE: minE, zb: groundW, zg: ground * groundW}
	for _, E := range all {
		if E >= limit {
			continue
		}
		w := math.Exp(-E/m.RT + Shift)
		res.zb += w
		res.zg += E * w
		res.states++
	}
	return res
}
