package sidd

const uncoveredG = -10000.0

type cell struct {
	sumG, sumB float64
}

type cellMatrix [][]cell

// newCellMatrix allocates the per-length, per-start accumulator grid for the single-run profile.
func newCellMatrix(m *Model) cellMatrix {
	cells := make(cellMatrix, m.MaxLen+1)
	for n := m.MinLen; n <= m.MaxLen; n++ {
		cells[n] = make([]cell, m.EngLen)
	}
	return cells
}

// add credits a state's weighted energy and weight to its cell using plain summation.
func (c cellMatrix) add(s, n int, E, w float64) {
	c[n][s].sumG = c[n][s].sumG + float64(E*w)
	c[n][s].sumB = c[n][s].sumB + w
}

type profile struct {
	sumG, sumB []float64
}

// accumulate spreads each cell's sums across the bases its window covers, wrapping modulo EngLen.
func (c cellMatrix) accumulate(m *Model) profile {
	p := profile{sumG: make([]float64, m.EngLen), sumB: make([]float64, m.EngLen)}
	for n := m.MinLen; n <= m.MaxLen; n++ {
		for s, cl := range c[n] {
			if cl.sumB == 0 && cl.sumG == 0 {
				continue
			}
			for i := s; i < s+n; i++ {
				x := i % m.EngLen
				p.sumG[x] = p.sumG[x] + cl.sumG
				p.sumB[x] = p.sumB[x] + cl.sumB
			}
		}
	}
	return p
}

// finish converts the per-base sums into P and G on res, filling uncovered bases with the maximum G and shifting G so its maximum is GMaxAdj.
func (p profile) finish(m *Model, zb, zg float64, res *Result) {
	aveG := zg / zb
	px := make([]float64, m.EngLen)
	gx := make([]float64, m.EngLen)
	for x := range gx {
		px[x] = p.sumB[x] / zb
		if p.sumB[x] == 0 {
			gx[x] = uncoveredG
		} else {
			gx[x] = p.sumG[x]/p.sumB[x] - aveG
		}
	}

	maxG := gx[0]
	for _, g := range gx[1:] {
		if maxG < g {
			maxG = g
		}
	}
	for x, g := range gx {
		if g == uncoveredG {
			gx[x] = maxG
		}
	}

	for _, v := range px {
		res.MeltedBases += v
	}
	adj := maxG - GMaxAdj
	for x := range gx {
		gx[x] = gx[x] - adj
	}
	res.P = px[:m.SeqLen]
	res.G = gx[:m.SeqLen]
}
