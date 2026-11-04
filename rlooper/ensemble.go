package rlooper

import (
	"cmp"
	"math"
	"slices"
	"sync"
)

type EnsembleResult struct {
	LogZ, LogGround, PGround float64
	MFE                      StructureRecord
	NumStructures            int
}

type StructureRecord struct {
	// TODO: add type to support multi-structure competition
	Start, End int
	Energy     float64
	Prob       float64
}

type logSumExp struct {
	max, sum float64
}

func newLogSumExp() logSumExp { return logSumExp{max: math.Inf(-1)} }

func (l *logSumExp) add(x float64) {
	if x <= l.max {
		l.sum += math.Exp(x - l.max)
		return
	}
	l.sum = l.sum*math.Exp(l.max-x) + 1
	l.max = x
}

func (l *logSumExp) merge(o logSumExp) {
	if math.IsInf(o.max, -1) {
		return
	}
	if math.IsInf(l.max, -1) {
		*l = o
		return
	}
	if o.max <= l.max {
		l.sum += o.sum * math.Exp(o.max-l.max)
		return
	}
	l.sum = l.sum*math.Exp(l.max-o.max) + o.sum
	l.max = o.max
}

func (l *logSumExp) value() float64 {
	if math.IsInf(l.max, -1) {
		return math.Inf(-1)
	}
	return l.max + math.Log(l.sum)
}

func compareStructures(a, b StructureRecord) int {
	return cmp.Or(
		cmp.Compare(a.Energy, b.Energy),
		cmp.Compare(a.Start, b.Start),
		cmp.Compare(a.End, b.End),
	)
}

func windowEnergy(m *EnergyModel, p []float64, s, mLen int) float64 {
	return m.SupercoilingEnergy(mLen) + (p[s+mLen-1] - p[s])
}

type Windows struct {
	MinLength int
	Circular  bool
	Workers   int
}

func maxLen(n, s int, circular bool) int {
	if circular {
		return n - 1
	}
	return n - s
}

func windowEnd(s, mLen, n int, circular bool) int {
	e := s + mLen - 1
	if circular {
		e %= n
	}
	return e
}

func clampWorkers(workers, rows int) int {
	if workers > rows {
		return rows
	}
	return workers
}

type ensembleAccum struct {
	lse   logSumExp
	mfe   StructureRecord
	count int
}

func newEnsembleAccum() ensembleAccum {
	return ensembleAccum{lse: newLogSumExp()}
}

func (a *ensembleAccum) addWindow(s, e int, g, rt float64) {
	a.lse.add(-g / rt)
	rec := StructureRecord{Start: s, End: e, Energy: g}
	if a.count == 0 || compareStructures(rec, a.mfe) < 0 {
		a.mfe = rec
	}
	a.count++
}

func (a *ensembleAccum) merge(o ensembleAccum) {
	a.lse.merge(o.lse)
	if o.count == 0 {
		return
	}
	if a.count == 0 || compareStructures(o.mfe, a.mfe) < 0 {
		a.mfe = o.mfe
	}
	a.count += o.count
}

func ComputeEnsemble(seq []byte, m *EnergyModel, w Windows) EnsembleResult {
	p := prefixBP(seq, m, w.Circular)
	rt := m.RT()
	n := len(seq)
	workers := clampWorkers(w.Workers, n)

	acc := newEnsembleAccum()
	if workers <= 1 {
		for s := 0; s < n; s++ {
			for mLen := w.MinLength; mLen <= maxLen(n, s, w.Circular); mLen++ {
				acc.addWindow(s, windowEnd(s, mLen, n, w.Circular), windowEnergy(m, p, s, mLen), rt)
			}
		}
	} else {
		results := make([]ensembleAccum, workers)
		var wg sync.WaitGroup
		for c := 0; c < workers; c++ {
			wg.Add(1)
			go func(c int) {
				defer wg.Done()
				chunk := newEnsembleAccum()
				for s := c; s < n; s += workers {
					for mLen := w.MinLength; mLen <= maxLen(n, s, w.Circular); mLen++ {
						chunk.addWindow(s, windowEnd(s, mLen, n, w.Circular), windowEnergy(m, p, s, mLen), rt)
					}
				}
				results[c] = chunk
			}(c)
		}
		wg.Wait()
		for _, r := range results {
			acc.merge(r)
		}
	}

	logGround := -m.GroundStateEnergy() / rt
	acc.lse.add(logGround)
	logZ := acc.lse.value()

	if acc.count > 0 {
		acc.mfe.Prob = math.Exp(-acc.mfe.Energy/rt - logZ)
	}

	return EnsembleResult{
		LogZ:          logZ,
		LogGround:     logGround,
		PGround:       math.Exp(logGround - logZ),
		MFE:           acc.mfe,
		NumStructures: acc.count,
	}
}

func AllStructures(seq []byte, m *EnergyModel, w Windows, ens EnsembleResult) []StructureRecord {
	p := prefixBP(seq, m, w.Circular)
	rt := m.RT()
	n := len(seq)
	recs := make([]StructureRecord, 0, ens.NumStructures)

	for s := 0; s < n; s++ {
		for mLen := w.MinLength; mLen <= maxLen(n, s, w.Circular); mLen++ {
			g := windowEnergy(m, p, s, mLen)
			recs = append(recs, StructureRecord{
				Start:  s,
				End:    windowEnd(s, mLen, n, w.Circular),
				Energy: g,
				Prob:   math.Exp(-g/rt - ens.LogZ),
			})
		}
	}

	slices.SortFunc(recs, compareStructures)

	return recs
}

func TopStructures(seq []byte, m *EnergyModel, w Windows, ens EnsembleResult, n int) []StructureRecord {
	if n <= 0 {
		return nil
	}
	p := prefixBP(seq, m, w.Circular)
	rt := m.RT()
	seqLen := len(seq)
	top := make([]StructureRecord, 0, n)

	for s := 0; s < seqLen; s++ {
		for mLen := w.MinLength; mLen <= maxLen(seqLen, s, w.Circular); mLen++ {
			g := windowEnergy(m, p, s, mLen)
			r := StructureRecord{Start: s, End: windowEnd(s, mLen, seqLen, w.Circular), Energy: g, Prob: math.Exp(-g/rt - ens.LogZ)}
			if len(top) == n && compareStructures(r, top[n-1]) >= 0 {
				continue
			}
			i := len(top)
			if i == n {
				i--
			} else {
				top = append(top, StructureRecord{})
			}
			for i > 0 && compareStructures(top[i-1], r) > 0 {
				top[i] = top[i-1]
				i--
			}
			top[i] = r
		}
	}

	return top
}
