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

func windowEnergy(m *EnergyModel, p []float64, s, e int) float64 {
	return m.SupercoilingEnergy(e-s+1) + (p[e] - p[s])
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

func ComputeEnsemble(seq []byte, m *EnergyModel, minLength, workers int) EnsembleResult {
	p := PrefixBPEnergies(seq, m)
	rt := m.RT()
	n := len(seq)
	workers = clampWorkers(workers, n)

	acc := newEnsembleAccum()
	if workers <= 1 {
		for s := 0; s < n; s++ {
			for e := s + minLength - 1; e < n; e++ {
				acc.addWindow(s, e, windowEnergy(m, p, s, e), rt)
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
					for e := s + minLength - 1; e < n; e++ {
						chunk.addWindow(s, e, windowEnergy(m, p, s, e), rt)
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

func AllStructures(seq []byte, m *EnergyModel, minLength int, ens EnsembleResult) []StructureRecord {
	p := PrefixBPEnergies(seq, m)
	rt := m.RT()
	recs := make([]StructureRecord, 0, ens.NumStructures)

	for s := 0; s < len(seq); s++ {
		for e := s + minLength - 1; e < len(seq); e++ {
			g := windowEnergy(m, p, s, e)
			recs = append(recs, StructureRecord{
				Start:  s,
				End:    e,
				Energy: g,
				Prob:   math.Exp(-g/rt - ens.LogZ),
			})
		}
	}

	slices.SortFunc(recs, compareStructures)

	return recs
}
