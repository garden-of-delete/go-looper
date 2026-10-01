package rlooper

import (
	"cmp"
	"math"
	"slices"
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

func ComputeEnsemble(seq []byte, m *EnergyModel, minLength, workers int) EnsembleResult {
	p := PrefixBPEnergies(seq, m)
	rt := m.RT()
	lse := newLogSumExp()
	var mfe StructureRecord
	numStructures := 0

	for s := 0; s < len(seq); s++ {
		for e := s + minLength - 1; e < len(seq); e++ {
			g := windowEnergy(m, p, s, e)
			lse.add(-g / rt)
			rec := StructureRecord{Start: s, End: e, Energy: g}
			if numStructures == 0 || compareStructures(rec, mfe) < 0 {
				mfe = rec
			}
			numStructures++
		}
	}

	logGround := -m.GroundStateEnergy() / rt
	lse.add(logGround)
	logZ := lse.value()

	if numStructures > 0 {
		mfe.Prob = math.Exp(-mfe.Energy/rt - logZ)
	}

	return EnsembleResult{
		LogZ:          logZ,
		LogGround:     logGround,
		PGround:       math.Exp(logGround - logZ),
		MFE:           mfe,
		NumStructures: numStructures,
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
