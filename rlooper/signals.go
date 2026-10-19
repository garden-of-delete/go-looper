package rlooper

import (
	"math"
	"sync"
)

type Signals struct {
	BPProb, AvgG, MFE []float64
}

type Residuals struct {
	Twist, LinkingDifference, Superhelicity, WrappingAbsorption, StrandSeparation float64
}

type signalAccum struct {
	bp, num []float64
}

func newSignalAccum(n int, withAvgG bool) signalAccum {
	acc := signalAccum{bp: make([]float64, n)}
	if withAvgG {
		acc.num = make([]float64, n)
	}
	return acc
}

func (a *signalAccum) addRow(m *EnergyModel, p []float64, s, minLength int, rt, logZ float64) {
	n := len(a.bp)
	first := s + minLength - 1
	var sumProb, sumProbG float64
	for e := n - 1; e > s; e-- {
		if e >= first {
			g := windowEnergy(m, p, s, e)
			prob := math.Exp(-g/rt - logZ)
			sumProb += prob
			sumProbG += prob * g
		}
		a.bp[e-1] += sumProb
		if a.num != nil {
			a.num[e-1] += sumProbG
		}
	}
}

func (a *signalAccum) merge(o signalAccum) {
	for i := range a.bp {
		a.bp[i] += o.bp[i]
	}
	for i := range a.num {
		a.num[i] += o.num[i]
	}
}

func ComputeSignals(seq []byte, m *EnergyModel, minLength int, ens EnsembleResult, withAvgG bool, workers int) Signals {
	n := len(seq)
	p := PrefixBPEnergies(seq, m)
	rt := m.RT()
	workers = clampWorkers(workers, n)

	acc := newSignalAccum(n, withAvgG)
	if workers <= 1 {
		for s := 0; s < n; s++ {
			acc.addRow(m, p, s, minLength, rt, ens.LogZ)
		}
	} else {
		results := make([]signalAccum, workers)
		var wg sync.WaitGroup
		for c := 0; c < workers; c++ {
			wg.Add(1)
			go func(c int) {
				defer wg.Done()
				chunk := newSignalAccum(n, withAvgG)
				for s := c; s < n; s += workers {
					chunk.addRow(m, p, s, minLength, rt, ens.LogZ)
				}
				results[c] = chunk
			}(c)
		}
		wg.Wait()
		for _, r := range results {
			acc.merge(r)
		}
	}

	sig := Signals{BPProb: acc.bp, MFE: make([]float64, n)}
	if withAvgG {
		sig.AvgG = make([]float64, n)
		for i, den := range acc.bp {
			if den > 0 {
				sig.AvgG[i] = acc.num[i] / den
			}
		}
	}

	if ens.NumStructures > 0 {
		for i := ens.MFE.Start; i < ens.MFE.End; i++ {
			sig.MFE[i] = 1
		}
	}

	return sig
}

func StructureResiduals(m *EnergyModel, r StructureRecord) (twist, linkingDifference float64) {
	l := float64(r.End - r.Start)
	fourPiSqC := 4 * math.Pi * math.Pi * m.C
	denom := fourPiSqC + m.K*l
	linkingDifference = (fourPiSqC / denom) * (m.Alpha + l*m.A)
	twist = (2 * math.Pi * m.K / denom) * (m.Alpha + l*m.A)
	return twist, linkingDifference
}

type residualAccum struct {
	twist, linkingDifference, wrappingAbsorption, strandSeparation float64
}

func (a *residualAccum) addRow(m *EnergyModel, p []float64, s, n, minLength int, rt, logZ float64) {
	for e := s + minLength - 1; e < n; e++ {
		g := windowEnergy(m, p, s, e)
		prob := math.Exp(-g/rt - logZ)
		l := float64(e - s)
		rtw, rld := StructureResiduals(m, StructureRecord{Start: s, End: e})
		a.twist += rtw * prob
		a.linkingDifference += rld * prob
		a.wrappingAbsorption += (m.Alpha - rld + l*m.A) * prob
		a.strandSeparation -= l * m.A * prob
	}
}

func (a *residualAccum) merge(o residualAccum) {
	a.twist += o.twist
	a.linkingDifference += o.linkingDifference
	a.wrappingAbsorption += o.wrappingAbsorption
	a.strandSeparation += o.strandSeparation
}

func ComputeResiduals(seq []byte, m *EnergyModel, minLength int, ens EnsembleResult, workers int) Residuals {
	var res Residuals

	if ens.NumStructures > 0 {
		n := len(seq)
		p := PrefixBPEnergies(seq, m)
		rt := m.RT()
		workers = clampWorkers(workers, n)

		var acc residualAccum
		if workers <= 1 {
			for s := 0; s < n; s++ {
				acc.addRow(m, p, s, n, minLength, rt, ens.LogZ)
			}
		} else {
			results := make([]residualAccum, workers)
			var wg sync.WaitGroup
			for c := 0; c < workers; c++ {
				wg.Add(1)
				go func(c int) {
					defer wg.Done()
					var chunk residualAccum
					for s := c; s < n; s += workers {
						chunk.addRow(m, p, s, n, minLength, rt, ens.LogZ)
					}
					results[c] = chunk
				}(c)
			}
			wg.Wait()
			for _, r := range results {
				acc.merge(r)
			}
		}

		res.Twist = acc.twist
		res.LinkingDifference = acc.linkingDifference
		res.WrappingAbsorption = acc.wrappingAbsorption
		res.StrandSeparation = acc.strandSeparation
	}

	res.Twist += (m.K * m.C * m.Alpha / (2 * math.Pi)) * ens.PGround
	res.LinkingDifference += m.Alpha * ens.PGround
	res.Superhelicity = res.LinkingDifference / (m.N * m.A)

	return res
}
