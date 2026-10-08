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
	bpDiff, numDiff, denDiff []float64
}

func newSignalAccum(n int, withAvgG bool) signalAccum {
	acc := signalAccum{bpDiff: make([]float64, n+1)}
	if withAvgG {
		acc.numDiff = make([]float64, n+1)
		acc.denDiff = make([]float64, n+1)
	}
	return acc
}

func (a *signalAccum) addWindow(s, e int, g, prob float64, withAvgG bool) {
	a.bpDiff[s] += prob
	a.bpDiff[e] -= prob
	if withAvgG {
		a.numDiff[s] += prob * g
		a.numDiff[e] -= prob * g
		a.denDiff[s] += prob
		a.denDiff[e] -= prob
	}
}

func (a *signalAccum) merge(o signalAccum) {
	for i := range a.bpDiff {
		a.bpDiff[i] += o.bpDiff[i]
	}
	if a.numDiff == nil {
		return
	}
	for i := range a.numDiff {
		a.numDiff[i] += o.numDiff[i]
		a.denDiff[i] += o.denDiff[i]
	}
}

func ComputeSignals(seq []byte, m *EnergyModel, minLength int, ens EnsembleResult, withAvgG bool, workers int) Signals {
	n := len(seq)
	sig := Signals{BPProb: make([]float64, n), MFE: make([]float64, n)}
	if withAvgG {
		sig.AvgG = make([]float64, n)
	}

	p := PrefixBPEnergies(seq, m)
	rt := m.RT()
	workers = clampWorkers(workers, n)

	acc := newSignalAccum(n, withAvgG)
	if workers <= 1 {
		for s := 0; s < n; s++ {
			for e := s + minLength - 1; e < n; e++ {
				g := windowEnergy(m, p, s, e)
				prob := math.Exp(-g/rt - ens.LogZ)
				acc.addWindow(s, e, g, prob, withAvgG)
			}
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
					for e := s + minLength - 1; e < n; e++ {
						g := windowEnergy(m, p, s, e)
						prob := math.Exp(-g/rt - ens.LogZ)
						chunk.addWindow(s, e, g, prob, withAvgG)
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

	var bpAcc, numAcc, denAcc float64
	for i := 0; i < n; i++ {
		bpAcc += acc.bpDiff[i]
		sig.BPProb[i] = bpAcc
		if withAvgG {
			numAcc += acc.numDiff[i]
			denAcc += acc.denDiff[i]
			if denAcc > 0 {
				sig.AvgG[i] = numAcc / denAcc
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

func ComputeResiduals(seq []byte, m *EnergyModel, minLength int, ens EnsembleResult) Residuals {
	var res Residuals

	if ens.NumStructures > 0 {
		p := PrefixBPEnergies(seq, m)
		rt := m.RT()
		for s := 0; s < len(seq); s++ {
			for e := s + minLength - 1; e < len(seq); e++ {
				g := windowEnergy(m, p, s, e)
				prob := math.Exp(-g/rt - ens.LogZ)
				l := float64(e - s)
				rtw, rld := StructureResiduals(m, StructureRecord{Start: s, End: e})
				res.Twist += rtw * prob
				res.LinkingDifference += rld * prob
				res.WrappingAbsorption += (m.Alpha - rld + l*m.A) * prob
				res.StrandSeparation -= l * m.A * prob
			}
		}
	}

	res.Twist += (m.K * m.C * m.Alpha / (2 * math.Pi)) * ens.PGround
	res.LinkingDifference += m.Alpha * ens.PGround
	res.Superhelicity = res.LinkingDifference / (m.N * m.A)

	return res
}
