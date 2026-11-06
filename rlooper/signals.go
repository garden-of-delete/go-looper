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

// newSignalAccum returns a zeroed per-worker accumulator; the energy-weighted numerator is allocated only when withAvgG is set.
func newSignalAccum(n int, withAvgG bool) signalAccum {
	acc := signalAccum{bp: make([]float64, n)}
	if withAvgG {
		acc.num = make([]float64, n)
	}
	return acc
}

// addRow adds the probability (and energy-weighted probability) of every window starting at s to bases s..s+mLen-2 (mLen-1 bases; the window's last base is never credited), using a running sum from the longest window down so each base is written once.
func (a *signalAccum) addRow(m *EnergyModel, p []float64, s, minLength int, rt, logZ float64, circular bool) {
	n := len(a.bp)
	var sumProb, sumProbG float64
	for mLen := maxLen(n, s, circular); mLen >= 2; mLen-- {
		if mLen >= minLength {
			g := windowEnergy(m, p, s, mLen)
			prob := math.Exp(-g/rt - logZ)
			sumProb += prob
			sumProbG += prob * g
		}
		idx := s + mLen - 2
		if circular {
			idx %= n
		}
		a.bp[idx] += sumProb
		if a.num != nil {
			a.num[idx] += sumProbG
		}
	}
}

// merge adds another worker's partial signal sums into a.
func (a *signalAccum) merge(o signalAccum) {
	for i := range a.bp {
		a.bp[i] += o.bp[i]
	}
	for i := range a.num {
		a.num[i] += o.num[i]
	}
}

// ComputeSignals returns the per-base base-pairing probability, average energy and MFE-membership tracks for the ensemble ens.
func ComputeSignals(seq []byte, m *EnergyModel, w Windows, ens EnsembleResult, withAvgG bool) Signals {
	n := len(seq)
	p := prefixBP(seq, m, w.Circular)
	rt := m.RT()
	workers := clampWorkers(w.Workers, n)

	acc := newSignalAccum(n, withAvgG)
	if workers <= 1 {
		for s := 0; s < n; s++ {
			acc.addRow(m, p, s, w.MinLength, rt, ens.LogZ, w.Circular)
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
					chunk.addRow(m, p, s, w.MinLength, rt, ens.LogZ, w.Circular)
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
		length := ens.MFE.Length(n)
		for j := 0; j <= length-2; j++ {
			sig.MFE[(ens.MFE.Start+j)%n] = 1
		}
	}

	return sig
}

// Length returns the number of bases spanned by r in a sequence of length n, accounting for circular wrap when End < Start.
func (r StructureRecord) Length(n int) int {
	if r.End < r.Start {
		return r.End - r.Start + 1 + n
	}
	return r.End - r.Start + 1
}

// StructureResiduals returns the twist and linking-difference residuals of a single structure spanning mLen bases.
func StructureResiduals(m *EnergyModel, mLen int) (twist, linkingDifference float64) {
	l := float64(mLen - 1)
	fourPiSqC := 4 * math.Pi * math.Pi * m.C
	denom := fourPiSqC + m.K*l
	linkingDifference = (fourPiSqC / denom) * (m.Alpha + l*m.A)
	twist = (2 * math.Pi * m.K / denom) * (m.Alpha + l*m.A)
	return twist, linkingDifference
}

type residualAccum struct {
	twist, linkingDifference, wrappingAbsorption, strandSeparation float64
}

// addRow adds the probability-weighted residual contributions of every window starting at s.
func (a *residualAccum) addRow(m *EnergyModel, p []float64, s, n, minLength int, rt, logZ float64, circular bool) {
	for mLen := minLength; mLen <= maxLen(n, s, circular); mLen++ {
		g := windowEnergy(m, p, s, mLen)
		prob := math.Exp(-g/rt - logZ)
		l := float64(mLen - 1)
		rtw, rld := StructureResiduals(m, mLen)
		a.twist += rtw * prob
		a.linkingDifference += rld * prob
		a.wrappingAbsorption += (m.Alpha - rld + l*m.A) * prob
		a.strandSeparation -= l * m.A * prob
	}
}

// merge adds another worker's partial residual sums into a.
func (a *residualAccum) merge(o residualAccum) {
	a.twist += o.twist
	a.linkingDifference += o.linkingDifference
	a.wrappingAbsorption += o.wrappingAbsorption
	a.strandSeparation += o.strandSeparation
}

// ComputeResiduals returns the ensemble-averaged twist, linking difference, superhelicity, wrapping and strand-separation residuals, including the ground-state contribution.
func ComputeResiduals(seq []byte, m *EnergyModel, w Windows, ens EnsembleResult) Residuals {
	var res Residuals

	if ens.NumStructures > 0 {
		n := len(seq)
		p := prefixBP(seq, m, w.Circular)
		rt := m.RT()
		workers := clampWorkers(w.Workers, n)

		var acc residualAccum
		if workers <= 1 {
			for s := 0; s < n; s++ {
				acc.addRow(m, p, s, n, w.MinLength, rt, ens.LogZ, w.Circular)
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
						chunk.addRow(m, p, s, n, w.MinLength, rt, ens.LogZ, w.Circular)
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
