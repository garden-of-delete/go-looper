package rlooper

import (
	"math"
)

type ModelParams struct {
	N                   float64
	A                   float64
	C                   float64
	T                   float64
	k                   float64
	a                   float64
	sigma               float64
	alpha               float64
	bpEnergies          BasePairEnergies
	homopolymerOverride bool
	overrideEnergy      float64
}

func NewParamsReasonableDefaults() ModelParams {

	p := ModelParams{
		N:     1500,
		A:     1 / 10.4,
		C:     1.8,
		T:     310,
		a:     10,
		sigma: -0.07,
	}

	p.k = (2200 * 0.0019858775 * p.T) / p.N
	p.alpha = p.N * p.sigma * p.A
	p.bpEnergies = NewBpEnergiesReasonableDefaults()
	p.homopolymerOverride = false
	return p
}

// setters should be used where provided to update related biophysical quantities
func (p *ModelParams) SetN(N float64) {
	p.N = N
	p.alpha = p.N * p.sigma * p.A
	p.k = (2200 * 0.0019858775 * p.T) / p.N
}

func (p *ModelParams) SetA(A float64) {
	p.A = A
	p.alpha = p.N * p.sigma * p.A
}

func (p *ModelParams) SetSuperhelicity(sigma float64) {
	p.sigma = sigma
	p.alpha = p.N * p.sigma * p.A
}

func (p *ModelParams) SetT(T float64) {
	p.T = T
	p.k = (2200 * 0.0019858775 * p.T) / p.N
}

func (p *ModelParams) SetHomopolymerOverride(energy float64) {
	p.homopolymerOverride = true
	p.overrideEnergy = energy
}

func computeBoltzmannFactor(E float64, T float64) float64 {
	R := 0.0019858775
	return math.Exp(-1 * E / R * T)
}

// computeBpsInterval takes two characters representing DNA bases and returns the energy
// differential between the rloop and non rloop states in terms of energy.
func (p *ModelParams) computeBpsInterval(first rune, second rune) float64 {
	if p.homopolymerOverride {
		return p.overrideEnergy
	}
	if first == 'C' {
		if second == 'C' { //CC
			return p.bpEnergies.rGG_dCC
		} else if second == 'G' { //CG
			return p.bpEnergies.rCG_dGC
		} else if second == 'T' { //CT
			return p.bpEnergies.rGA_dCT
		} else { //CA
			return p.bpEnergies.rGU_dCA
		}
	} else if first == 'G' {
		if second == 'C' { //GC
			return p.bpEnergies.rCG_dGC
		} else if second == 'G' { //GG
			return p.bpEnergies.rCC_dGG
		} else if second == 'T' { //GT
			return p.bpEnergies.rCA_dGT
		} else { //GA
			return p.bpEnergies.rCU_dGA
		}
	} else if first == 'T' {
		if second == 'C' { //TC
			return p.bpEnergies.rAG_dTC
		} else if second == 'G' { //TG
			return p.bpEnergies.rAC_dTG
		} else if second == 'T' { //TT
			return p.bpEnergies.rAA_dTT
		} else { //TA
			return p.bpEnergies.rAU_dTA
		}
	} else {
		if second == 'C' { //AC
			return p.bpEnergies.rUG_dAC
		} else if second == 'G' { //AG
			return p.bpEnergies.rUC_dAG
		} else if second == 'T' { //AT
			return p.bpEnergies.rUA_dAT
		} else { //AA
			return p.bpEnergies.rUU_dAA
		}
	}
}

// ComputeStructure computes the free energy of a structure
// handles structures that cross circular boundaries automatically
func (p *ModelParams) ComputeStructure(seq []rune, w Window, structure *Structure) {
	var nBases int
	if w.End > w.Start { // if structure doesn't cross circular boundary
		nBases = w.End - w.Start + 1
	} else { // else structure includes the boundary of a circular dna piece
		nBases = len(seq) - w.Start + w.End // TODO: test
	}

	freeEnergy := 2 * math.Pow(math.Pi, 2) * p.C * p.k * math.Pow(p.alpha+float64(nBases)+p.A, 2) /
		(4*math.Pow(math.Pi, 2)*p.C + p.k*float64(nBases))

	var bpEnergy float64
	for i := w.Start; i != w.End; { // TODO: test
		// Handle the last base pair separately to avoid index out of range
		if i == len(seq)-1 {
			b0, b1 := seq[i], seq[0]
			bpEnergy += p.computeBpsInterval(b0, b1)
			i = 0
		} else {
			b0, b1 := seq[i], seq[i+1]
			bpEnergy += p.computeBpsInterval(b0, b1)
			i++
		}
	}

	structure.FreeEnergy = freeEnergy + bpEnergy
	structure.BoltzmannFactor = computeBoltzmannFactor(freeEnergy, p.T)
}

func (p *ModelParams) GroundStateFactor() float64 {
	return computeBoltzmannFactor(p.k*math.Pow(p.alpha, 2)-p.a, p.T)
}

func (p *ModelParams) GroundStateEnergy() float64 {
	return p.k*math.Pow(p.alpha, 2) - p.a
}
