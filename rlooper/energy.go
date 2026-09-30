package rlooper

import "math"

const GasConstant = 0.0019858775

type EnergyModel struct {
	N, A, C, T, K, Nucleation, Sigma, Alpha float64
	BP                                      BasePairEnergies
	Homopolymer                             bool
	OverrideEnergy                          float64
	Unconstrained                           bool
}

func NewDefaultModel() *EnergyModel {
	m := &EnergyModel{N: 1500, A: 1 / 10.4, C: 1.8, T: 310, Nucleation: 10, Sigma: -0.07,
		BP: NewBpEnergiesReasonableDefaults()}
	m.recompute()
	return m
}

func (m *EnergyModel) recompute() {
	m.K = 2200 * GasConstant * m.T / m.N
	m.Alpha = m.N * m.Sigma * m.A
}

func (m *EnergyModel) SetN(n float64)     { m.N = n; m.recompute() }
func (m *EnergyModel) SetSigma(s float64) { m.Sigma = s; m.recompute() }
func (m *EnergyModel) SetT(t float64)     { m.T = t; m.recompute() }
func (m *EnergyModel) SetA(a float64)     { m.A = a; m.recompute() }
func (m *EnergyModel) RT() float64        { return GasConstant * m.T }

func (m *EnergyModel) SupercoilingEnergy(bases int) float64 {
	if m.Unconstrained {
		return 0
	}
	n := float64(bases)
	q := m.Alpha + n*m.A
	return 2 * math.Pi * math.Pi * m.C * m.K * q * q / (4*math.Pi*math.Pi*m.C + m.K*n)
}

func (m *EnergyModel) GroundStateEnergy() float64 {
	return m.K*m.Alpha*m.Alpha/2 - m.Nucleation
}

func (m *EnergyModel) StepEnergy(b0, b1 byte) float64 {
	if m.Homopolymer {
		return m.OverrideEnergy
	}
	switch b0 {
	case 'C':
		switch b1 {
		case 'C':
			return m.BP.rGG_dCC
		case 'G':
			return m.BP.rGC_dCG
		case 'T':
			return m.BP.rGA_dCT
		default:
			return m.BP.rGU_dCA
		}
	case 'G':
		switch b1 {
		case 'C':
			return m.BP.rCG_dGC
		case 'G':
			return m.BP.rCC_dGG
		case 'T':
			return m.BP.rCA_dGT
		default:
			return m.BP.rCU_dGA
		}
	case 'T':
		switch b1 {
		case 'C':
			return m.BP.rAG_dTC
		case 'G':
			return m.BP.rAC_dTG
		case 'T':
			return m.BP.rAA_dTT
		default:
			return m.BP.rAU_dTA
		}
	default:
		switch b1 {
		case 'C':
			return m.BP.rUG_dAC
		case 'G':
			return m.BP.rUC_dAG
		case 'T':
			return m.BP.rUA_dAT
		default:
			return m.BP.rUU_dAA
		}
	}
}

func PrefixBPEnergies(seq []byte, m *EnergyModel) []float64 {
	p := make([]float64, len(seq))
	for i := 1; i < len(seq); i++ {
		p[i] = p[i-1] + m.StepEnergy(seq[i-1], seq[i])
	}
	return p
}
