package rlooper

type BasePairEnergies struct {
	rGG_dCC              float64
	rGC_dCG              float64
	rGA_dCT              float64
	rGU_dCA              float64
	rCG_dGC              float64
	rCC_dGG              float64
	rCA_dGT              float64
	rCU_dGA              float64
	rAG_dTC              float64
	rAC_dTG              float64
	rAA_dTT              float64
	rAU_dTA              float64
	rUG_dAC              float64
	rUC_dAG              float64
	rUA_dAT              float64
	rUU_dAA              float64
	homopolymer_override bool
	unconstrained        bool
	override_energy      float64
}

func NewBpEnergiesReasonableDefaults() BasePairEnergies {

	return BasePairEnergies{
		rGG_dCC:              -0.36,
		rGC_dCG:              -0.16,
		rGA_dCT:              -0.1,
		rGU_dCA:              -0.06,
		rCG_dGC:              0.97,
		rCC_dGG:              0.34,
		rCA_dGT:              0.45,
		rCU_dGA:              0.38,
		rAG_dTC:              -0.12,
		rAC_dTG:              -0.16,
		rAA_dTT:              0.6,
		rAU_dTA:              -0.12,
		rUG_dAC:              .45,
		rUC_dAG:              .5,
		rUA_dAT:              .28,
		rUU_dAA:              .8,
		homopolymer_override: false,
		unconstrained:        false,
		override_energy:      0.0,
	}
}
