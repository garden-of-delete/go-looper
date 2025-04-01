package model

type Params struct {
	N          float64
	A          float64
	C          float64
	T          float64
	k          float64
	a          float64
	sigma      float64
	alpha      float64
	bpEnergies BasePairEnergies
}

func NewParamsReasonableDefaults() Params {

	p := Params{
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
	return p
}

// setters should be used where provided to update related biophysical quantities
func (p *Params) setN(N float64) {
	p.N = N
	p.alpha = p.N * p.sigma * p.A
	p.k = (2200 * 0.0019858775 * p.T) / p.N
}

func (p *Params) setA(A float64) {
	p.A = A
	p.alpha = p.N * p.sigma * p.A
}

func (p *Params) setSuperhelicity(sigma float64) {
	p.sigma = sigma
	p.alpha = p.N * p.sigma * p.A
}

func (p *Params) setT(T float64) {
	p.T = T
	p.k = (2200 * 0.0019858775 * p.T) / p.N
}
