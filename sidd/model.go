package sidd

import (
	"errors"
	"fmt"
	"math"
)

type Geometry int

const (
	GeomLinear Geometry = iota
	GeomClamp
	GeomCircular
)

const (
	C          = 3.6
	A          = 10.4
	Nucleation = 10.84
	Shift      = 300.0
	GMaxAdj    = 12.0

	clampPad = 50
	linking  = 2220.0
)

var pi = math.Pi

type Params struct {
	Sigma, Temperature, Salt, Threshold float64
	MinLen, MaxLen                      int
	N                                   float64
	NearestNeighbor                     bool
	Geometry                            Geometry
	MaxRuns                             int
}

func DefaultParams() Params {
	return Params{Sigma: -0.06, Temperature: 310, Salt: 0.01, Threshold: 12, MinLen: 1, MaxLen: 250, MaxRuns: 4}
}

type Model struct {
	Params
	Seq                       []byte
	SeqLen, EngLen            int
	N, RT, K, Alpha, BAT, BGC float64
	perBase                   []float64
}

var baseCode = [256]int8{'A': 1, 'C': 2, 'G': 3, 'T': 4}

func NewModel(seq []byte, p Params) (*Model, error) {
	switch {
	case len(seq) == 0:
		return nil, errors.New("sidd: empty sequence")
	case p.MinLen < 1:
		return nil, fmt.Errorf("sidd: minimum length %d must be at least 1", p.MinLen)
	case p.MaxLen < p.MinLen:
		return nil, fmt.Errorf("sidd: maximum length %d is below minimum length %d", p.MaxLen, p.MinLen)
	case p.MinLen > len(seq):
		return nil, fmt.Errorf("sidd: minimum length %d exceeds the sequence length %d", p.MinLen, len(seq))
	case p.NearestNeighbor && p.Geometry == GeomLinear:
		return nil, errors.New("sidd: nearest-neighbor energetics require clamp or circular geometry")
	case !isFinite(p.Sigma) || !isFinite(p.Temperature) || !isFinite(p.Salt) || !isFinite(p.Threshold):
		return nil, errors.New("sidd: sigma, temperature, salt and threshold must be finite numbers")
	case p.Salt <= 0 || p.Temperature <= 0:
		return nil, errors.New("sidd: salt and temperature must be positive")
	case p.Threshold < 0:
		return nil, fmt.Errorf("sidd: threshold %v must not be negative", p.Threshold)
	case p.Sigma > 0:
		return nil, fmt.Errorf("sidd: positive sigma %v (positive supercoiling) is not supported; SIST's positive-magnitude -s convention corresponds to a negative --sigma here", p.Sigma)
	case p.MaxRuns < 1 || p.MaxRuns > 4:
		return nil, fmt.Errorf("sidd: maximum runs %d must be between 1 and 4", p.MaxRuns)
	}
	for i, b := range seq {
		if baseCode[b] == 0 {
			return nil, fmt.Errorf("sidd: invalid base %q at position %d", b, i)
		}
	}

	m := &Model{Params: p, SeqLen: len(seq), EngLen: len(seq)}
	m.Seq = append([]byte(nil), seq...)
	if p.Geometry == GeomClamp {
		m.EngLen += clampPad
		for i := 0; i < clampPad; i++ {
			m.Seq = append(m.Seq, 'G')
		}
	}
	if p.Geometry != GeomLinear && p.MaxLen >= m.EngLen {
		return nil, fmt.Errorf("sidd: maximum length %d must be below the sequence length %d (plus %d clamp bases for clamp geometry)", p.MaxLen, len(seq), m.EngLen-len(seq))
	}
	if p.Geometry == GeomLinear {
		m.MaxLen = min(p.MaxLen, m.SeqLen)
	}

	m.N = p.N
	if m.N == 0 {
		m.N = float64(m.EngLen)
	}
	if m.N <= 0 {
		return nil, fmt.Errorf("sidd: N %v must be positive", m.N)
	}

	tmAT := 354.65 + 16.6*math.Log10(p.Salt)
	tmGC := tmAT + 41.0
	m.RT = 1.9872 * p.Temperature / 1000.0
	m.BAT = 7.2464 * (1.0 - p.Temperature/tmAT)
	m.BGC = 9.0172 * (1.0 - p.Temperature/tmGC)
	m.K = linking * m.RT / m.N
	m.Alpha = p.Sigma * m.N / A

	if p.NearestNeighbor {
		m.perBase = m.nearestNeighborEnergies()
	}
	return m, nil
}

func (m *Model) nearestNeighborEnergies() []float64 {
	enthalpy := [4][4]float64{
		{32.3, 35.8, 36.0, 31.2},
		{35.8, 35.1, 39.8, 36.0},
		{36.0, 39.8, 35.1, 35.8},
		{31.2, 36.0, 35.8, 32.3},
	}
	entropy := [4][4]float64{
		{95.6, 100.4, 102.4, 95.3},
		{100.4, 93.9, 102.4, 102.4},
		{102.4, 102.4, 93.9, 100.4},
		{95.3, 102.4, 100.4, 95.6},
	}
	var dg [4][4]float64
	for i := range dg {
		for j := range dg[i] {
			h := enthalpy[i][j] / 4.184
			s := entropy[i][j] / 4184.0
			s = 1.0 / (16.6*math.Log10(m.Salt/0.1)/h + 1.0/s)
			dg[i][j] = h - float64(m.Temperature*s)
		}
	}

	perBase := make([]float64, m.EngLen)
	for i := range perBase {
		left := baseCode[m.Seq[(i+m.EngLen-1)%m.EngLen]] - 1
		cur := baseCode[m.Seq[i]] - 1
		right := baseCode[m.Seq[(i+1)%m.EngLen]] - 1
		perBase[i] = (dg[left][cur] + dg[cur][right]) / 2.0
	}
	return perBase
}

func (m *Model) WindowEnergy(s, n int) float64 {
	if m.perBase != nil {
		sum := 0.0
		for i := s; i < s+n; i++ {
			sum += m.perBase[i%m.EngLen]
		}
		return sum
	}
	var at, gc int
	for i := s; i < s+n; i++ {
		switch m.Seq[i%m.EngLen] {
		case 'A', 'T':
			at++
		default:
			gc++
		}
	}
	return float64(m.BAT*float64(at)) + float64(m.BGC*float64(gc))
}

func (m *Model) GRes(n int) float64 {
	q := m.Alpha + float64(n)/A
	t1 := q * q
	t2 := 2 * pi * pi * C * m.K * t1
	return t2 / (4*pi*pi*C + float64(m.K*float64(n)))
}

func (m *Model) GroundEnergy() float64 {
	return 0.5 * m.K * m.Alpha * m.Alpha
}
