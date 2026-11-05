package sidd

import (
	"bytes"
	"fmt"
	"math"
	"os"
	"strings"
	"testing"
)

func loadPBR322(t *testing.T) []byte {
	t.Helper()
	data, err := os.ReadFile("../testdata/fixtures/pbr322.toy.fa")
	if err != nil {
		t.Fatal(err)
	}
	var seq []byte
	for _, line := range bytes.Split(data, []byte("\n")) {
		if !bytes.HasPrefix(line, []byte(">")) {
			seq = append(seq, bytes.TrimSpace(line)...)
		}
	}
	if len(seq) != 4361 {
		t.Fatalf("fixture length = %d", len(seq))
	}
	return seq
}

func relDiff(a, b float64) float64 {
	return math.Abs(a-b) / math.Max(math.Abs(a), math.Abs(b))
}

func TestModelDefaults(t *testing.T) {
	p := DefaultParams()
	want := Params{Sigma: -0.06, Temperature: 310, Salt: 0.01, Threshold: 12, MinLen: 1, MaxLen: 250, MaxRuns: 4}
	if p != want {
		t.Fatalf("DefaultParams = %+v, want %+v", p, want)
	}
}

func TestModelPinnedCoefficients(t *testing.T) {
	p := DefaultParams()
	p.Geometry = GeomClamp
	m, err := NewModel(loadPBR322(t), p)
	if err != nil {
		t.Fatal(err)
	}
	for name, c := range map[string][2]string{
		"BAT":   {fmt.Sprintf("%.6g", m.BAT), "0.258116"},
		"BGC":   {fmt.Sprintf("%.6g", m.BGC), "1.30488"},
		"Alpha": {fmt.Sprintf("%.6g", m.Alpha), "-25.4481"},
	} {
		if c[0] != c[1] {
			t.Errorf("%s = %s, want %s", name, c[0], c[1])
		}
	}
}

func TestModelClampGeometry(t *testing.T) {
	seq := loadPBR322(t)
	p := DefaultParams()
	p.Geometry = GeomClamp
	m, err := NewModel(seq, p)
	if err != nil {
		t.Fatal(err)
	}
	if m.SeqLen != 4361 || m.EngLen != 4411 || len(m.Seq) != 4411 {
		t.Fatalf("lengths = %d/%d/%d", m.SeqLen, m.EngLen, len(m.Seq))
	}
	if string(m.Seq[m.SeqLen:]) != strings.Repeat("G", 50) {
		t.Errorf("pad = %q", m.Seq[m.SeqLen:])
	}
	if m.N != 4411 {
		t.Errorf("auto N = %v, want 4411", m.N)
	}

	p.N = 5000
	m, err = NewModel(seq, p)
	if err != nil {
		t.Fatal(err)
	}
	if m.N != 5000 || m.EngLen != 4411 {
		t.Errorf("override N = %v EngLen = %d", m.N, m.EngLen)
	}
}

func TestModelCircularGeometry(t *testing.T) {
	seq := loadPBR322(t)
	p := DefaultParams()
	p.Geometry = GeomCircular
	m, err := NewModel(seq, p)
	if err != nil {
		t.Fatal(err)
	}
	if m.EngLen != 4361 || m.N != 4361 {
		t.Errorf("EngLen = %d N = %v", m.EngLen, m.N)
	}
	p.MaxLen = 4361
	if _, err := NewModel(seq, p); err == nil {
		t.Error("expected error when MaxLen >= EngLen")
	}
}

func TestModelValidation(t *testing.T) {
	seq := []byte("ACGTACGTAC")
	cases := map[string]func(*Params) []byte{
		"linear NN":            func(p *Params) []byte { p.NearestNeighbor = true; return seq },
		"empty":                func(p *Params) []byte { return nil },
		"min zero":             func(p *Params) []byte { p.MinLen = 0; return seq },
		"max<min":              func(p *Params) []byte { p.MinLen, p.MaxLen = 5, 4; return seq },
		"negative N":           func(p *Params) []byte { p.N = -1; return seq },
		"bad base":             func(p *Params) []byte { return []byte("ACGX") },
		"clamp max":            func(p *Params) []byte { p.Geometry = GeomClamp; p.MaxLen = 60; return seq },
		"positive sigma":       func(p *Params) []byte { p.Sigma = 0.06; return seq },
		"NaN sigma":            func(p *Params) []byte { p.Sigma = math.NaN(); return seq },
		"infinite temperature": func(p *Params) []byte { p.Temperature = math.Inf(1); return seq },
		"NaN salt":             func(p *Params) []byte { p.Salt = math.NaN(); return seq },
		"infinite threshold":   func(p *Params) []byte { p.Threshold = math.Inf(1); return seq },
		"negative threshold":   func(p *Params) []byte { p.Threshold = -1; return seq },
		"min>seq":              func(p *Params) []byte { p.MinLen, p.MaxLen = 11, 20; return seq },
	}
	for name, mutate := range cases {
		p := DefaultParams()
		s := mutate(&p)
		if _, err := NewModel(s, p); err == nil {
			t.Errorf("%s: expected error", name)
		}
	}
}

func TestModelGResMatchesGround(t *testing.T) {
	p := DefaultParams()
	p.Geometry = GeomClamp
	m, err := NewModel(loadPBR322(t), p)
	if err != nil {
		t.Fatal(err)
	}
	if d := relDiff(m.GRes(0), m.GroundEnergy()); d >= 1e-12 {
		t.Errorf("GRes(0) = %v, GroundEnergy = %v, rel diff %g", m.GRes(0), m.GroundEnergy(), d)
	}
	if m.GRes(10) == m.GRes(0) {
		t.Error("GRes should vary with n")
	}
}

func TestModelCopolymericWindow(t *testing.T) {
	p := DefaultParams()
	p.Geometry = GeomCircular
	p.MaxLen = 20
	m, err := NewModel([]byte("AAAAAAAAAACCCCCCCCCCATGC"), p)
	if err != nil {
		t.Fatal(err)
	}
	if got := m.WindowEnergy(0, 10); got != 10*m.BAT {
		t.Errorf("all-A window = %v, want %v", got, 10*m.BAT)
	}
	if got := m.WindowEnergy(10, 10); got != 10*m.BGC {
		t.Errorf("all-C window = %v, want %v", got, 10*m.BGC)
	}
	if got, want := m.WindowEnergy(20, 6), 4*m.BAT+2*m.BGC; relDiff(got, want) > 1e-12 {
		t.Errorf("wrapped window = %v, want ~%v", got, want)
	}
}

func TestModelNearestNeighbor(t *testing.T) {
	seq := []byte("ACGTTGCAAGCT")
	p := DefaultParams()
	p.Geometry = GeomCircular
	p.NearestNeighbor = true
	p.MaxLen = 8
	m, err := NewModel(seq, p)
	if err != nil {
		t.Fatal(err)
	}
	dh := [4][4]float64{{32.3, 35.8, 36.0, 31.2}, {35.8, 35.1, 39.8, 36.0}, {36.0, 39.8, 35.1, 35.8}, {31.2, 36.0, 35.8, 32.3}}
	ds := [4][4]float64{{95.6, 100.4, 102.4, 95.3}, {100.4, 93.9, 102.4, 102.4}, {102.4, 102.4, 93.9, 100.4}, {95.3, 102.4, 100.4, 95.6}}
	code := map[byte]int{'A': 0, 'C': 1, 'G': 2, 'T': 3}
	dg := func(a, b byte) float64 {
		h := dh[code[a]][code[b]] / 4.184
		s := ds[code[a]][code[b]] / 4184.0
		s = 1 / (16.6*math.Log10(p.Salt/0.1)/h + 1/s)
		return h - p.Temperature*s
	}
	n := len(seq)
	sum := 0.0
	for i := 0; i < n; i++ {
		left, right := seq[(i+n-1)%n], seq[(i+1)%n]
		want := (dg(left, seq[i]) + dg(seq[i], right)) / 2
		if relDiff(m.perBase[i], want) >= 1e-12 {
			t.Errorf("perBase[%d] = %v, want %v", i, m.perBase[i], want)
		}
	}
	for _, i := range []int{10, 11, 0, 1} {
		sum += m.perBase[i]
	}
	if got := m.WindowEnergy(10, 4); got != sum {
		t.Errorf("wrapped NN window = %v, want %v", got, sum)
	}
}

func TestModelPositiveSigmaMessage(t *testing.T) {
	p := DefaultParams()
	p.Sigma = 0.06
	_, err := NewModel([]byte("ACGTACGTAC"), p)
	if err == nil || !strings.Contains(err.Error(), "negative --sigma") {
		t.Fatalf("error = %v", err)
	}
}

func TestModelLinearMaxLenCap(t *testing.T) {
	seq := randomSeq(400, 4)
	p := DefaultParams()
	p.Geometry = GeomLinear
	p.MaxLen = 1000000
	p.MaxRuns = 1
	m, err := NewModel(seq, p)
	if err != nil {
		t.Fatal(err)
	}
	if m.MaxLen != len(seq) {
		t.Errorf("MaxLen = %d, want %d", m.MaxLen, len(seq))
	}
	res, err := Run(m)
	if err != nil {
		t.Fatal(err)
	}
	if len(res.P) != len(seq) {
		t.Errorf("len(P) = %d, want %d", len(res.P), len(seq))
	}
}
