package sidd

import (
	"bytes"
	"math"
	"os"
	"testing"

	"golooper/internal/goldencmp"
)

// TestFormatSIST pins FormatSIST to C++ default float output, including nan and inf.
func TestFormatSIST(t *testing.T) {
	tests := []struct {
		in   float64
		want string
	}{
		{0.2581159999, "0.258116"},
		{1.304881234, "1.30488"},
		{-25.448076923076923, "-25.4481"},
		{4.347821234e-10, "4.34782e-10"},
		{1.00410004e-06, "1.0041e-06"},
		{1.01700001e-08, "1.017e-08"},
		{0.000115, "0.000115"},
		{0.0021658712, "0.00216587"},
		{11.98304, "11.983"},
		{12, "12"},
		{0, "0"},
		{10.84, "10.84"},
		{-0.06, "-0.06"},
		{4908539.0, "4.90854e+06"},
		{123456.4, "123456"},
		{1e100, "1e+100"},
		{2.226931e-60, "2.22693e-60"},
		{2769779, "2.76978e+06"},
		{2.191011e-60, "2.19101e-60"},
		{300, "300"},
		{math.NaN(), "nan"},
		{math.Inf(1), "inf"},
		{math.Inf(-1), "-inf"},
	}
	for _, tt := range tests {
		if got := FormatSIST(tt.in); got != tt.want {
			t.Errorf("FormatSIST(%v) = %q, want %q", tt.in, got, tt.want)
		}
	}
}

// assertReportParity asserts a single-run report is byte-identical to the SIST golden stdout for goldenCase.
func assertReportParity(t *testing.T, goldenCase string, mutate func(*Params)) {
	t.Helper()
	p := DefaultParams()
	p.MaxRuns = 1
	mutate(&p)
	m, err := NewModel(loadPBR322(t), p)
	if err != nil {
		t.Fatal(err)
	}
	res, err := Run(m)
	if err != nil {
		t.Fatal(err)
	}
	var got bytes.Buffer
	if err := WriteReport(&got, "pbr322.toy.fa", m, res); err != nil {
		t.Fatal(err)
	}
	want, err := os.ReadFile("../testdata/golden/sidd/" + goldenCase + "/stdout.txt")
	if err != nil {
		t.Fatal(err)
	}
	if bytes.Equal(got.Bytes(), want) {
		return
	}
	gotLines := bytes.Split(got.Bytes(), []byte("\n"))
	wantLines := bytes.Split(want, []byte("\n"))
	for i := 0; i < len(gotLines) && i < len(wantLines); i++ {
		if !bytes.Equal(gotLines[i], wantLines[i]) {
			t.Fatalf("line %d differs:\n golden: %q\n    got: %q", i+1, wantLines[i], gotLines[i])
		}
	}
	t.Fatalf("line counts differ: golden %d, got %d", len(wantLines), len(gotLines))
}

// TestReportGoldenParity_ClampCopolyDefaults pins single-run byte parity for clamp copolymeric defaults.
func TestReportGoldenParity_ClampCopolyDefaults(t *testing.T) {
	assertReportParity(t, "clamp-copoly-defaults", func(p *Params) {
		p.Geometry = GeomClamp
	})
}

// TestReportGoldenParity_CircularCopolyDefaults pins single-run byte parity for circular copolymeric defaults.
func TestReportGoldenParity_CircularCopolyDefaults(t *testing.T) {
	assertReportParity(t, "circular-copoly-defaults", func(p *Params) {
		p.Geometry = GeomCircular
	})
}

// TestReportGoldenParity_ClampCopolyAlt pins single-run byte parity for clamp copolymeric with non-default temperature, salt, sigma and threshold.
func TestReportGoldenParity_ClampCopolyAlt(t *testing.T) {
	assertReportParity(t, "clamp-copoly-alt", func(p *Params) {
		p.Geometry = GeomClamp
		p.Temperature = 315
		p.Salt = 0.1
		p.Sigma = -0.04
		p.Threshold = 10
	})
}

// assertMultiRunParity asserts a multi-run report matches the patched-SIST golden within a 1e-6 relative float tolerance.
func assertMultiRunParity(t *testing.T, goldenCase, fixture string, mutate func(*Params)) {
	t.Helper()
	p := DefaultParams()
	p.MaxRuns = 4
	mutate(&p)
	data, err := os.ReadFile("../testdata/fixtures/" + fixture)
	if err != nil {
		t.Fatal(err)
	}
	var seq []byte
	for _, line := range bytes.Split(data, []byte("\n")) {
		if !bytes.HasPrefix(line, []byte(">")) {
			seq = append(seq, bytes.TrimSpace(line)...)
		}
	}
	m, err := NewModel(seq, p)
	if err != nil {
		t.Fatal(err)
	}
	res, err := Run(m)
	if err != nil {
		t.Fatal(err)
	}
	var got bytes.Buffer
	if err := WriteReport(&got, fixture, m, res); err != nil {
		t.Fatal(err)
	}
	want, err := os.ReadFile("../testdata/golden/sidd-multirun/" + goldenCase + "/stdout.txt")
	if err != nil {
		t.Fatal(err)
	}
	for _, mm := range goldencmp.Compare(want, got.Bytes(), goldencmp.Tolerances{FloatRel: 1e-6, MinMagnitude: 1e-12}) {
		t.Error(mm)
	}
}

// TestReportMultiRunParity_ClampCopolyDefaults pins multi-run parity for clamp copolymeric defaults.
func TestReportMultiRunParity_ClampCopolyDefaults(t *testing.T) {
	assertMultiRunParity(t, "mr-clamp-copoly-defaults", "pbr322.toy.fa", func(p *Params) {
		p.Geometry = GeomClamp
	})
}

// TestReportMultiRunParity_Hot pins multi-run parity at elevated temperature and threshold.
func TestReportMultiRunParity_Hot(t *testing.T) {
	assertMultiRunParity(t, "mr-clamp-copoly-hot", "pbr322.toy.fa", func(p *Params) {
		p.Geometry = GeomClamp
		p.Temperature = 320
		p.Threshold = 15
		p.Sigma = -0.04
	})
}

// TestReportMultiRunParity_ClampNNDefaults pins multi-run parity for clamp nearest-neighbor defaults.
func TestReportMultiRunParity_ClampNNDefaults(t *testing.T) {
	assertMultiRunParity(t, "mr-clamp-nn-defaults", "pbr322.toy.fa", func(p *Params) {
		p.Geometry = GeomClamp
		p.NearestNeighbor = true
	})
}

// TestReportMultiRunParity_CircularCopolyDefaults pins multi-run parity for circular copolymeric defaults.
func TestReportMultiRunParity_CircularCopolyDefaults(t *testing.T) {
	assertMultiRunParity(t, "mr-circular-copoly-defaults", "pbr322.toy.fa", func(p *Params) {
		p.Geometry = GeomCircular
	})
}

// TestReportMultiRunParity_CircularNNSynth pins multi-run parity for circular nearest-neighbor on the synthetic 1600 sequence.
func TestReportMultiRunParity_CircularNNSynth(t *testing.T) {
	assertMultiRunParity(t, "mr-circular-nn-synth", "sidd_synth1600.fa", func(p *Params) {
		p.Geometry = GeomCircular
		p.NearestNeighbor = true
	})
}

// TestReportMultiRunParity_CircularNNDefaults pins multi-run parity for circular nearest-neighbor defaults.
func TestReportMultiRunParity_CircularNNDefaults(t *testing.T) {
	assertMultiRunParity(t, "mr-circular-nn-defaults", "pbr322.toy.fa", func(p *Params) {
		p.Geometry = GeomCircular
		p.NearestNeighbor = true
	})
}

// TestReportMultiRunParity_ClampCopolySynth pins multi-run parity for clamp copolymeric on the synthetic 1600 sequence.
func TestReportMultiRunParity_ClampCopolySynth(t *testing.T) {
	assertMultiRunParity(t, "mr-clamp-copoly-synth", "sidd_synth1600.fa", func(p *Params) {
		p.Geometry = GeomClamp
	})
}

// TestReportMultiRunParity_OriginStraddle pins multi-run parity on a rotated sequence so runs straddle the origin.
func TestReportMultiRunParity_OriginStraddle(t *testing.T) {
	assertMultiRunParity(t, "mr-origin-straddle", "pbr322_rot310.toy.fa", func(p *Params) {
		p.Geometry = GeomCircular
	})
}

// TestReportMultiRunParity_ClampCopolyAlt pins multi-run parity for clamp copolymeric with non-default parameters.
func TestReportMultiRunParity_ClampCopolyAlt(t *testing.T) {
	assertMultiRunParity(t, "mr-clamp-copoly-alt", "pbr322.toy.fa", func(p *Params) {
		p.Geometry = GeomClamp
		p.Temperature = 315
		p.Salt = 0.1
		p.Sigma = -0.04
		p.Threshold = 10
	})
}
