package sidd

import (
	"bytes"
	"os"
	"testing"
)

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
	}
	for _, tt := range tests {
		if got := FormatSIST(tt.in); got != tt.want {
			t.Errorf("FormatSIST(%v) = %q, want %q", tt.in, got, tt.want)
		}
	}
}

func assertReportParity(t *testing.T, goldenCase string, mutate func(*Params)) {
	t.Helper()
	p := DefaultParams()
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

func TestReportGoldenParity_ClampCopolyDefaults(t *testing.T) {
	assertReportParity(t, "clamp-copoly-defaults", func(p *Params) {
		p.Geometry = GeomClamp
	})
}

func TestReportGoldenParity_CircularCopolyDefaults(t *testing.T) {
	assertReportParity(t, "circular-copoly-defaults", func(p *Params) {
		p.Geometry = GeomCircular
	})
}

func TestReportGoldenParity_ClampCopolyAlt(t *testing.T) {
	assertReportParity(t, "clamp-copoly-alt", func(p *Params) {
		p.Geometry = GeomClamp
		p.Temperature = 315
		p.Salt = 0.1
		p.Sigma = -0.04
		p.Threshold = 10
	})
}
