package cmd

import "testing"

// TestParseSuperhelicityDomain pins that --N accepts "auto" or a positive integer and rejects zero, negatives, non-digits, decimals and the empty string.
func TestParseSuperhelicityDomain(t *testing.T) {
	cases := []struct {
		in        string
		wantAuto  bool
		wantValue int
		wantErr   bool
	}{
		{"auto", true, 0, false},
		{"1500", false, 1500, false},
		{"1", false, 1, false},
		{"0", false, 0, true},
		{"-5", false, 0, true},
		{"15OO", false, 0, true},
		{"2000.5", false, 0, true},
		{"", false, 0, true},
	}

	for _, c := range cases {
		auto, value, err := parseSuperhelicityDomain(c.in)
		if (err != nil) != c.wantErr {
			t.Errorf("parseSuperhelicityDomain(%q) error = %v, wantErr %v", c.in, err, c.wantErr)
			continue
		}
		if c.wantErr {
			continue
		}
		if auto != c.wantAuto || value != c.wantValue {
			t.Errorf("parseSuperhelicityDomain(%q) = (%v, %d), want (%v, %d)", c.in, auto, value, c.wantAuto, c.wantValue)
		}
	}
}

// TestValidateMinLength pins the boundary at 2: 2 and above pass, 1 and below fail.
func TestValidateMinLength(t *testing.T) {
	cases := []struct {
		n       int
		wantErr bool
	}{
		{2, false},
		{3, false},
		{1, true},
		{0, true},
		{-1, true},
	}

	for _, c := range cases {
		if err := validateMinLength(c.n); (err != nil) != c.wantErr {
			t.Errorf("validateMinLength(%d) error = %v, wantErr %v", c.n, err, c.wantErr)
		}
	}
}

// TestCircularFlagAccepted pins that --circular parses and sets cfg.Circular.
func TestCircularFlagAccepted(t *testing.T) {
	defer func() { cfg.Circular = false }()
	rootCmd.SetArgs([]string{"show-config", "--circular"})
	defer rootCmd.SetArgs(nil)

	if err := rootCmd.Execute(); err != nil {
		t.Fatalf("Execute with --circular = %v, want nil", err)
	}
	if !cfg.Circular {
		t.Error("cfg.Circular = false after --circular, want true")
	}
}

// TestValidateSigmaRange pins rejection of zero and negative steps and of min > max, and acceptance of min == max.
func TestValidateSigmaRange(t *testing.T) {
	cases := []struct {
		name           string
		min, max, step float64
		wantErr        bool
	}{
		{"ok", -0.10, -0.02, 0.01, false},
		{"zero step", -0.10, -0.02, 0, true},
		{"negative step", -0.10, -0.02, -0.01, true},
		{"min greater than max", -0.02, -0.10, 0.01, true},
		{"min equals max", -0.05, -0.05, 0.01, false},
	}

	for _, c := range cases {
		if err := validateSigmaRange(c.min, c.max, c.step); (err != nil) != c.wantErr {
			t.Errorf("%s: validateSigmaRange(%g, %g, %g) error = %v, wantErr %v", c.name, c.min, c.max, c.step, err, c.wantErr)
		}
	}
}
