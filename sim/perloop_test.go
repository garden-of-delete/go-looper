package sim

import (
	"bufio"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"golooper/config"
)

func parsePerloopRows(t *testing.T, path string) (sigmas, probs []float64) {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()

	sc := bufio.NewScanner(f)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) != 2 {
			t.Fatalf("%s: malformed row %q", path, line)
		}
		sigma, err := strconv.ParseFloat(fields[0], 64)
		if err != nil {
			t.Fatal(err)
		}
		prob, err := strconv.ParseFloat(fields[1], 64)
		if err != nil {
			t.Fatal(err)
		}
		sigmas = append(sigmas, sigma)
		probs = append(probs, prob)
	}
	if err := sc.Err(); err != nil {
		t.Fatal(err)
	}
	return sigmas, probs
}

func checkSigma(t *testing.T, got, want float64) {
	t.Helper()
	if diff := got - want; diff < -1e-6 || diff > 1e-6 {
		t.Fatalf("sigma got %g, want %g", got, want)
	}
}

func TestSigmaGrid(t *testing.T) {
	grid := sigmaGrid(-0.10, -0.02, 0.01)
	if len(grid) != 9 {
		t.Fatalf("len(grid) = %d, want 9", len(grid))
	}
	checkSigma(t, grid[len(grid)-1], -0.02)

	grid = sigmaGrid(-0.10, -0.025, 0.01)
	if len(grid) != 8 {
		t.Fatalf("len(grid) = %d, want 8", len(grid))
	}
	checkSigma(t, grid[len(grid)-1], -0.03)
}

func TestPerloopHonorsAutoN(t *testing.T) {
	dir := t.TempDir()

	cfgDefault := &config.Config{
		InfileName:  "../testdata/fixtures/synthetic200.fa",
		OutfileName: filepath.Join(dir, "default"),
	}
	if err := RunPerloop(cfgDefault, -0.10, -0.02, 0.01, 2); err != nil {
		t.Fatalf("RunPerloop (default N): %v", err)
	}
	_, defaultProbs := parsePerloopRows(t, filepath.Join(dir, "default_perloop.txt"))

	cfgAuto := &config.Config{
		InfileName:     "../testdata/fixtures/synthetic200.fa",
		OutfileName:    filepath.Join(dir, "auto"),
		AutoDomainSize: true,
	}
	if err := RunPerloop(cfgAuto, -0.10, -0.02, 0.01, 2); err != nil {
		t.Fatalf("RunPerloop (auto N): %v", err)
	}
	_, autoProbs := parsePerloopRows(t, filepath.Join(dir, "auto_perloop.txt"))

	if len(defaultProbs) != len(autoProbs) {
		t.Fatalf("row count mismatch: %d vs %d", len(defaultProbs), len(autoProbs))
	}

	const tol = 1e-6
	differs := false
	for i := range defaultProbs {
		if diff := defaultProbs[i] - autoProbs[i]; diff > tol || diff < -tol {
			differs = true
			break
		}
	}
	if !differs {
		t.Fatalf("expected --N auto to change perloop output under synthetic200, got identical probabilities: %v", autoProbs)
	}
}

func TestPerloopCircularDiffersFromLinear(t *testing.T) {
	dir := t.TempDir()
	pinnedN := 1500

	linCfg := &config.Config{
		InfileName:          "../testdata/fixtures/synthetic200.fa",
		OutfileName:         filepath.Join(dir, "linear"),
		SuperhelicityDomain: &pinnedN,
	}
	if err := RunPerloop(linCfg, -0.10, -0.02, 0.01, 2); err != nil {
		t.Fatalf("RunPerloop (linear): %v", err)
	}
	_, linearProbs := parsePerloopRows(t, filepath.Join(dir, "linear_perloop.txt"))

	circCfg := &config.Config{
		InfileName:          "../testdata/fixtures/synthetic200.fa",
		OutfileName:         filepath.Join(dir, "circular"),
		SuperhelicityDomain: &pinnedN,
		Circular:            true,
	}
	if err := RunPerloop(circCfg, -0.10, -0.02, 0.01, 2); err != nil {
		t.Fatalf("RunPerloop (circular): %v", err)
	}
	_, circularProbs := parsePerloopRows(t, filepath.Join(dir, "circular_perloop.txt"))

	if len(linearProbs) != len(circularProbs) {
		t.Fatalf("row count mismatch: %d vs %d", len(linearProbs), len(circularProbs))
	}

	differs := false
	for i := range linearProbs {
		tol := math.Max(1e-9, 1e-5*math.Abs(linearProbs[i]))
		if diff := circularProbs[i] - linearProbs[i]; diff > tol || diff < -tol {
			differs = true
			break
		}
	}
	if !differs {
		t.Fatalf("expected --circular to change perloop output on synthetic200, got identical probabilities: linear=%v circular=%v", linearProbs, circularProbs)
	}
}

func TestRunPerloopAgainstGolden(t *testing.T) {
	dir := t.TempDir()
	cfg := &config.Config{
		InfileName:  "../testdata/fixtures/synthetic200.fa",
		OutfileName: filepath.Join(dir, "out"),
	}

	if err := RunPerloop(cfg, -0.10, -0.02, 0.01, 2); err != nil {
		t.Fatalf("RunPerloop: %v", err)
	}

	gotSigmas, gotProbs := parsePerloopRows(t, filepath.Join(dir, "out_perloop.txt"))
	wantSigmas, wantProbs := parsePerloopRows(t, "../testdata/golden/perloop/perloop_out.txt")

	if len(gotSigmas) != len(wantSigmas) {
		t.Fatalf("row count got %d, want %d", len(gotSigmas), len(wantSigmas))
	}

	for i := range gotSigmas {
		checkSigma(t, gotSigmas[i], wantSigmas[i])
		checkTol(t, gotProbs[i], wantProbs[i])
	}
}
