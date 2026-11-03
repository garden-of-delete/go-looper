package sim

import (
	"bufio"
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
