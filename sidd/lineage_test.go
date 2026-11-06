package sidd

import (
	"bufio"
	"crypto/sha256"
	"encoding/hex"
	"io/fs"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

const repoRoot = ".."

// readManifest parses a sha256sum-format manifest into a path-to-digest map, failing on malformed, duplicate or empty input.
func readManifest(t *testing.T, path string) map[string]string {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	entries := map[string]string{}
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		sum, rel, ok := strings.Cut(scanner.Text(), "  ")
		if !ok || len(sum) != sha256.Size*2 || rel == "" {
			t.Fatalf("%s: malformed line %q", path, scanner.Text())
		}
		if _, dup := entries[rel]; dup {
			t.Fatalf("%s: %s listed twice", path, rel)
		}
		entries[rel] = sum
	}
	if err := scanner.Err(); err != nil {
		t.Fatal(err)
	}
	if len(entries) == 0 {
		t.Fatalf("%s is empty", path)
	}
	return entries
}

// TestGoldenLineage pins every file listed in the manifest to its sha256 entry and requires only the stdout.txt files under the two sidd golden dirs to be listed; the six v1 rlooper/sim golden dirs are not covered. Each manifest must also list its SIST patches.
func TestGoldenLineage(t *testing.T) {
	for _, dir := range []string{"testdata/golden/sidd", "testdata/golden/sidd-multirun"} {
		t.Run(dir, func(t *testing.T) {
			entries := readManifest(t, filepath.Join(repoRoot, dir, "MANIFEST.sha256"))
			patches := []string{"scripts/sist-1r.patch"}
			if dir == "testdata/golden/sidd-multirun" {
				patches = []string{"scripts/sist-dedup-fix.patch", "scripts/sist-probe.patch"}
			}
			for _, patch := range patches {
				if _, ok := entries[patch]; !ok {
					t.Errorf("%s is not listed in the manifest", patch)
				}
			}
			for rel, want := range entries {
				data, err := os.ReadFile(filepath.Join(repoRoot, rel))
				if err != nil {
					t.Errorf("%v", err)
					continue
				}
				sum := sha256.Sum256(data)
				if got := hex.EncodeToString(sum[:]); got != want {
					t.Errorf("%s: sha256 %s, manifest has %s", rel, got, want)
				}
			}
			err := filepath.WalkDir(filepath.Join(repoRoot, dir), func(path string, d fs.DirEntry, err error) error {
				if err != nil || d.IsDir() || d.Name() != "stdout.txt" {
					return err
				}
				rel, err := filepath.Rel(repoRoot, path)
				if err != nil {
					return err
				}
				if _, ok := entries[filepath.ToSlash(rel)]; !ok {
					t.Errorf("%s is not listed in the manifest", rel)
				}
				return nil
			})
			if err != nil {
				t.Fatal(err)
			}
		})
	}
}
