package main

import (
	"bytes"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

const (
	fixturesDir = "../../testdata/fixtures"
	goldenDir   = "../../testdata/golden/sidd"
)

func chdir(t *testing.T, dir string) {
	t.Helper()
	prev, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	if err := os.Chdir(dir); err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() {
		if err := os.Chdir(prev); err != nil {
			t.Fatal(err)
		}
	})
}

func runIn(t *testing.T, dir string, args ...string) (string, error) {
	t.Helper()
	chdir(t, dir)
	var out bytes.Buffer
	cmd := newRootCmd()
	cmd.SetOut(&out)
	cmd.SetErr(&bytes.Buffer{})
	cmd.SetArgs(args)
	err := cmd.Execute()
	return out.String(), err
}

func stageContent(t *testing.T, name string, content []byte) string {
	t.Helper()
	dir := t.TempDir()
	if err := os.WriteFile(filepath.Join(dir, name), content, 0o644); err != nil {
		t.Fatal(err)
	}
	return dir
}

func stageFixture(t *testing.T, name string) string {
	t.Helper()
	data, err := os.ReadFile(filepath.Join(fixturesDir, name))
	if err != nil {
		t.Fatal(err)
	}
	return stageContent(t, name, data)
}

func TestFlagValidation(t *testing.T) {
	short := ">s\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
	tests := []struct {
		name    string
		content string
		args    []string
		wantErr string
	}{
		{"clamp and circular", short, []string{"--clamp", "--circular"}, "clamp"},
		{"linear nearest-neighbor", short, []string{"--nearest-neighbor"}, "clamp or circular"},
		{"multi-record", short + ">t\nACGT\n", []string{"--clamp"}, "exactly one"},
		{"ambiguous base", ">s\nACGTNACGT\n", []string{"--clamp"}, "position 5"},
		{"positive sigma", short, []string{"--clamp", "--sigma", "0.06"}, "negative --sigma"},
		{"bad N", short, []string{"--clamp", "--N", "-3"}, "--N"},
		{"leading blank line fasta", "\n" + short, []string{"--clamp", "--maxlength", "20"}, ""},
		{"comment-led fasta", ";note\n" + short, []string{"--clamp", "--maxlength", "20"}, ""},
		{"lowercase accepted", ">s\nacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt\n", []string{"--clamp", "--maxlength", "20"}, ""},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			dir := stageContent(t, "in.fa", []byte(tt.content))
			_, err := runIn(t, dir, append([]string{"-f", "in.fa"}, tt.args...)...)
			if tt.wantErr == "" {
				if err != nil {
					t.Fatalf("unexpected error: %v", err)
				}
				return
			}
			if err == nil || !strings.Contains(err.Error(), tt.wantErr) {
				t.Fatalf("error = %v, want containing %q", err, tt.wantErr)
			}
		})
	}
}

func TestEndToEndParityMatrix(t *testing.T) {
	cases := []struct {
		name, fixture string
		flags         []string
	}{
		{"clamp-copoly-defaults", "pbr322.toy.fa", []string{"--clamp"}},
		{"clamp-nn-defaults", "pbr322.toy.fa", []string{"--clamp", "--nearest-neighbor"}},
		{"circular-copoly-defaults", "pbr322.toy.fa", []string{"--circular"}},
		{"circular-nn-defaults", "pbr322.toy.fa", []string{"--circular", "--nearest-neighbor"}},
		{"clamp-copoly-alt", "pbr322.toy.fa", []string{"--clamp", "--temperature", "315", "--salt", "0.1", "--sigma", "-0.04", "--threshold", "10"}},
		{"clamp-copoly-synth", "sidd_synth1600.fa", []string{"--clamp"}},
		{"circular-nn-synth", "sidd_synth1600.fa", []string{"--circular", "--nearest-neighbor"}},
		{"circular-nn-alt", "sidd_synth1600.fa", []string{"--circular", "--nearest-neighbor", "--temperature", "300", "--salt", "0.05", "--sigma", "-0.07", "--threshold", "11"}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			want, err := os.ReadFile(filepath.Join(goldenDir, tc.name, "stdout.txt"))
			if err != nil {
				t.Fatal(err)
			}
			dir := stageFixture(t, tc.fixture)
			got, err := runIn(t, dir, append([]string{"-f", tc.fixture}, tc.flags...)...)
			if err != nil {
				t.Fatal(err)
			}
			if got != string(want) {
				t.Fatalf("output differs from golden\n%s", firstDiff(string(want), got))
			}
		})
	}
}

func firstDiff(want, got string) string {
	w, g := strings.Split(want, "\n"), strings.Split(got, "\n")
	for i := 0; i < len(w) || i < len(g); i++ {
		var wl, gl string
		if i < len(w) {
			wl = w[i]
		}
		if i < len(g) {
			gl = g[i]
		}
		if wl != gl {
			return "line " + strconv.Itoa(i+1) + "\nwant: " + wl + "\ngot:  " + gl
		}
	}
	return "no differing line"
}

func TestTracks(t *testing.T) {
	raw, err := os.ReadFile(filepath.Join(fixturesDir, "pbr322.toy.fa"))
	if err != nil {
		t.Fatal(err)
	}
	seq := strings.Join(strings.Fields(string(raw)), "")
	n := strconv.Itoa(len(seq))
	tests := []struct {
		name, file, content   string
		wantChrom, wantBrowse string
		wantStart             int
	}{
		{"fixture", "pbr322.toy.fa", string(raw), "pbr322.toy.fa", "1-" + n, 1},
		{"headerless raw sequence", "my seq.txt", seq + "\n", "my_seq.txt", "1-" + n, 1},
		{"name header without range", "in.fa", ">pBR322 toy plasmid\n" + seq + "\n", "pBR322", "1-" + n, 1},
		{"range header untouched", "in.fa", ">pBR322 range=chrZ:101-" + strconv.Itoa(100+len(seq)) + "\n" + seq + "\n", "chrZ", "101-" + strconv.Itoa(100+len(seq)), 101},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			dir := stageContent(t, tt.file, []byte(tt.content))
			report, err := runIn(t, dir, "-f", tt.file, "--clamp", "--tracks", "trk")
			if err != nil {
				t.Fatal(err)
			}
			var wantP, wantG []string
			inTable := false
			for _, line := range strings.Split(report, "\n") {
				if strings.HasPrefix(line, "Position\t") {
					inTable = true
					continue
				}
				if !inTable || line == "" {
					continue
				}
				f := strings.Split(line, "\t")
				wantP = append(wantP, f[2])
				wantG = append(wantG, f[3])
			}
			if len(wantP) == 0 {
				t.Fatal("no table rows parsed")
			}
			for suffix, want := range map[string][]string{"_sidd_bpprob.wig": wantP, "_sidd_gx.wig": wantG} {
				data, err := os.ReadFile(filepath.Join(dir, "trk"+suffix))
				if err != nil {
					t.Fatal(err)
				}
				lines := strings.Split(strings.TrimRight(string(data), "\n"), "\n")
				if len(lines) < 3 {
					t.Fatalf("%s: too short", suffix)
				}
				wantBrowser := "browser position " + tt.wantChrom + ":" + tt.wantBrowse
				if lines[0] != wantBrowser {
					t.Errorf("%s: browser line = %q, want %q", suffix, lines[0], wantBrowser)
				}
				wantFixed := "fixedStep chrom=" + tt.wantChrom + " start=" + strconv.Itoa(tt.wantStart) + " step=1"
				if lines[2] != wantFixed {
					t.Errorf("%s: fixedStep line = %q, want %q", suffix, lines[2], wantFixed)
				}
				if strings.Join(lines[3:], "\n") != strings.Join(want, "\n") {
					t.Fatalf("%s values differ from report table", suffix)
				}
			}
		})
	}
}

func TestOutfile(t *testing.T) {
	dir := stageFixture(t, "pbr322.toy.fa")
	want, err := runIn(t, dir, "-f", "pbr322.toy.fa", "--clamp")
	if err != nil {
		t.Fatal(err)
	}
	stdout, err := runIn(t, dir, "-f", "pbr322.toy.fa", "--clamp", "-o", "report.txt")
	if err != nil {
		t.Fatal(err)
	}
	if stdout != "" {
		t.Errorf("stdout = %q, want empty with -o", stdout)
	}
	got, err := os.ReadFile(filepath.Join(dir, "report.txt"))
	if err != nil {
		t.Fatal(err)
	}
	if string(got) != want {
		t.Fatalf("outfile differs from stdout report\n%s", firstDiff(want, string(got)))
	}
}

func TestVersion(t *testing.T) {
	out, err := runIn(t, t.TempDir(), "--version")
	if err != nil {
		t.Fatal(err)
	}
	if out != "golooper-sidd v0.4.0\n" {
		t.Fatalf("version output = %q", out)
	}
}
