package goldencmp

import (
	"bytes"
	"math"
	"os"
	"strings"
	"testing"
)

const goldenPath = "../../testdata/golden/sidd-multirun/mr-clamp-copoly-defaults/stdout.txt"

var canonical = Tolerances{FloatRel: 1e-6, MinMagnitude: 1e-12}

// loadGolden reads the committed multi-run golden used as the comparator's fixture.
func loadGolden(t *testing.T) []byte {
	t.Helper()
	data, err := os.ReadFile(goldenPath)
	if err != nil {
		t.Fatal(err)
	}
	return data
}

// replaceLine returns data with line index rewritten by edit.
func replaceLine(t *testing.T, data []byte, index int, edit func(string) string) []byte {
	t.Helper()
	lines := strings.Split(string(data), "\n")
	lines[index] = edit(lines[index])
	return []byte(strings.Join(lines, "\n"))
}

// replaceOnce replaces old in line, failing the test unless old occurs exactly once.
func replaceOnce(t *testing.T, line, old, replacement string) string {
	t.Helper()
	if strings.Count(line, old) != 1 {
		t.Fatalf("%q does not occur exactly once in %q", old, line)
	}
	return strings.Replace(line, old, replacement, 1)
}

// single returns the only mismatch in mm, failing the test if there is not exactly one.
func single(t *testing.T, mm []Mismatch) Mismatch {
	t.Helper()
	if len(mm) != 1 {
		t.Fatalf("got %d mismatches, want 1: %v", len(mm), mm)
	}
	return mm[0]
}

const (
	lineMinMax     = 13
	lineTwoRun     = 15
	lineThreeRun   = 16
	lineScaling    = 20
	lineMelted     = 23
	lineFirstTable = 25
)

// TestCompareSelf pins that a golden compared against itself yields no mismatches.
func TestCompareSelf(t *testing.T) {
	golden := loadGolden(t)
	if got := Compare(golden, golden, canonical); len(got) != 0 {
		t.Fatalf("self comparison mismatched: %v", got)
	}
}

// TestCompareDoctored pins that doctored goldens are caught: tolerated 7th-digit noise passes, while real value, type, label, whitespace, field-count and line-count changes fail with the right line and field.
func TestCompareDoctored(t *testing.T) {
	golden := loadGolden(t)

	t.Run("seventh digit of P passes", func(t *testing.T) {
		got := replaceLine(t, golden, lineFirstTable-1, func(s string) string {
			return replaceOnce(t, s, "e-07", "1e-07")
		})
		if mm := Compare(golden, got, canonical); len(mm) != 0 {
			t.Fatalf("unexpected mismatches: %v", mm)
		}
	})

	t.Run("fifth digit of P fails with line, field and rel diff", func(t *testing.T) {
		got := replaceLine(t, golden, lineFirstTable-1, func(s string) string {
			return replaceOnce(t, s, "2.43341e-07", "2.43391e-07")
		})
		m := single(t, Compare(golden, got, canonical))
		want := math.Abs(2.43391e-07-2.43341e-07) / 2.43341e-07
		if m.Line != lineFirstTable || m.Field != "P(x)" || m.Golden != "2.43341e-07" || m.Got != "2.43391e-07" || math.Abs(m.RelDiff-want) > 1e-12 {
			t.Fatalf("unexpected mismatch %+v, want rel diff %v", m, want)
		}
	})

	t.Run("two-run count off by one fails exactly with rel diff", func(t *testing.T) {
		got := replaceLine(t, golden, lineTwoRun-1, func(s string) string {
			return replaceOnce(t, s, "4854993", "4854994")
		})
		m := single(t, Compare(golden, got, canonical))
		want := 1.0 / 4854993
		if m.Line != lineTwoRun || m.Field != "Number of two-run states" || math.Abs(m.RelDiff-want) > 1e-12 {
			t.Fatalf("unexpected mismatch %+v", m)
		}
	})

	t.Run("integer printed as double fails", func(t *testing.T) {
		got := replaceLine(t, golden, lineThreeRun-1, func(s string) string {
			return replaceOnce(t, s, "2.76978e+06", "2769779")
		})
		m := single(t, Compare(golden, got, canonical))
		if m.Line != lineThreeRun || m.Field != "Number of three-run states" || m.Golden != "2.76978e+06" || m.Got != "2769779" {
			t.Fatalf("unexpected mismatch %+v", m)
		}
	})

	t.Run("double printed as integer fails", func(t *testing.T) {
		got := replaceLine(t, golden, lineTwoRun-1, func(s string) string {
			return replaceOnce(t, s, "4854993", "4.85499e+06")
		})
		m := single(t, Compare(golden, got, canonical))
		if m.Line != lineTwoRun {
			t.Fatalf("unexpected mismatch %+v", m)
		}
	})

	t.Run("label text mismatch is named", func(t *testing.T) {
		got := replaceLine(t, golden, lineTwoRun-1, func(s string) string {
			return replaceOnce(t, s, "two-run", "2-run")
		})
		m := single(t, Compare(golden, got, canonical))
		if m.Line != lineTwoRun || m.Field != "text token 2" || m.Golden != "two-run" || m.Got != "2-run" {
			t.Fatalf("unexpected mismatch %+v", m)
		}
	})

	t.Run("inserted equals sign changes field count", func(t *testing.T) {
		got := replaceLine(t, golden, lineScaling-1, func(s string) string {
			return replaceOnce(t, s, "factor 300", "factor = 300")
		})
		mm := Compare(golden, got, canonical)
		if len(mm) == 0 || mm[0].Line != lineScaling || mm[0].Field != "field count" {
			t.Fatalf("unexpected mismatches: %v", mm)
		}
	})

	t.Run("extra internal space is detected", func(t *testing.T) {
		got := replaceLine(t, golden, 8, func(s string) string {
			return replaceOnce(t, s, "= 0.01 M", "=  0.01 M")
		})
		if mm := Compare(golden, got, canonical); len(mm) == 0 {
			t.Fatal("extra internal space passed")
		}
	})

	t.Run("trailing space is detected", func(t *testing.T) {
		got := replaceLine(t, golden, 7, func(s string) string { return s + " " })
		if mm := Compare(golden, got, canonical); len(mm) == 0 {
			t.Fatal("trailing space passed")
		}
	})

	t.Run("trailing carriage return is detected", func(t *testing.T) {
		got := replaceLine(t, golden, lineScaling-1, func(s string) string { return s + "\r" })
		m := single(t, Compare(golden, got, canonical))
		if m.Line != lineScaling || m.Got != "300\r" {
			t.Fatalf("unexpected mismatch %+v", m)
		}
	})

	t.Run("table row field count", func(t *testing.T) {
		got := replaceLine(t, golden, lineFirstTable-1, func(s string) string { return s + "\t1" })
		mm := Compare(golden, got, canonical)
		if len(mm) != 1 || mm[0].Line != lineFirstTable || mm[0].Field != "field count" {
			t.Fatalf("unexpected mismatches: %v", mm)
		}
	})

	t.Run("fields on the MinE line are named by their own label", func(t *testing.T) {
		got := replaceLine(t, golden, lineMinMax-1, func(s string) string {
			return replaceOnce(t, s, "99.6981", "99.7981")
		})
		m := single(t, Compare(golden, got, canonical))
		if m.Line != lineMinMax || m.Field != "Max_E" {
			t.Fatalf("unexpected mismatch %+v", m)
		}
	})

	t.Run("BGC is not named BAT", func(t *testing.T) {
		got := replaceLine(t, golden, lineMinMax-2, func(s string) string {
			return replaceOnce(t, s, "1.30488", "1.31488")
		})
		m := single(t, Compare(golden, got, canonical))
		if m.Field != "BGC" {
			t.Fatalf("unexpected mismatch %+v", m)
		}
	})

	t.Run("base letter is exact", func(t *testing.T) {
		got := replaceLine(t, golden, lineFirstTable-1, func(s string) string {
			parts := strings.Split(s, "\t")
			parts[1] = "N"
			return strings.Join(parts, "\t")
		})
		m := single(t, Compare(golden, got, canonical))
		if m.Field != "Base" {
			t.Fatalf("unexpected mismatch %+v", m)
		}
	})

	t.Run("deleted last line reports a line count mismatch", func(t *testing.T) {
		trimmed := bytes.TrimSuffix(golden, []byte("\n"))
		cut := trimmed[:bytes.LastIndexByte(trimmed, '\n')+1]
		mm := Compare(golden, cut, canonical)
		if len(mm) == 0 || mm[0].Field != "line count" || mm[0].Golden != "4386" || mm[0].Got != "4385" {
			t.Fatalf("unexpected mismatches: %v", mm)
		}
		last := mm[len(mm)-1]
		if last.Line != 4385 || last.Field != "Position" || last.Golden != "4361" || last.Got != "" {
			t.Fatalf("the missing row is not named: %v", mm)
		}
	})
}

// TestCompareMinMagnitude pins that the magnitude floor lets two sub-floor table values pass but not a sub-floor against an above-floor value.
func TestCompareMinMagnitude(t *testing.T) {
	golden := []byte("12\tA\t1e-13\t5.5\n")
	if mm := Compare(golden, []byte("12\tA\t2e-13\t5.5\n"), canonical); len(mm) != 0 {
		t.Fatalf("table values both below floor must pass: %v", mm)
	}
	if mm := Compare(golden, []byte("12\tA\t1e-9\t5.5\n"), canonical); len(mm) != 1 {
		t.Fatalf("table value above floor must fail: %v", mm)
	}
}

// TestCompareHeaderFloatsIgnoreFloor pins that a header float such as the partition function is always compared, even when far below the floor.
func TestCompareHeaderFloatsIgnoreFloor(t *testing.T) {
	golden := loadGolden(t)
	const prefix = "Total Partition Function = "
	idx := -1
	for i, l := range strings.Split(string(golden), "\n") {
		if strings.HasPrefix(l, prefix) {
			idx = i
		}
	}
	if idx < 0 {
		t.Fatal("no Total Partition Function line")
	}
	for name, value := range map[string]string{
		"scaled by 7.4":      "1.62135e-59",
		"below floor":        "9.99999e-13",
		"order of magnitude": "2.19101e-59",
	} {
		t.Run(name, func(t *testing.T) {
			got := replaceLine(t, golden, idx, func(string) string { return prefix + value })
			m := single(t, Compare(golden, got, canonical))
			if m.Line != idx+1 || m.Field != "Total Partition Function" {
				t.Fatalf("unexpected mismatch %+v", m)
			}
		})
	}
}
