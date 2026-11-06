package goldencmp

import (
	"fmt"
	"math"
	"strconv"
	"strings"
)

type Tolerances struct {
	FloatRel     float64
	MinMagnitude float64
}

type Mismatch struct {
	Line               int
	Field, Golden, Got string
	RelDiff            float64
}

// String formats the mismatch with its line, field, both values and relative difference.
func (m Mismatch) String() string {
	return fmt.Sprintf("line %d, field %s: golden %q, got %q, rel diff %.3g", m.Line, m.Field, m.Golden, m.Got, m.RelDiff)
}

var tableColumns = []string{"Position", "Base", "P(x)", "G(x)"}

// Compare returns every difference between golden and got: labels and integers must match exactly, floats within tol.FloatRel, and a line-count difference is reported alongside the per-line mismatches.
func Compare(golden, got []byte, tol Tolerances) []Mismatch {
	g := strings.Split(string(golden), "\n")
	r := strings.Split(string(got), "\n")
	var out []Mismatch
	if len(g) != len(r) {
		out = append(out, Mismatch{
			Line:    min(len(g), len(r)) + 1,
			Field:   "line count",
			Golden:  strconv.Itoa(len(g)),
			Got:     strconv.Itoa(len(r)),
			RelDiff: math.Inf(1),
		})
	}
	for i := 0; i < min(len(g), len(r)); i++ {
		if g[i] != r[i] {
			out = append(out, compareLine(i+1, g[i], r[i], tol)...)
		}
	}
	return out
}

// compareLine tokenizes a differing line by tab (per-base table rows when it has four columns) or by space (header text), applying the magnitude floor only to table rows.
func compareLine(line int, golden, got string, tol Tolerances) []Mismatch {
	if strings.Contains(golden, "\t") {
		isTable := strings.Count(golden, "\t") == len(tableColumns)-1
		return compareTokens(line, strings.Split(golden, "\t"), strings.Split(got, "\t"), func(_ []string, i int) string {
			if isTable {
				return tableColumns[i]
			}
			return fmt.Sprintf("column %d", i)
		}, tol, isTable)
	}
	return compareTokens(line, strings.Split(golden, " "), strings.Split(got, " "), textFieldName, tol, false)
}

// textFieldName names a numeric header token by the label words preceding it, so a mismatch reports its field rather than an index.
func textFieldName(tokens []string, i int) string {
	if !isNumber(tokens[i]) {
		return fmt.Sprintf("text token %d", i)
	}
	var words []string
	for j := i - 1; j >= 0 && !isNumber(tokens[j]); j-- {
		if strings.HasSuffix(tokens[j], ";") {
			break
		}
		if tokens[j] != "=" {
			words = append([]string{tokens[j]}, words...)
		}
	}
	if len(words) == 0 {
		return fmt.Sprintf("token %d", i)
	}
	return strings.Join(words, " ")
}

// compareTokens reports a field-count mismatch and each differing token pair, naming fields via name; floor enables the MinMagnitude rule.
func compareTokens(line int, golden, got []string, name func([]string, int) string, tol Tolerances, floor bool) []Mismatch {
	var out []Mismatch
	if len(golden) != len(got) {
		out = append(out, Mismatch{
			Line:    line,
			Field:   "field count",
			Golden:  fmt.Sprintf("%d fields", len(golden)),
			Got:     fmt.Sprintf("%d fields", len(got)),
			RelDiff: math.Inf(1),
		})
	}
	for i := 0; i < min(len(golden), len(got)); i++ {
		if m, ok := compareToken(golden[i], got[i], tol, floor); !ok {
			m.Line = line
			m.Field = name(golden, i)
			out = append(out, m)
		}
	}
	return out
}

// isNumber reports whether s parses as a float and contains a digit, excluding Inf and NaN spellings.
func isNumber(s string) bool {
	if !strings.ContainsAny(s, "0123456789") {
		return false
	}
	_, err := strconv.ParseFloat(s, 64)
	return err == nil
}

// isInteger reports whether s parses as a base-10 integer.
func isInteger(s string) bool {
	_, err := strconv.ParseInt(s, 10, 64)
	return err == nil
}

// relDiff returns |got-golden|/|golden|, or +Inf if either side is not a number.
func relDiff(golden, got string) float64 {
	if !isNumber(golden) || !isNumber(got) {
		return math.Inf(1)
	}
	gv, _ := strconv.ParseFloat(golden, 64)
	rv, _ := strconv.ParseFloat(got, 64)
	return math.Abs(rv-gv) / math.Abs(gv)
}

// compareToken passes identical tokens; otherwise only two non-integer floats pass, within tol.FloatRel, and when floor is set values at or below MinMagnitude pass only if both are, so an integer/float or label change never passes silently.
func compareToken(golden, got string, tol Tolerances, floor bool) (Mismatch, bool) {
	m := Mismatch{Golden: golden, Got: got, RelDiff: relDiff(golden, got)}
	if golden == got {
		return m, true
	}
	if !isNumber(golden) || isInteger(golden) || !isNumber(got) || isInteger(got) {
		return m, false
	}
	gv, _ := strconv.ParseFloat(golden, 64)
	rv, _ := strconv.ParseFloat(got, 64)
	if floor && (math.Abs(gv) <= tol.MinMagnitude || math.Abs(rv) <= tol.MinMagnitude) {
		return m, math.Abs(gv) <= tol.MinMagnitude && math.Abs(rv) <= tol.MinMagnitude
	}
	return m, m.RelDiff <= tol.FloatRel
}
