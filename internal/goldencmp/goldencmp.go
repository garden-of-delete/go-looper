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

func (m Mismatch) String() string {
	return fmt.Sprintf("line %d, field %s: golden %q, got %q, rel diff %.3g", m.Line, m.Field, m.Golden, m.Got, m.RelDiff)
}

var tableColumns = []string{"Position", "Base", "P(x)", "G(x)"}

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

func compareLine(line int, golden, got string, tol Tolerances) []Mismatch {
	if strings.Contains(golden, "\t") {
		isTable := strings.Count(golden, "\t") == len(tableColumns)-1
		return compareTokens(line, strings.Split(golden, "\t"), strings.Split(got, "\t"), func(_ []string, i int) string {
			if isTable {
				return tableColumns[i]
			}
			return fmt.Sprintf("column %d", i)
		}, tol)
	}
	return compareTokens(line, strings.Split(golden, " "), strings.Split(got, " "), textFieldName, tol)
}

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

func compareTokens(line int, golden, got []string, name func([]string, int) string, tol Tolerances) []Mismatch {
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
		if m, ok := compareToken(golden[i], got[i], tol); !ok {
			m.Line = line
			m.Field = name(golden, i)
			out = append(out, m)
		}
	}
	return out
}

func isNumber(s string) bool {
	if !strings.ContainsAny(s, "0123456789") {
		return false
	}
	_, err := strconv.ParseFloat(s, 64)
	return err == nil
}

func isInteger(s string) bool {
	_, err := strconv.ParseInt(s, 10, 64)
	return err == nil
}

func relDiff(golden, got string) float64 {
	if !isNumber(golden) || !isNumber(got) {
		return math.Inf(1)
	}
	gv, _ := strconv.ParseFloat(golden, 64)
	rv, _ := strconv.ParseFloat(got, 64)
	return math.Abs(rv-gv) / math.Abs(gv)
}

func compareToken(golden, got string, tol Tolerances) (Mismatch, bool) {
	m := Mismatch{Golden: golden, Got: got, RelDiff: relDiff(golden, got)}
	if golden == got {
		return m, true
	}
	if !isNumber(golden) || isInteger(golden) || !isNumber(got) || isInteger(got) {
		return m, false
	}
	gv, _ := strconv.ParseFloat(golden, 64)
	rv, _ := strconv.ParseFloat(got, 64)
	if math.Abs(gv) <= tol.MinMagnitude || math.Abs(rv) <= tol.MinMagnitude {
		return m, math.Abs(gv) <= tol.MinMagnitude && math.Abs(rv) <= tol.MinMagnitude
	}
	return m, m.RelDiff <= tol.FloatRel
}
