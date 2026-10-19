package rlooper

import (
	"bufio"
	"bytes"
	"math"
	"os"
	"strconv"
	"strings"
	"testing"
)

func TestEnsembleMatchesDirectSum(t *testing.T) {
	m := NewDefaultModel()
	seq := []byte("GATTACA")
	ens := ComputeEnsemble(seq, m, 2, 1)
	if ens.NumStructures != 21 {
		t.Fatalf("count=%d", ens.NumStructures)
	}
	p := PrefixBPEnergies(seq, m)
	var z float64
	for s := 0; s < len(seq); s++ {
		for e := s + 1; e < len(seq); e++ {
			g := m.SupercoilingEnergy(e-s+1) + p[e] - p[s]
			z += math.Exp(-g / m.RT())
		}
	}
	z += math.Exp(-m.GroundStateEnergy() / m.RT())
	if rel(ens.LogZ, math.Log(z)) > 1e-12 {
		t.Fatalf("logZ %g vs %g", ens.LogZ, math.Log(z))
	}
	recs := AllStructures(seq, m, 2, ens)
	var sum float64
	for _, r := range recs {
		sum += r.Prob
	}
	if math.Abs(sum+ens.PGround-1) > 1e-9 {
		t.Fatalf("probability sum %g", sum+ens.PGround)
	}
	if recs[0] != ens.MFE {
		t.Fatalf("MFE %+v, want recs[0] %+v", ens.MFE, recs[0])
	}
}

func TestEnsembleNoWindows(t *testing.T) {
	m := NewDefaultModel()
	ens := ComputeEnsemble([]byte("ACGT"), m, 10, 1)
	if ens.NumStructures != 0 || math.Abs(ens.PGround-1) > 1e-12 || math.IsNaN(ens.LogZ) {
		t.Fatalf("%+v", ens)
	}
	if ens.MFE != (StructureRecord{}) {
		t.Fatalf("MFE %+v, want zero value", ens.MFE)
	}
}

func TestEnsembleNoOverflow(t *testing.T) {
	m := NewDefaultModel()
	m.Homopolymer, m.OverrideEnergy = true, -5
	seq := bytes.Repeat([]byte("A"), 5000)
	ens := ComputeEnsemble(seq, m, 2, 1)
	if math.IsInf(ens.LogZ, 0) || math.IsNaN(ens.LogZ) || ens.PGround < 0 {
		t.Fatalf("%+v", ens)
	}
}

func rel(a, b float64) float64 { return math.Abs(a-b) / math.Max(math.Abs(b), 1e-300) }

func TestCompareStructures(t *testing.T) {
	tests := []struct {
		name string
		a, b StructureRecord
		want int
	}{
		{
			name: "lower energy sorts first regardless of position",
			a:    StructureRecord{Start: 10, End: 20, Energy: 1},
			b:    StructureRecord{Start: 0, End: 5, Energy: 2},
			want: -1,
		},
		{
			name: "equal energy, lower Start first",
			a:    StructureRecord{Start: 0, End: 5, Energy: 1},
			b:    StructureRecord{Start: 1, End: 5, Energy: 1},
			want: -1,
		},
		{
			name: "equal energy and Start, lower End first",
			a:    StructureRecord{Start: 0, End: 5, Energy: 1},
			b:    StructureRecord{Start: 0, End: 6, Energy: 1},
			want: -1,
		},
		{
			name: "identical records",
			a:    StructureRecord{Start: 0, End: 5, Energy: 1},
			b:    StructureRecord{Start: 0, End: 5, Energy: 1},
			want: 0,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := compareStructures(tt.a, tt.b); got != tt.want {
				t.Fatalf("compareStructures(%+v, %+v) = %d, want %d", tt.a, tt.b, got, tt.want)
			}
		})
	}
}

func TestEnsembleMFETieBreak(t *testing.T) {
	m := NewDefaultModel()
	m.Homopolymer, m.OverrideEnergy = true, 0.5
	seq := bytes.Repeat([]byte("A"), 150)
	ens := ComputeEnsemble(seq, m, 2, 1)
	recs := AllStructures(seq, m, 2, ens)

	if recs[1].Energy != recs[0].Energy {
		t.Fatalf("setup does not exercise a tie: recs[0].Energy=%g recs[1].Energy=%g", recs[0].Energy, recs[1].Energy)
	}
	if ens.MFE.Start != 0 {
		t.Fatalf("MFE.Start = %d, want 0 (first-in-enumeration should win the tie)", ens.MFE.Start)
	}
	if recs[0] != ens.MFE {
		t.Fatalf("MFE %+v, want recs[0] %+v", ens.MFE, recs[0])
	}
}

func TestLogSumExpMergeEmptyEmptyThenNonEmpty(t *testing.T) {
	a := newLogSumExp()
	b := newLogSumExp()
	a.merge(b)
	x := newLogSumExp()
	x.add(3.5)
	a.merge(x)
	if a.value() != 3.5 {
		t.Fatalf("value=%g want 3.5", a.value())
	}
}

func TestLogSumExpMergeEmptyTakesNonEmpty(t *testing.T) {
	empty := newLogSumExp()
	x := newLogSumExp()
	x.add(1)
	x.add(2)
	empty.merge(x)
	if empty.value() != x.value() {
		t.Fatalf("value=%g want %g", empty.value(), x.value())
	}
}

func TestLogSumExpMergeNonEmptyWithEmptyUnchanged(t *testing.T) {
	x := newLogSumExp()
	x.add(1)
	x.add(2)
	want := x.value()
	x.merge(newLogSumExp())
	if x.value() != want {
		t.Fatalf("value=%g want %g", x.value(), want)
	}
}

func TestTopStructuresMatchesAllStructures(t *testing.T) {
	genes := readFixture(t, "../testdata/fixtures/synthetic200.fa")
	g := genes[0]
	g.ApplyOrientation(false, false, false)
	m := NewDefaultModel()
	ens := ComputeEnsemble(g.Sequence, m, 2, 1)

	all := AllStructures(g.Sequence, m, 2, ens)
	top := TopStructures(g.Sequence, m, 2, ens, 5)

	if len(top) != 5 {
		t.Fatalf("len(top) = %d, want 5", len(top))
	}
	for i := range top {
		if top[i] != all[i] {
			t.Fatalf("TopStructures[%d] = %+v, want AllStructures[%d] = %+v", i, top[i], i, all[i])
		}
	}
}

func TestEnsembleAgainstGoldenDump(t *testing.T) {
	genes := readFixture(t, "../testdata/fixtures/synthetic200.fa")
	g := genes[0]
	g.ApplyOrientation(false, false, false)
	m := NewDefaultModel()
	ens := ComputeEnsemble(g.Sequence, m, 2, 1)
	recs := AllStructures(g.Sequence, m, 2, ens)
	want := parseDump(t, "../testdata/golden/defaults/outSYNTH200_dump.txt")
	if len(recs) != len(want) {
		t.Fatalf("count %d vs %d", len(recs), len(want))
	}
	for _, r := range recs {
		w, ok := want[[2]int64{int64(r.Start) + g.Loc.Start, int64(r.End) + g.Loc.Start}]
		if !ok {
			t.Fatalf("missing %d-%d", r.Start, r.End)
		}
		checkTol(t, r.Energy, w.energy)
		checkTol(t, r.Prob, w.prob)
	}
}

func readFixture(t *testing.T, path string) []*Gene {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	genes, err := ReadGenes(f)
	if err != nil {
		t.Fatal(err)
	}
	return genes
}

type goldenRow struct {
	energy, prob float64
}

func parseDump(t *testing.T, path string) map[[2]int64]goldenRow {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	want := make(map[[2]int64]goldenRow)
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64*1024), 1024*1024)
	first := true
	for sc.Scan() {
		if first {
			first = false
			continue
		}
		fields := strings.Fields(sc.Text())
		if len(fields) != 4 {
			t.Fatalf("malformed dump line: %q", sc.Text())
		}
		start, err := strconv.ParseInt(fields[0], 10, 64)
		if err != nil {
			t.Fatal(err)
		}
		stop, err := strconv.ParseInt(fields[1], 10, 64)
		if err != nil {
			t.Fatal(err)
		}
		if start == 0 && stop == 0 {
			continue
		}
		energy, err := strconv.ParseFloat(fields[2], 64)
		if err != nil {
			t.Fatal(err)
		}
		prob, err := strconv.ParseFloat(fields[3], 64)
		if err != nil {
			t.Fatal(err)
		}
		want[[2]int64{start, stop}] = goldenRow{energy: energy, prob: prob}
	}
	if err := sc.Err(); err != nil {
		t.Fatal(err)
	}
	return want
}

func checkTol(t *testing.T, got, want float64) {
	t.Helper()
	tol := math.Max(1e-9, 1e-5*math.Abs(want))
	if math.Abs(got-want) > tol {
		t.Fatalf("got %g want %g (tol %g)", got, want, tol)
	}
}
