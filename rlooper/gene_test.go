package rlooper

import (
	"os"
	"strings"
	"testing"
)

// TestReadGenesFixture pins header parsing of name, chromosome, strand and range plus sequence length for synthetic200.
func TestReadGenesFixture(t *testing.T) {
	f, err := os.Open("../testdata/fixtures/synthetic200.fa")
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	genes, err := ReadGenes(f)
	if err != nil {
		t.Fatal(err)
	}
	g := genes[0]
	if g.Name != "SYNTH200" || g.Loc.Chromosome != "chrT" || g.Loc.Strand != "+" ||
		g.Loc.Start != 1000 || g.Loc.End != 1199 || len(g.Sequence) != 200 {
		t.Fatalf("%+v len=%d", g.Loc, len(g.Sequence))
	}
}

// TestReadGenesMultiRecordLowercaseStrand pins multi-record parsing and uppercasing of lowercase bases.
func TestReadGenesMultiRecordLowercaseStrand(t *testing.T) {
	in := ">a range=chr1:10-13 strand=-\nacgt\n>b range=chr2:5-8 strand=+\nGGCC\n"
	genes, err := ReadGenes(strings.NewReader(in))
	if err != nil || len(genes) != 2 {
		t.Fatalf("genes=%d err=%v", len(genes), err)
	}
	if string(genes[0].Sequence) != "ACGT" || genes[0].Loc.Strand != "-" {
		t.Fatalf("%s %s", genes[0].Sequence, genes[0].Loc.Strand)
	}
}

// TestReadGenesRejectsInvalidBase pins that a non-ACGT base is an error naming the base, record and 1-based position.
func TestReadGenesRejectsInvalidBase(t *testing.T) {
	_, err := ReadGenes(strings.NewReader(">x range=c:1-4 strand=+\nACNT\n"))
	if err == nil || !strings.Contains(err.Error(), "N") {
		t.Fatalf("want invalid-base error, got %v", err)
	}
	if !strings.Contains(err.Error(), "x") || !strings.Contains(err.Error(), "3") {
		t.Fatalf("want record name and position in error, got %v", err)
	}
}

// TestReadGenesRejectsInvalidBaseReportsOriginalByte pins that the error reports the original lowercase byte, not its uppercased form.
func TestReadGenesRejectsInvalidBaseReportsOriginalByte(t *testing.T) {
	_, err := ReadGenes(strings.NewReader(">x range=c:1-4 strand=+\nACnT\n"))
	if err == nil || !strings.Contains(err.Error(), "n") {
		t.Fatalf("want invalid-base error reporting original lowercase byte, got %v", err)
	}
	if strings.Contains(err.Error(), "N") {
		t.Fatalf("want original lowercase byte, not uppercased, got %v", err)
	}
	if !strings.Contains(err.Error(), "x") || !strings.Contains(err.Error(), "3") {
		t.Fatalf("want record name and position in error, got %v", err)
	}
}

// TestReadGenesLongLine pins that a single 120000-base line is read without a line-length limit.
func TestReadGenesLongLine(t *testing.T) {
	long := strings.Repeat("ACGT", 30000)
	genes, err := ReadGenes(strings.NewReader(">x range=c:1-120000 strand=+\n" + long + "\n"))
	if err != nil || len(genes[0].Sequence) != 120000 {
		t.Fatalf("len=%d err=%v", len(genes[0].Sequence), err)
	}
}

// TestOrientation pins strand-based orientation (+ complements, - reverses) and the + strand with invert case.
func TestOrientation(t *testing.T) {
	g := &Gene{Sequence: []byte("AACG"), Loc: Loci{Strand: "+"}}
	g.ApplyOrientation(false, false, false)
	if string(g.Sequence) != "TTGC" {
		t.Fatal(string(g.Sequence))
	}
	g = &Gene{Sequence: []byte("AACG"), Loc: Loci{Strand: "-"}}
	g.ApplyOrientation(false, false, false)
	if string(g.Sequence) != "GCAA" {
		t.Fatal(string(g.Sequence))
	}
	g = &Gene{Sequence: []byte("AACG"), Loc: Loci{Strand: "+"}}
	g.ApplyOrientation(false, false, true)
	if string(g.Sequence) != "GCAA" {
		t.Fatal(string(g.Sequence))
	}
}

// TestApplyOrientationInvertFlagCombinesComplementAndReverse pins that invert alone reverse-complements a strandless gene.
func TestApplyOrientationInvertFlagCombinesComplementAndReverse(t *testing.T) {
	g := &Gene{Sequence: []byte("AACG"), Loc: Loci{Strand: ""}}
	g.ApplyOrientation(false, false, true)
	if string(g.Sequence) != "CGTT" {
		t.Fatal(string(g.Sequence))
	}
}

// TestApplyOrientationComplementFlagIdempotentWithInvert pins the C++ flag interaction where complement plus invert on a + strand complements twice and reverses once.
func TestApplyOrientationComplementFlagIdempotentWithInvert(t *testing.T) {
	g := &Gene{Sequence: []byte("AACG"), Loc: Loci{Strand: "+"}}
	g.ApplyOrientation(true, false, true)
	want := "GCAA"
	if string(g.Sequence) != want {
		t.Fatalf("ApplyOrientation(complement=true, reverse=false, invert=true) on strand \"+\" AACG = %q, want %q (traced: strand \"+\" complements AACG->TTGC; complement_flag, set true by both --complement and --invert, complements once more TTGC->AACG; reverse_flag, set true by --invert, reverses once AACG->GCAA)", string(g.Sequence), want)
	}
}

// TestReadGenesHeaderless pins that a header with no range or strand attributes yields zero-valued Loci.
func TestReadGenesHeaderless(t *testing.T) {
	genes, err := ReadGenes(strings.NewReader(">x\nACGT\n"))
	if err != nil || len(genes) != 1 {
		t.Fatalf("genes=%d err=%v", len(genes), err)
	}
	g := genes[0]
	if g.Name != "x" || g.Loc.Chromosome != "" || g.Loc.Strand != "" || g.Loc.Start != 0 || g.Loc.End != 0 || string(g.Sequence) != "ACGT" {
		t.Fatalf("%+v %s", g.Loc, g.Sequence)
	}
}
