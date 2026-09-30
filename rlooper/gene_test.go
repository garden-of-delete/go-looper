package rlooper

import (
	"os"
	"strings"
	"testing"
)

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

func TestReadGenesRejectsInvalidBase(t *testing.T) {
	_, err := ReadGenes(strings.NewReader(">x range=c:1-4 strand=+\nACNT\n"))
	if err == nil || !strings.Contains(err.Error(), "N") {
		t.Fatalf("want invalid-base error, got %v", err)
	}
}

func TestReadGenesLongLine(t *testing.T) {
	long := strings.Repeat("ACGT", 30000)
	genes, err := ReadGenes(strings.NewReader(">x range=c:1-120000 strand=+\n" + long + "\n"))
	if err != nil || len(genes[0].Sequence) != 120000 {
		t.Fatalf("len=%d err=%v", len(genes[0].Sequence), err)
	}
}

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

func TestApplyOrientationInvertFlagCombinesComplementAndReverse(t *testing.T) {
	g := &Gene{Sequence: []byte("AACG"), Loc: Loci{Strand: ""}}
	g.ApplyOrientation(false, false, true)
	if string(g.Sequence) != "CGTT" {
		t.Fatal(string(g.Sequence))
	}
}

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
