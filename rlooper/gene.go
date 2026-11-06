package rlooper

import (
	"bufio"
	"fmt"
	"io"
	"regexp"
	"strconv"
	"strings"
)

type Loci struct {
	Chromosome string
	Strand     string
	Start      int64
	End        int64
}

type Gene struct {
	Name     string
	Header   string
	Loc      Loci
	Sequence []byte
}

var (
	headerRangeRegex  = regexp.MustCompile(`range=([^:\s]+):(\d+)-(\d+)`)
	headerStrandRegex = regexp.MustCompile(`(?i)strand=([-+])`)
)

// ReadGenes parses FASTA records from r, uppercasing bases and skipping spaces and tabs, and returns an error naming the record and position of the first non-ACGT base.
func ReadGenes(r io.Reader) ([]*Gene, error) {
	br := bufio.NewReader(r)
	var genes []*Gene
	var current *Gene
	seqLen := 0

	for {
		line, readErr := br.ReadString('\n')
		line = strings.TrimRight(line, "\r\n")

		if strings.HasPrefix(line, ">") {
			current = parseHeader(line)
			genes = append(genes, current)
			seqLen = 0
		} else if current != nil {
			for i := 0; i < len(line); i++ {
				c := line[i]
				if c == ' ' || c == '\t' {
					continue
				}
				base, ok := upperBase(c)
				if !ok {
					return nil, fmt.Errorf("record %s: invalid base %q at position %d", current.Name, c, seqLen+1)
				}
				current.Sequence = append(current.Sequence, base)
				seqLen++
			}
		}

		if readErr != nil {
			if readErr == io.EOF {
				break
			}
			return nil, readErr
		}
	}

	return genes, nil
}

// parseHeader builds a Gene from a FASTA header line, taking the name from the first field and optional range=chr:start-end and strand=+/- attributes.
func parseHeader(line string) *Gene {
	g := &Gene{Header: line}
	body := strings.TrimPrefix(line, ">")
	if fields := strings.Fields(body); len(fields) > 0 {
		g.Name = fields[0]
	}
	if m := headerRangeRegex.FindStringSubmatch(body); m != nil {
		g.Loc.Chromosome = m[1]
		g.Loc.Start, _ = strconv.ParseInt(m[2], 10, 64)
		g.Loc.End, _ = strconv.ParseInt(m[3], 10, 64)
	}
	if m := headerStrandRegex.FindStringSubmatch(body); m != nil {
		g.Loc.Strand = m[1]
	}
	return g
}

// upperBase uppercases c and reports whether the result is one of A, C, G, T; the returned byte is the uppercased input even when invalid.
func upperBase(c byte) (byte, bool) {
	if c >= 'a' && c <= 'z' {
		c -= 'a' - 'A'
	}
	switch c {
	case 'A', 'T', 'C', 'G':
		return c, true
	}
	return c, false
}

var complementBase = map[byte]byte{'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

// Complement replaces each base of the sequence with its Watson-Crick complement in place.
func (g *Gene) Complement() {
	for i, b := range g.Sequence {
		g.Sequence[i] = complementBase[b]
	}
}

// Invert reverses the sequence in place without complementing it.
func (g *Gene) Invert() {
	for i, j := 0, len(g.Sequence)-1; i < j; i, j = i+1, j-1 {
		g.Sequence[i], g.Sequence[j] = g.Sequence[j], g.Sequence[i]
	}
}

// ApplyOrientation first orients the sequence by Loc.Strand (complement for +, reverse for -), then applies the complement, reverse and invert flags, where invert implies both.
func (g *Gene) ApplyOrientation(complement, reverse, invert bool) {
	switch g.Loc.Strand {
	case "+":
		g.Complement()
	case "-":
		g.Invert()
	}
	if complement || invert {
		g.Complement()
	}
	if reverse || invert {
		g.Invert()
	}
}
