package genome

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"regexp"
	"strconv"
	"strings"
)

// Gene bioinformatic representation of a gene
type Gene struct {
	GeneName string
	Header   string
	Pos      Loci
	Sequence []rune
	// moved vector<Structure> and ground_state_energy to ensemble
}

type FastaHeader struct {
	GeneName      string
	BasePairRange string
	Start         int64
	End           int64
	FivePad       int
	ThreePad      int
	Strand        string
	RepeatMasking string
}

func fileLineScanner(filename string) []string {

	var values []string
	file, err := os.Open(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer func(file *os.File) {
		err := file.Close()
		if err != nil {
			log.Fatal(err)
		}
	}(file)

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if line == "" {
			continue
		}
		values = append(values, line)
	}
	if len(values) == 0 {
		log.Fatal("ERROR: unable to read lines from input file")
	}
	return values
}

func atoiToInt64(s string) int64 {
	n, err := strconv.Atoi(s)
	if err != nil {
		log.Print("WARN: invalid strconv.Atoi result: ", n)
	}
	return int64(n)
}

func parseHeader(header string) (*FastaHeader, error) {

	header = strings.TrimPrefix(header, ">")
	fields := strings.Fields(header)

	if len(fields) < 2 {
		return nil, fmt.Errorf("invalid FASTA Header format")
	}

	geneName := fields[0]
	parsed := &FastaHeader{GeneName: geneName}

	rangeRegex := regexp.MustCompile(`range=([^:]+):(\d+)-(\d+)`)
	padRegex := regexp.MustCompile(`([53]'?pad)=(\d+)`)
	strandRegex := regexp.MustCompile(`strand=([-+])`)
	repeatMaskRegex := regexp.MustCompile(`repeatMasking=([a-zA-Z0-9]+)`) // Adjusted to be more flexible

	for _, field := range fields[1:] {
		if matches := rangeRegex.FindStringSubmatch(field); matches != nil {
			parsed.BasePairRange = matches[1] + ":" + matches[2] + "-" + matches[3]
			parsed.Start = atoiToInt64(matches[2])
			parsed.End = atoiToInt64(matches[3])
		}
		if matches := padRegex.FindStringSubmatch(field); matches != nil {
			val, _ := strconv.Atoi(matches[2])
			if strings.Contains(matches[1], "5'") {
				parsed.FivePad = val
			} else {
				parsed.ThreePad = val
			}
		}
		if matches := strandRegex.FindStringSubmatch(field); matches != nil {
			parsed.Strand = matches[1]
		}
		if matches := repeatMaskRegex.FindStringSubmatch(field); matches != nil {
			parsed.RepeatMasking = matches[1]
		}
	}

	return parsed, nil
}

func NewGene(filename string) *Gene {

	lines := fileLineScanner(filename)
	header, err := parseHeader(lines[0])
	if err != nil {
		log.Fatal("ERROR: unable to parse FASTA Header for input: ", filename)
	}

	var seq []rune
	for _, line := range lines[1:] {
		for i := 0; i < len(line); i++ {
			c := rune(strings.ToUpper(line)[i])
			if c == 'A' || c == 'T' || c == 'C' || c == 'G' {
				seq = append(seq, c)
			} else if c == '\n' || c == ' ' || c == '\t' || c == '\r' {
				continue // TODO: allow multi-gene FASTA files
			} else {
				log.Println("WARN: unrecognized character in input file: ", c)
			}
		}
	}

	return &Gene{
		GeneName: header.GeneName,
		Header:   lines[0],
		Pos: Loci{
			chromosome: "", // TODO: parse from Header
			strand:     header.Strand,
			startPos:   header.Start,
			endPos:     header.End,
		},
		Sequence: seq,
	}
}
