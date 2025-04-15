package rlooper

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

func fileLineScanner(filename string) ([]string, error) {

	var values []string
	file, err := os.Open(filename)
	if err != nil {
		return values, err
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
	return values, nil
}

func atoiToInt64(s string) (int64, error) {
	n, err := strconv.Atoi(s)
	if err != nil {
		return 0, fmt.Errorf("invalid strconv.Atoi result: %v", err)
	}
	return int64(n), nil
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
	strandRegex := regexp.MustCompile(`Strand=([-+])`)
	repeatMaskRegex := regexp.MustCompile(`repeatMasking=([a-zA-Z0-9]+)`) // Adjusted to be more flexible

	for _, field := range fields[1:] {
		if matches := rangeRegex.FindStringSubmatch(field); matches != nil {
			parsed.BasePairRange = matches[1] + ":" + matches[2] + "-" + matches[3]
			start, err := atoiToInt64(matches[2])
			if err != nil {
				return nil, fmt.Errorf("invalid start position: %v", err)
			}
			end, err := atoiToInt64(matches[3])
			if err != nil {
				return nil, fmt.Errorf("invalid end position: %v", err)
			}
			parsed.Start = start
			parsed.End = end
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

// TODO: support for reading multiple Genes out of a single FASTA file
func NewGene(filename string) *Gene {

	lines, err := fileLineScanner(filename)
	if err != nil {
		log.Fatal("ERROR: unable to read lines from input file: ", filename)
	}
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
	if len(seq) < 2 { // from here on, seq should always be initialized and len(seq) > 1
		log.Fatal("ERROR: can't construct a gene with an empty sequence")
	}

	return &Gene{
		GeneName: header.GeneName,
		Header:   lines[0],
		Pos: Loci{
			Chromosome: "", // TODO: parse from Header
			Strand:     header.Strand,
			StartPos:   header.Start,
			EndPos:     header.End,
		},
		Sequence: seq,
	}
}

func (g *Gene) printGene() {
	fmt.Println(g.Header)
	for _, v := range g.Sequence {
		fmt.Print(v)
	}
	fmt.Print('\n')
}

// computeStructuresSerial computes structures the rlooper2 way, which is to say serially in a single thread
// for performance comparison with computeStructures
func (g *Gene) computeStructuresSerial(model *ModelParams, minLoopLength int, circular bool) []Structure {

	windows := FromLinearWindows(g.Sequence, minLoopLength)
	if circular {
		windows = append(windows, FromCircularWindows(g.Sequence, minLoopLength)...)
	}
	var result []Structure
	for _, w := range windows {
		structure := Structure{
			Pos: Loci{
				g.Pos.Chromosome,
				g.Pos.Strand,
				int64(w.Start), // TODO: loci Pos is in terms of genomic coordinates in rlooper2
				int64(w.End),
			},
			FreeEnergy:      0,
			BoltzmannFactor: 0,
			Probability:     0,
		}
		model.ComputeStructure(g.Sequence, w, &structure)
		result = append(result, structure)
	}
	return result
}

func (g *Gene) computeStructuresConcurrent(ec *ExecutionContext, model *ModelParams, minLoopLength int, circular bool) []Structure {
	windows := FromLinearWindows(g.Sequence, minLoopLength)
	if circular {
		windows = append(windows, FromCircularWindows(g.Sequence, minLoopLength)...)
	}

	// Handle case where no threads are requested
	if ec.NumThreads <= 0 {
		// Fall back to serial computation
		return g.computeStructuresSerial(model, minLoopLength, circular)
	}

	// Calculate block size once outside the loop, edge case where n windows are less than numThreads
	blockSize := len(windows) / ec.NumThreads
	if blockSize == 0 {
		blockSize = 1
	}

	// Create a channel with a reasonable buffer size
	structureChan := make(chan []Structure, ec.NumThreads)

	// Create an empty slice for results
	results := make([]Structure, 0)

	ec.WaitGroup.Add(ec.NumThreads)
	for i := 0; i < ec.NumThreads; i++ { //compute structures in parallel
		go func(i int) {
			defer ec.WaitGroup.Done()

			var chunk []Window
			if i != ec.NumThreads-1 {
				chunk = windows[i*blockSize : (i+1)*blockSize]
			} else {
				chunk = windows[i*blockSize:]
			}

			// Collect structures for this chunk
			chunkStructures := make([]Structure, 0, len(chunk))
			for _, w := range chunk {
				structure := Structure{
					Pos: Loci{
						g.Pos.Chromosome,
						g.Pos.Strand,
						int64(w.Start), // TODO: loci Pos is in terms of genomic coordinates in rlooper2
						int64(w.End),
					},
					FreeEnergy:      0,
					BoltzmannFactor: 0,
					Probability:     0,
				}
				model.ComputeStructure(g.Sequence, w, &structure)
				chunkStructures = append(chunkStructures, structure)
			}
			structureChan <- chunkStructures
		}(i)
	}

	// Start a goroutine to close the channel when all workers are done
	go func() {
		ec.WaitGroup.Wait()
		close(structureChan)
	}()

	// Collect all structures from the channel
	for chunk := range structureChan {
		results = append(results, chunk...)
	}

	return results
}
