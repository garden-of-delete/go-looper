package sidd

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"strconv"
)

func FormatSIST(v float64) string {
	switch {
	case math.IsNaN(v):
		return "nan"
	case math.IsInf(v, 1):
		return "inf"
	case math.IsInf(v, -1):
		return "-inf"
	}
	return strconv.FormatFloat(v, 'g', 6, 64)
}

func WriteReport(w io.Writer, inputName string, m *Model, r *Result) error {
	bw := bufio.NewWriter(w)
	transition := "Transition: Melting - copolymeric energetics"
	if m.NearestNeighbor {
		transition = "Energetics: Melting - nearest neighbor energetics"
	}
	molecule := "linear DNA"
	if m.Geometry == GeomCircular {
		molecule = "circular DNA"
	}
	fmt.Fprintf(bw, "DNA sequence file: %s\n", inputName)
	fmt.Fprintf(bw, "Sequence Length = %d\n", m.SeqLen)
	fmt.Fprintf(bw, "%s\n", transition)
	fmt.Fprintf(bw, "Molecule type: %s\n", molecule)
	fmt.Fprintf(bw, "Stress Level = %s\n", FormatSIST(m.Sigma))
	fmt.Fprintf(bw, "linking difference = %s TURNS\n", FormatSIST(m.Alpha))
	fmt.Fprintf(bw, "Threshold = %s kcal/mol\n", FormatSIST(m.Threshold))
	fmt.Fprintf(bw, "Temperature = %s K\n", FormatSIST(m.Temperature))
	fmt.Fprintf(bw, "Salt Concentration = %s M\n", FormatSIST(m.Salt))
	fmt.Fprintf(bw, "C = %s kcal-bp/rad^2\n", FormatSIST(C))
	fmt.Fprintf(bw, "Nucleation energy = %s kcal/mol\n", FormatSIST(Nucleation))
	fmt.Fprintf(bw, "BAT = %s kcal/mol; BGC = %s kcal/mol\n", FormatSIST(m.BAT), FormatSIST(m.BGC))
	fmt.Fprintf(bw, "MinE = %s kcal/mol; Max_E = %s kcal/mol\n", FormatSIST(r.MinE), FormatSIST(r.MaxE))
	fmt.Fprintf(bw, "Number of one-run states = %d\n", r.OneRunStates)
	fmt.Fprintf(bw, "Number of melted bases = %s\n", FormatSIST(r.MeltedBases))
	fmt.Fprint(bw, "Position\tBase\tP(x)\tG(x)\n")
	for i := 0; i < m.SeqLen; i++ {
		fmt.Fprintf(bw, "%d\t%c\t%s\t%s\n", i+1, m.Seq[i], FormatSIST(r.P[i]), FormatSIST(r.G[i]))
	}
	return bw.Flush()
}
