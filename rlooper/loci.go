package rlooper

type Loci struct {
	Chromosome string
	Strand     string
	StartPos   int64
	EndPos     int64
}

func (L *Loci) getLength() int64 {
	return L.EndPos - L.StartPos
}
