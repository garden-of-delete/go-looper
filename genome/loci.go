package genome

type Loci struct {
	chromosome string
	strand     string
	startPos   int64
	endPos     int64
}

func (L *Loci) getLength() int64 {
	return L.endPos - L.startPos
}
