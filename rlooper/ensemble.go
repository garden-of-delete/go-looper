package rlooper

import (
	"cmp"
	"math"
	"slices"
	"sync"
)

type EnsembleResult struct {
	LogZ, LogGround, PGround float64
	MFE                      StructureRecord
	NumStructures            int
}

type StructureRecord struct {
	// TODO: add type to support multi-structure competition
	Start, End int
	Energy     float64
	Prob       float64
}

type logSumExp struct {
	max, sum float64
}

// newLogSumExp returns an empty accumulator whose max is -Inf, representing log(0).
func newLogSumExp() logSumExp { return logSumExp{max: math.Inf(-1)} }

// add folds x into the running log-sum-exp, rescaling the sum when x exceeds the current max.
func (l *logSumExp) add(x float64) {
	if x <= l.max {
		l.sum += math.Exp(x - l.max)
		return
	}
	l.sum = l.sum*math.Exp(l.max-x) + 1
	l.max = x
}

// merge folds another accumulator into l; it combines per-worker partial sums and treats an empty (-Inf) operand as identity.
func (l *logSumExp) merge(o logSumExp) {
	if math.IsInf(o.max, -1) {
		return
	}
	if math.IsInf(l.max, -1) {
		*l = o
		return
	}
	if o.max <= l.max {
		l.sum += o.sum * math.Exp(o.max-l.max)
		return
	}
	l.sum = l.sum*math.Exp(l.max-o.max) + o.sum
	l.max = o.max
}

// value returns log of the accumulated sum, or -Inf when nothing was added.
func (l *logSumExp) value() float64 {
	if math.IsInf(l.max, -1) {
		return math.Inf(-1)
	}
	return l.max + math.Log(l.sum)
}

// compareStructures orders records by energy, then Start, then End; the Start/End tie-breaks make the MFE choice match the serial enumeration order.
func compareStructures(a, b StructureRecord) int {
	return cmp.Or(
		cmp.Compare(a.Energy, b.Energy),
		cmp.Compare(a.Start, b.Start),
		cmp.Compare(a.End, b.End),
	)
}

// windowEnergy returns the free energy of the mLen-base window starting at s as supercoiling plus the prefix-sum difference p[s+mLen-1]-p[s]; the subtraction is grouped before the add, reassociating it shifts results by ulps, and no test covers that.
func windowEnergy(m *EnergyModel, p []float64, s, mLen int) float64 {
	return m.SupercoilingEnergy(mLen) + (p[s+mLen-1] - p[s])
}

type Windows struct {
	MinLength int
	Circular  bool
	Workers   int
}

// maxLen returns the longest window length permitted at start s: n-1 on a circle, n-s on a line.
func maxLen(n, s int, circular bool) int {
	if circular {
		return n - 1
	}
	return n - s
}

// windowEnd returns the inclusive end index of the mLen-base window starting at s, wrapping modulo n when circular.
func windowEnd(s, mLen, n int, circular bool) int {
	e := s + mLen - 1
	if circular {
		e %= n
	}
	return e
}

// clampWorkers caps the worker count at the number of rows so no worker is created without work.
func clampWorkers(workers, rows int) int {
	if workers > rows {
		return rows
	}
	return workers
}

type ensembleAccum struct {
	lse   logSumExp
	mfe   StructureRecord
	count int
}

// newEnsembleAccum returns an empty per-worker ensemble accumulator.
func newEnsembleAccum() ensembleAccum {
	return ensembleAccum{lse: newLogSumExp()}
}

// addWindow records one window with free energy g, adding its Boltzmann weight to the partition sum and updating the MFE under compareStructures.
func (a *ensembleAccum) addWindow(s, e int, g, rt float64) {
	a.lse.add(-g / rt)
	rec := StructureRecord{Start: s, End: e, Energy: g}
	if a.count == 0 || compareStructures(rec, a.mfe) < 0 {
		a.mfe = rec
	}
	a.count++
}

// merge combines a worker's partial accumulator into a, preserving the same MFE tie-break as serial enumeration.
func (a *ensembleAccum) merge(o ensembleAccum) {
	a.lse.merge(o.lse)
	if o.count == 0 {
		return
	}
	if a.count == 0 || compareStructures(o.mfe, a.mfe) < 0 {
		a.mfe = o.mfe
	}
	a.count += o.count
}

// ComputeEnsemble enumerates every window permitted by w and returns the partition function (including the ground state), ground-state probability, and the MFE structure with its probability.
func ComputeEnsemble(seq []byte, m *EnergyModel, w Windows) EnsembleResult {
	p := prefixBP(seq, m, w.Circular)
	rt := m.RT()
	n := len(seq)
	workers := clampWorkers(w.Workers, n)

	acc := newEnsembleAccum()
	if workers <= 1 {
		for s := 0; s < n; s++ {
			for mLen := w.MinLength; mLen <= maxLen(n, s, w.Circular); mLen++ {
				acc.addWindow(s, windowEnd(s, mLen, n, w.Circular), windowEnergy(m, p, s, mLen), rt)
			}
		}
	} else {
		results := make([]ensembleAccum, workers)
		var wg sync.WaitGroup
		for c := 0; c < workers; c++ {
			wg.Add(1)
			go func(c int) {
				defer wg.Done()
				chunk := newEnsembleAccum()
				for s := c; s < n; s += workers {
					for mLen := w.MinLength; mLen <= maxLen(n, s, w.Circular); mLen++ {
						chunk.addWindow(s, windowEnd(s, mLen, n, w.Circular), windowEnergy(m, p, s, mLen), rt)
					}
				}
				results[c] = chunk
			}(c)
		}
		wg.Wait()
		for _, r := range results {
			acc.merge(r)
		}
	}

	logGround := -m.GroundStateEnergy() / rt
	acc.lse.add(logGround)
	logZ := acc.lse.value()

	if acc.count > 0 {
		acc.mfe.Prob = math.Exp(-acc.mfe.Energy/rt - logZ)
	}

	return EnsembleResult{
		LogZ:          logZ,
		LogGround:     logGround,
		PGround:       math.Exp(logGround - logZ),
		MFE:           acc.mfe,
		NumStructures: acc.count,
	}
}

// AllStructures returns every enumerated structure with its Boltzmann probability under ens, sorted by compareStructures.
func AllStructures(seq []byte, m *EnergyModel, w Windows, ens EnsembleResult) []StructureRecord {
	p := prefixBP(seq, m, w.Circular)
	rt := m.RT()
	n := len(seq)
	recs := make([]StructureRecord, 0, ens.NumStructures)

	for s := 0; s < n; s++ {
		for mLen := w.MinLength; mLen <= maxLen(n, s, w.Circular); mLen++ {
			g := windowEnergy(m, p, s, mLen)
			recs = append(recs, StructureRecord{
				Start:  s,
				End:    windowEnd(s, mLen, n, w.Circular),
				Energy: g,
				Prob:   math.Exp(-g/rt - ens.LogZ),
			})
		}
	}

	slices.SortFunc(recs, compareStructures)

	return recs
}

// TopStructures returns the n lowest-energy structures in compareStructures order using a bounded insertion sort, without materializing all structures.
func TopStructures(seq []byte, m *EnergyModel, w Windows, ens EnsembleResult, n int) []StructureRecord {
	if n <= 0 {
		return nil
	}
	p := prefixBP(seq, m, w.Circular)
	rt := m.RT()
	seqLen := len(seq)
	top := make([]StructureRecord, 0, n)

	for s := 0; s < seqLen; s++ {
		for mLen := w.MinLength; mLen <= maxLen(seqLen, s, w.Circular); mLen++ {
			g := windowEnergy(m, p, s, mLen)
			r := StructureRecord{Start: s, End: windowEnd(s, mLen, seqLen, w.Circular), Energy: g, Prob: math.Exp(-g/rt - ens.LogZ)}
			if len(top) == n && compareStructures(r, top[n-1]) >= 0 {
				continue
			}
			i := len(top)
			if i == n {
				i--
			} else {
				top = append(top, StructureRecord{})
			}
			for i > 0 && compareStructures(top[i-1], r) > 0 {
				top[i] = top[i-1]
				i--
			}
			top[i] = r
		}
	}

	return top
}
