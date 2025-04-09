package src

type Window struct {
	Start int
	End   int
}

// FromLinearWindows generates all index ranges as Window >= min window length (mll).
// all possible structures >= minLoopLength
func FromLinearWindows(seq []rune, minLoopLength int) []Window {
	var result []Window
	for i := range seq {
		for j := i + minLoopLength - 1; j < len(seq); j++ {
			result = append(result, Window{Start: i, End: j})
		}
	}
	if len(result) > 0 {
		return result
	} else {
		return []Window{}
	}
}

// FromCircularWindows generates all index ranges that include the circular boundary between
// the beginning and end of the input sequence.
// FromCircularWindows union FromLinearWindows should produce all possible
// windows > minLoopLength on some input sequence.
func FromCircularWindows(seq []rune, minLoopLength int) []Window {
	var result []Window
	//for i := len(seq) - 1 - minLoopLength; i < len(seq); i++ {
	//	for j := i + minLoopLength - 1; j < minLoopLength; j++ {
	//		result = append(result, Window{Start: i, End: minLoopLength % j})
	//	}
	//}
	for i := 1; i < len(seq); i++ {
		for j := 0; j < i-minLoopLength+2; j++ { // TODO: does this work?
			length := (len(seq) - i) + (j + 1)
			if !(length < minLoopLength) {
				result = append(result, Window{i, j})
			}
			//fmt.Println(j)
			//if j == len(seq) { // if over the end of the sequence, return to start
			//	j = 0
			//}
		}
	}
	return result
}
