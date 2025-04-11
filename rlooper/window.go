package rlooper

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
	for i := 1; i < len(seq); i++ {
		for j := 0; j < i; j++ {
			length := (len(seq) - i) + j + 1
			if length >= minLoopLength {
				result = append(result, Window{Start: i, End: j})
			}
		}
	}
	if len(result) > 0 {
		return result
	} else {
		return []Window{}
	}
}

func WindowToString(seq []rune, w Window) string {
	var result string
	for i := w.Start; i != w.End+1; {
		result += string(seq[i])
		if i == len(seq)-1 {
			if w.End == len(seq)-1 {
				break // Exit if we've reached the end
			}
			i = 0 // Loop back to start of sequence
		} else {
			i++
		}
	}
	return result
}

func PrintWindows(seq []rune, windows []Window) []string {
	var result []string
	for _, w := range windows {
		result = append(result, WindowToString(seq, w))
	}
	return result
}
