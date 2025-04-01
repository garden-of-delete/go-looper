package genome

// Windower is an interface for
type Windower interface {
	NextWindow() (int, int, bool)
	GetState() *WindowerState
	GetConfig() *WindowerConfig
}

type WindowerConfig struct {
	Sequence      []rune
	MinWindowSize int
	IsCircular    bool
}

type WindowerState struct {
	CurrentStartIndex int
	CurrentStopIndex  int
}

type FromAllWindows struct {
	conf  WindowerConfig
	state WindowerState
}

type SlidingWindow struct {
	conf  WindowerConfig
	state WindowerState
}

type WindowerOption func(windower *Windower) // functional option pattern

func (w *FromAllWindows) getState() *WindowerState {
	return &w.state
}

func (w *FromAllWindows) getConfig() *WindowerConfig {
	return &w.conf
}

//func WithMinWindowSize(size int) WindowerOption {
//	return func(w *Windower) {
//		if size < 2 {
//			log.Fatal("ERROR: minWindowSize must be >= 2")
//		}
//
//		w.MinWindowSize = size
//		w.state.CurrentStopIndex = w.state.CurrentStartIndex + w.conf.MinWindowSize - 1
//	}
//}

func WithCircular(v bool) WindowerOption {
	return func(w *FromAllWindows) {
		w.conf.IsCircular = v
	}
}

func NewWindowerFromAllWindows(s []rune, options ...WindowerOption) *FromAllWindows {
	wc := WindowerConfig{
		Sequence:      s,
		MinWindowSize: 2,
		IsCircular:    false,
	}
	ws := WindowerState{
		CurrentStartIndex: 0,
		CurrentStopIndex:  wc.MinWindowSize - 1,
	}
	w := FromAllWindows{
		conf:  wc,
		state: ws,
	}
	for _, opt := range options {
		opt(&w)
	}
	return &w
}

//func (w *FromAllWindows) NextWindow() (newStart int, newEnd int, hasNextWindow bool) {
//
//	if !w.isCircular { // if sequence not circular
//		if w.currentStartIndex == w.currentStopIndex-w.minWindowSize+1 &&
//			w.currentStopIndex == len(w.Sequence)-1 { // if not at final window
//			return // hasNextWindow = false by default
//		} else { // not circular and next window exists
//			if w.currentStopIndex < len(w.Sequence)-1
//		}
//	}
//}
