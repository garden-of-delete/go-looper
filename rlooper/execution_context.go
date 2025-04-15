package rlooper

import "sync"

type ExecutionContext struct {
	NumThreads int
	WaitGroup  *sync.WaitGroup
}

type SafeFloat64 struct {
    mu sync.Mutex
    val float64
}
