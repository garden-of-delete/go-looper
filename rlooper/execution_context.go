package rlooper

import "sync"

type ExecutionContext struct {
	NumThreads int
	WaitGroup  *sync.WaitGroup
	Concurrency bool

}

type SafeFloat64 struct {
    mu sync.Mutex
    val float64
}
