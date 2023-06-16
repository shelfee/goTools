package mCount

import (
	"fmt"
	"github.com/shopspring/decimal"
	"sync"
	"testing"
	"time"
)

func TestToFloat(t *testing.T) {
	n, _ := decimal.NewFromString("12112.34325235")
	startTm := time.Now()
	w := sync.WaitGroup{}
	w.Add(10)
	for i := 0; i < 10; i++ {
		go func(i int) {
			for j := 0; j < 10000000; j++ {
				n.Float64()
			}
			w.Done()
			fmt.Println("done", i)

		}(i)
	}
	w.Wait()
	fmt.Println("latency:", time.Now().Sub(startTm).Milliseconds(), "ms")
}
