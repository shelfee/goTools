package mOKX

import (
	"github.com/EasyGolang/goTools/mCount"
)

// 按照 最低 价排序
func Sort_HLPer(data []TypeKd) []TypeKd {
	size := len(data)
	list := make([]TypeKd, size)
	copy(list, data)

	var swapped bool
	for i := size - 1; i > 0; i-- {
		swapped = false
		for j := 0; j < i; j++ {
			a := list[j+1].HLPer
			b := list[j].HLPer
			if mCount.Le(a, b) < 0 {
				list[j], list[j+1] = list[j+1], list[j]
				swapped = true
			}
		}
		if !swapped {
			break
		}
	}
	// 设置 Idx 并翻转
	listIDX := []TypeKd{}
	j := 0
	for i := len(list) - 1; i > -1; i-- {
		Kdata := list[i]
		Kdata.HLPer_idx = j + 1
		listIDX = append(listIDX, Kdata)
		j++
	}
	return listIDX
}
