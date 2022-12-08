package mTalib

import (
	"fmt"

	"github.com/EasyGolang/goTools/mCount"
	"github.com/EasyGolang/goTools/mTalib/talib"
)

func MA(opt CListOpt) string {
	n := opt.Period
	cLen := len(opt.CList)

	if cLen < n+1 {
		return "0"
	}

	dotNum := mCount.GetDecimal(opt.CList[0]) // 计算小数点位数
	var floatList []float64
	for _, val := range opt.CList {
		valDot := mCount.GetDecimal(val) // 计算当前的小数点位数
		if valDot > dotNum {             // 如果当前小数点位数大于现存小数点位数，则替换
			dotNum = valDot
		}
		floatVal := mCount.ToFloat(val, dotNum)
		floatList = append(floatList, floatVal) // 将数值完整的转化
	}

	// 计算 MA 指标
	pArr := talib.Sma(floatList, n)
	emaFloat := pArr[cLen-1]

	// 保留精确度，并转为字符串
	maStr := fmt.Sprintf("%f", emaFloat)
	maStr = mCount.CentRound(maStr, dotNum)

	return maStr
}

func (_this *ClistObj) MA() *ClistObj {
	if _this.CLen < _this.Period+1 {
		return _this
	}
	pArr := talib.Sma(_this.FList, _this.Period)
	_this.Result = pArr[_this.CLen-1]
	return _this
}

/**

// 这个库废弃了, 用table 重写

func MA(opt CListOpt) string {
	KDList := opt.CList
	n := opt.Cycle

	c_len := len(KDList) // K线总长
	c_n := n             // 长度
	if c_len < n {
		c_n = c_len
	}

	c_list := KDList[c_len-c_n:]
	ma_add := "0"

	if len(opt.Precision) < 1 {
		opt.Precision = KDList[0]
	}
	dotNum := mCount.GetDecimal(opt.Precision) // 计算小数点位数

	for _, KD := range c_list {
		// 数据源
		C := KD

		ma_add = mCount.Add(ma_add, C)
	}

	maRe := mCount.Div(ma_add, mStr.ToStr(c_n))

	maRe = mCount.CentRound(maRe, dotNum)

	return maRe
}

*/
