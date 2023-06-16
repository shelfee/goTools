package zDec

import (
	"bytes"
	"errors"
	"fmt"
	"github.com/shopspring/decimal"
	"math"
	"strconv"
	"strings"
)

type Number interface {
	Add(number Number) Number
	Sub(number Number) Number
	Mul(number Number) Number
	Div(number Number) Number
	Neg() Number
	String() string
}

func absOp(i int64) uint64 {
	if i < 0 {
		return uint64(-i)
	}
	return uint64(i)
}
func shrinkInt(v int64) (int64, int) {
	t := v
	i := 0
	for t%10 == 0 && t != 0 {
		t /= 10
		i++
	}
	return t, i
}

func shrinkBigInt(v decimal.Decimal) (decimal.Decimal, int) {
	i := 0
	if !v.IsInteger() {
		for !v.IsInteger() && i < NanoOffset {
			v = v.Shift(1)
			i -= 1
		}
		return v, i
	}
	dv := decimal.NewFromInt(10)
	for !v.IsZero() && v.Mod(dv).IsZero() {
		v = v.Shift(-1)
		i += 1
	}
	return v, i
}

type Float64 float64

func (f Float64) Add(n Number) Number {
	switch n.(type) {
	case DynNum:
		r := NewDynNumFromFloat64(float64(f))
		return r.Add(n)
	case Float64:
		v := n.(Float64)
		return f + v
	default:
		panic("Float64 type not support")
	}
	return f
}
func (f Float64) Sub(n Number) Number {
	switch n.(type) {
	case DynNum:
		r := NewDynNumFromFloat64(float64(f))
		return r.Sub(n)
	case Float64:
		v := n.(Float64)
		return f - v
	default:
		panic("Float64 type not support")
	}
	return f
}
func (f Float64) Mul(n Number) Number {
	switch n.(type) {
	case DynNum:
		return n.Mul(f)
	case Float64:
		v := n.(Float64)
		return Float64(float64(f) * float64(v))
	default:
		panic("Float64 type not support")
	}
}
func (f Float64) Div(n Number) Number {
	switch n.(type) {
	case Float64:
		v := n.(Float64)
		return Float64(float64(f) / float64(v))
	default:
		panic("Float64 type not support")
	}
}
func (f Float64) Neg() Number {
	return -f
}

func (f Float64) String() string {
	return fmt.Sprintf("%f", f)
}

// DynNum DynNum: support unit as nano and atto, nano is for normal number (>10e-2 ~ 10e10)
type DynNum struct {
	state int // 0: use nano, 1: use nanoBigNum, 2: use atto, 3: use attoBigNum

	nano       int64           // if num * 10^9 is integer, represent num * 10^9 and set state is 0
	nanoBigNum decimal.Decimal // if num * 10^9 is larger than 2^63, use decimal.Decimal instead of int

	atto       int64           // if num * 10^18 is integer, represent num * 10^18 and set state is 0
	attoBigNum decimal.Decimal // if num * 10^18 is larger than 2^63, use decimal.Decimal instead of int
}

const NanoOffset = 9
const NanoUnit int64 = 10e8
const BoardNanoVal int64 = 4611686018

var NanoBigNum = decimal.NewFromInt(NanoUnit)
var AttoBigNum = decimal.NewFromInt(NanoUnit).Mul(NanoBigNum)
var BoardBigNum = decimal.NewFromInt(math.MaxInt64 / 2)

var MinPreBigNum = decimal.NewFromInt(BoardNanoVal)

func (d DynNum) ShrinkInt() (interface{}, int) {
	if d.state == 0 {
		v, i := shrinkInt(d.nano)
		return v, i - NanoOffset
	}
	if d.state == 2 {
		v, i := shrinkInt(d.atto)
		return v, i - NanoOffset*2
	}
	if d.state == 1 {
		v, i := shrinkBigInt(d.nanoBigNum)
		return v, i - NanoOffset
	}
	v, i := shrinkBigInt(d.attoBigNum)
	return v, i - NanoOffset*2

}

func NewDynNumFromString(s string) (DynNum, error) {
	r := DynNum{}
	if strings.Trim(s, "0.") == "" {
		return r, nil
	}
	v, err := decimal.NewFromString(s)
	if err != nil {
		return r, err
	}
	r.nanoBigNum = v.Mul(NanoBigNum)
	if r.nanoBigNum.IsInteger() {
		r.state = 1
		r.ShrinkBigNum()
		return r, nil
	}
	r.attoBigNum = r.nanoBigNum.Mul(NanoBigNum)
	r.state = 3
	r.ShrinkBigNum()
	return r, nil
}

func NewDynNumFromInt(i int64) DynNum {
	if absOp(i) < 9223372036 {
		return DynNum{
			state: 0,
			nano:  i * NanoUnit,
		}
	}
	return DynNum{
		state:      1,
		nanoBigNum: decimal.NewFromInt(i).Mul(NanoBigNum),
	}
}

func NewDynNumFromFloat64(i float64) DynNum {
	if math.Abs(i) < 1e-6 {
		return DynNum{
			state: 2,
			nano:  int64(math.Round(i * float64(NanoUnit) * float64(NanoUnit))),
		}
	}
	if math.Abs(i) < 9223372036 {
		return DynNum{
			state: 0,
			nano:  int64(math.Round(i * float64(NanoUnit))),
		}
	}
	return DynNum{
		state:      1,
		nanoBigNum: decimal.NewFromFloat(i).Mul(NanoBigNum),
	}
}

func (d *DynNum) IsZero() bool {
	if d.state == 0 {
		return d.nano == 0
	}
	if d.state == 1 {
		return d.nanoBigNum.IsZero()
	}
	if d.state == 2 {
		return d.atto == 0
	}
	if d.state == 3 {
		return d.attoBigNum.IsZero()
	}
	panic("invalid state")
}

func (d *DynNum) IsPositive() bool {
	if d.state == 0 {
		return d.nano > 0
	}
	if d.state == 1 {
		return d.nanoBigNum.IsPositive()
	}
	if d.state == 2 {
		return d.atto > 0
	}
	if d.state == 3 {
		return d.attoBigNum.IsPositive()
	}
	panic("invalid state")
}

func (d *DynNum) IsNegative() bool {
	if d.state == 0 {
		return d.nano < 0
	}
	if d.state == 1 {
		return d.nanoBigNum.IsNegative()
	}
	if d.state == 2 {
		return d.atto < 0
	}
	if d.state == 3 {
		return d.attoBigNum.IsNegative()
	}
	panic("invalid state")
}

func (d *DynNum) Abs() Number {
	if d.IsNegative() {
		return d.Neg()
	}
	return d
}

// tryPow only support postive power now
func tryPow(v, p float64) (r float64, mov int) {
	if v == 0 {
		r = 0
		mov = 0
		return
	}
	if p == 1 {
		r = v
		mov = 0
		return
	}

	//flag := v > 0
	if math.Round(p)-p < 1.e-3 {

	}
	if math.Abs(p) > 10 {
		exp := math.Log10(v) * p
		mov = int(exp)
		r = math.Pow(10, exp-float64(mov))
		return
	}
	mov = 0
	unit := 4.
	for true {
		r = math.Pow(v, p)
		if !math.IsInf(r, 0) {
			return
		}
		move := unit * p
		mov += int(move)
		v = v / 1e4 * math.Pow(10, move-math.Round(move))
	}
	return
}

//func (d DynNum) Pow(n Number) Number {
//
//	r := DynNum{}
//	switch n.(type) {
//	case Float64:
//		f := n.(Float64)
//		if d.state == 0 || d.state == 2 {
//			v := d.atto
//			off := -2 * NanoOffset
//			if d.state == 0 {
//				v = d.nano
//				off = -NanoOffset
//			}
//			if float64(f) < 1 {
//				newOff := float64(off) * float64(f)
//				digitMov := int64(newOff)
//				extra := math.Pow(10, newOff-float64(digitMov))
//				newVal := math.Pow(float64(v), float64(f)) * extra
//
//				if uint64(math.MaxInt64/NanoUnit) > absOp(v) {
//					val := math.Pow(float64(v*NanoUnit), float64(f))
//					u := 1.
//					if d.state == 2 {
//						u = math.Pow(float64(NanoUnit), float64(f))
//					}
//					res := val * u
//					if math.Abs(res) < float64(BoardNanoVal) {
//						r.atto = int64(math.Round(res * float64(NanoUnit)))
//						r.state = 2
//						return r
//					}
//					if float64(math.MaxInt64/NanoUnit) > math.Abs(res) {
//						r.nano = int64(math.Round(res * float64(NanoUnit)))
//						r.state = 0
//						return r
//					}
//					r.nanoBigNum = decimal.NewFromFloat(res).Shift(NanoOffset)
//					r.state = 1
//					return r
//				} else {
//
//				}
//			}
//		}
//
//	}
//
//}

func (d *DynNum) AttoType() DynNum {
	// transfer to atto unit
	if d.state > 1 {
		return *d
	}
	r := *d
	if d.state == 0 && uint64(math.MaxInt64/NanoUnit) > absOp(d.nano) {
		r.atto = NanoUnit * d.nano
		r.state = 2
		return r
	} else if d.state == 0 {
		r.attoBigNum = decimal.NewFromInt(d.nano).Mul(NanoBigNum)
		r.state = 3
		return r
	} else {
		r.attoBigNum = d.nanoBigNum.Mul(NanoBigNum)
		r.state = 3
		return r
	}
}

func (d *DynNum) ShrinkBigNum() {
	if d.state == 3 {
		if d.attoBigNum.IsInteger() && d.attoBigNum.Mod(NanoBigNum).IsZero() {
			d.nanoBigNum = d.attoBigNum.Shift(-NanoOffset)
			d.state = 1
			d.ShrinkBigNum()
		}
		if d.attoBigNum.GreaterThan(BoardBigNum.Neg()) && d.attoBigNum.LessThan(BoardBigNum) && d.attoBigNum.IsInteger() {
			d.atto = d.attoBigNum.IntPart()
			d.state = 2
		} else {
			d.state = 3
		}
	} else if d.state == 1 {
		if d.nanoBigNum.GreaterThan(BoardBigNum.Neg()) && d.nanoBigNum.LessThan(BoardBigNum) {
			d.nano = d.nanoBigNum.IntPart()
			d.state = 0
		} else {
			d.state = 1
		}
	}
}

func (d *DynNum) FillBigNum() {
	if d.state == 1 || d.state == 3 {
		return
	}
	if d.state == 0 {
		if d.nanoBigNum.IsZero() {
			d.nanoBigNum = decimal.NewFromInt(d.nano)
		}
		return
	}
	if d.attoBigNum.IsZero() {
		d.attoBigNum = decimal.NewFromInt(d.atto)
	}
}

func (d DynNum) Add(n Number) Number {
	r := DynNum{}
	switch n.(type) {
	case Float64:
		v := n.(Float64)
		return d.Add(NewDynNumFromFloat64(float64(v)))
	case DynNum:
		v := n.(DynNum)
		if v.state == 0 && d.state == 0 {
			if uint64(math.MaxInt64)-absOp(d.nano) > absOp(v.nano) {
				r.state = 0
				r.nano = d.nano + v.nano
				return r
			}
			r.state = 1
			r.nanoBigNum = decimal.NewFromInt(d.nano).Add(decimal.NewFromInt(v.nano))
			return r
		}
		if (v.state == 1 && d.state < 2) || (d.state == 1 && d.state < 2) {
			v.FillBigNum()
			d.FillBigNum()
			r.nanoBigNum = v.nanoBigNum.Add(d.nanoBigNum)
			r.state = 1
			r.ShrinkBigNum()
			return r
		}
		// transfer to atto unit
		d1 := d.AttoType()
		v1 := v.AttoType()
		if d1.state == 2 && v1.state == 2 {
			if math.MaxInt64-absOp(d1.atto) > absOp(v1.atto) {
				r.state = 0
				r.atto = d1.atto + v1.atto
				return r
			}
			r.state = 1
			r.attoBigNum = decimal.NewFromInt(d1.atto).Add(decimal.NewFromInt(v1.atto))
			return r
		}
		v1.FillBigNum()
		d1.FillBigNum()
		r.attoBigNum = v1.attoBigNum.Add(d1.attoBigNum)
		r.state = 3
		r.ShrinkBigNum()
		return r
	default:
		panic("Only support DynNum type")
	}
	return r
}
func (d DynNum) Neg() Number {
	r := DynNum{}
	if d.state == 0 && d.nano != math.MinInt64 {
		r.state = 0
		r.nano = -d.nano
		return r
	}
	if d.nanoBigNum.IsZero() {
		d.nanoBigNum = decimal.NewFromInt(d.nano)
	}
	if d.state < 2 {
		r.state = 1
		r.nanoBigNum = d.nanoBigNum.Neg()
	}
	if d.state == 2 && d.atto != math.MinInt64 {
		r.state = 2
		r.atto = -d.atto
	}
	if d.attoBigNum.IsZero() {
		d.attoBigNum = decimal.NewFromInt(d.atto)
	}
	r.state = 3
	r.attoBigNum = d.attoBigNum.Neg()
	return r
}
func (d DynNum) Sub(n Number) Number {
	r := DynNum{}
	switch n.(type) {
	case Float64:
		v := n.(Float64)
		return d.Add(NewDynNumFromFloat64(float64(-v)))
	case DynNum:
		v := n.(DynNum)
		return v.Add(n.Neg())
	default:
		panic("Only support DynNum type")
	}
	return r
}

func (d DynNum) Mul(n Number) Number {
	r := DynNum{}
	switch n.(type) {
	case DynNum:
		v := n.(DynNum)
		if v.IsZero() || d.IsZero() {
			return r
		}
		dvt, dc := d.ShrinkInt()
		vvt, vc := d.ShrinkInt()
		var dBig, vBig decimal.Decimal
		switch dvt.(type) {
		case int64:
			dv := dvt.(int64)
			switch vvt.(type) {
			case int64:
				vv := vvt.(int64)
				if math.MaxInt64/absOp(dv) > absOp(vv) {
					z, zc := shrinkInt(dv * vv)
					totalOffset := dc + vc + zc - NanoOffset - NanoOffset
					//1. 退化成atto
					if dc+vc+zc < NanoOffset {
						expPart := int64(math.Pow(10, float64(totalOffset+2*NanoOffset)))
						if math.MaxInt64/z > expPart {
							r.atto = z * expPart
							r.state = 2
						} else {
							r.attoBigNum = decimal.NewFromInt(z).Mul(decimal.NewFromInt(expPart))
							r.state = 3
							r.ShrinkBigNum()
						}
						return r
					}
					expPart := int64(math.Pow(10, float64(totalOffset+NanoOffset)))
					if math.MaxInt64/z > expPart {
						r.nano = z * expPart
						r.state = 0
					} else {
						r.nanoBigNum = decimal.NewFromInt(z).Mul(decimal.NewFromInt(expPart))
						r.state = 1
						r.ShrinkBigNum()
					}
					return r
				}
				decimal.NewFromInt(dv).Mul(decimal.NewFromInt(vv))
				totalOffset := dc + vc - NanoOffset - NanoOffset
				if dc+vc < NanoOffset {
					expPart := int64(math.Pow(10, float64(totalOffset+2*NanoOffset)))
					r.attoBigNum = decimal.NewFromInt(dv).Mul(decimal.NewFromInt(vv)).Mul(decimal.NewFromInt(expPart))
					r.state = 3
					r.ShrinkBigNum()
					return r
				}
				expPart := int64(math.Pow(10, float64(totalOffset+NanoOffset)))
				r.nanoBigNum = decimal.NewFromInt(dv).Mul(decimal.NewFromInt(vv)).Mul(decimal.NewFromInt(expPart))
				r.state = 1
				r.ShrinkBigNum()
				return r
			default:
				vBig = vvt.(decimal.Decimal)
				dBig = decimal.NewFromInt(dv)
			}
		default:
			dBig = dvt.(decimal.Decimal)
			ok := false
			if vBig, ok = vvt.(decimal.Decimal); !ok {
				vBig = decimal.NewFromInt(vvt.(int64))
			}
		}
		r.attoBigNum = dBig.Mul(vBig).Mul(decimal.NewFromInt(10).Pow(decimal.NewFromInt(int64(dc + vc))))
		r.state = 3
		r.ShrinkBigNum()
		return r
	case Float64:
		v := n.(Float64)
		if v == 0 || d.IsZero() {
			return r
		}
		if d.state == 0 {
			if float64(math.MaxInt64) > math.Abs(float64(v))*math.Abs(float64(d.nano)) {
				r.state = 0
				r.nano = int64(math.Round(float64(d.nano) * float64(v)))
				return r
			}
			r.state = 1
			r.nanoBigNum = decimal.NewFromInt(d.nano).Mul(decimal.NewFromFloat(float64(v)))
		}
		if d.state == 1 {
			if d.state == 0 && !d.nanoBigNum.IsZero() {
				d.nanoBigNum = decimal.NewFromInt(d.nano)
			}
			r.nanoBigNum = d.nanoBigNum.Mul(decimal.NewFromFloat(float64(v)))
			r.state = 1
			r.ShrinkBigNum()
			return r
		}
		if d.state == 2 {
			if d.atto == 0 {
				return r
			}
			if float64(math.MaxInt64) > math.Abs(float64(v))*math.Abs(float64(d.atto)) {
				r.state = 2
				r.atto = int64(math.Round(float64(d.atto) * float64(v)))
				return r
			}
			r.state = 1
			r.attoBigNum = decimal.NewFromInt(d.atto).Mul(decimal.NewFromFloat(float64(v)))
		}
		if d.state == 3 {
			if d.state == 0 && !d.attoBigNum.IsZero() {
				d.attoBigNum = decimal.NewFromInt(d.nano)
			}
			r.attoBigNum = d.attoBigNum.Mul(decimal.NewFromFloat(float64(v)))
			r.state = 3
			r.ShrinkBigNum()
			return r
		}

	default:
		panic("type invalid")
	}
	return r
}

func (d DynNum) Div(n Number) Number {
	r := DynNum{}
	switch n.(type) {
	case DynNum:
		v := n.(DynNum)
		if v.IsZero() {
			panic("div zero")
		}
		if d.state == 0 && v.state == 0 {
			r.nano = int64(math.Round(float64(d.nano) / float64(d.nano) * float64(NanoUnit)))
			if r.nano < 1 {
				r.atto = int64(math.Round(float64(NanoUnit) * float64(d.nano) / float64(d.nano) * float64(NanoUnit)))
				r.state = 2
				return r
			}
			r.state = 0
			return r
		}
		if d.state < 2 && v.state < 2 {
			d.FillBigNum()
			v.FillBigNum()
			r.nanoBigNum = d.nanoBigNum.Div(v.nanoBigNum).Mul(NanoBigNum)
			r.state = 1
			r.ShrinkBigNum()
			return r
		}
		if d.state == 2 && v.state == 2 {
			r.nano = int64(math.Round(float64(d.atto) / float64(d.atto) * float64(NanoUnit)))
			if r.nano < 1 {
				r.atto = int64(math.Round(float64(NanoUnit) * float64(d.atto) / float64(d.atto) * float64(NanoUnit)))
				r.state = 2
				return r
			}
			r.state = 0
			return r
		}
		v1 := v.AttoType()
		d1 := d.AttoType()
		if d1.state > 1 && v1.state > 1 {
			d1.FillBigNum()
			v1.FillBigNum()
			r.nanoBigNum = d1.attoBigNum.Div(v1.attoBigNum).Mul(NanoBigNum)
			r.state = 1
			r.ShrinkBigNum()
			return r
		}
	case Float64:
		v := n.(Float64)
		if v == 0 {
			panic("div zero")
		}
		if d.state == 0 {
			if d.nano == 0 {
				return r
			}
			if float64(math.MaxInt64) > math.Abs(float64(d.nano))/math.Abs(float64(v)) {
				r.state = 0
				r.nano = int64(math.Round(float64(d.nano) / float64(v)))
				return r
			}
			r.state = 1
			r.nanoBigNum = decimal.NewFromInt(d.nano).Div(decimal.NewFromFloat(float64(v)))
		}
		if d.state == 1 {
			if d.state == 0 && !d.nanoBigNum.IsZero() {
				d.nanoBigNum = decimal.NewFromInt(d.nano)
			}
			r.nanoBigNum = d.nanoBigNum.Div(decimal.NewFromFloat(float64(v)))
			r.state = 1
			r.ShrinkBigNum()
			return r
		}
		if d.state == 2 {
			if d.atto == 0 {
				return r
			}
			if float64(math.MaxInt64) > math.Abs(float64(d.atto))/math.Abs(float64(v)) {
				r.state = 2
				r.atto = int64(math.Round(float64(d.atto) / float64(v)))
				return r
			}
			r.state = 1
			r.attoBigNum = decimal.NewFromInt(d.atto).Div(decimal.NewFromFloat(float64(v)))
		}
		if d.state == 3 {
			if d.state == 0 && !d.attoBigNum.IsZero() {
				d.attoBigNum = decimal.NewFromInt(d.nano)
			}
			r.attoBigNum = d.attoBigNum.Div(decimal.NewFromFloat(float64(v)))
			r.state = 3
			r.ShrinkBigNum()
			return r
		}
	}
	return d
}

func (d DynNum) String() string {
	op := func(v int64, digitMove int) string {
		if v == 0 {
			return "0"
		}
		var newByte []byte

		s := strconv.FormatInt(v, 10)
		pos := 0
		if v < 0 {
			newByte = make([]byte, len(s)+2)
			newByte[pos] = '-'
			pos++
			v = -v
		} else {
			newByte = make([]byte, len(s)+1)
		}
		// x * 10e-digitMove   10-3   1200000000 -> 1200000.000
		if len(s) <= digitMove {
			prefix := "0." + strings.Repeat("0", digitMove-len(s))
			return prefix + strings.TrimRight(s, "0")
		}
		lastNoZeroPos := -1
		cnt := len(s) - digitMove - 1
		for i := 0; i < len(s); i++ {
			newByte[pos] = s[i]
			pos++
			if lastNoZeroPos != -1 && s[i] != '0' {
				lastNoZeroPos = pos
			}
			if i == cnt {
				newByte[pos] = '.'
				lastNoZeroPos = pos
				pos++
			}
		}
		r := string(newByte[:lastNoZeroPos])
		return r
	}
	if d.state == 0 {
		return op(d.nano, 9)
	}
	if d.state == 1 {
		return d.nanoBigNum.Shift(-NanoOffset).String()
	}
	if d.state == 2 {
		return op(d.atto, 18)
	}
	if d.state == 3 {
		return d.attoBigNum.Shift(-2 * NanoOffset).String()
	}
	panic(errors.New(fmt.Sprintf("invalid num state: %d", d.state)))
}

func (d DynNum) MarshalJSON() ([]byte, error) {
	return []byte(d.String()), nil
}

func (d *DynNum) UnmarshalJSON(data []byte) error {
	if bytes.Equal(data, []byte{'n', 'u', 'l', 'l'}) {
		return nil
	}
	if len(data) == 0 {
		return nil
	}
	var err error
	s := strings.Trim(string(data), "\"' \t")

	if *d, err = NewDynNumFromString(s); err != nil {
		return err
	}
	return nil
}
