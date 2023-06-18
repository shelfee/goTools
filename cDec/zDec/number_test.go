package zDec

import (
	"fmt"
	"github.com/shopspring/decimal"
	"github.com/stretchr/testify/assert"
	"math"
	"math/rand"
	"reflect"
	"testing"
)

func mustBigInt(s string) decimal.Decimal {
	v, _ := decimal.NewFromString(s)
	return v
}
func mustDynNum(s string) DynNum {
	v, _ := NewDynNumFromString(s)
	return v
}

func Test_shrinkInt(t *testing.T) {
	type args struct {
		v int64
	}
	tests := []struct {
		name  string
		args  args
		want  int64
		want1 int
	}{
		{
			name:  "",
			args:  args{},
			want:  0,
			want1: 0,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, got1 := shrinkInt(tt.args.v)
			if got != tt.want {
				t.Errorf("shrinkInt() got = %v, want %v", got, tt.want)
			}
			if got1 != tt.want1 {
				t.Errorf("shrinkInt() got1 = %v, want %v", got1, tt.want1)
			}
		})
	}
}

func Test_shrinkBigInt(t *testing.T) {
	v, _ := decimal.NewFromString("112414124253161461614613451")
	r, c := shrinkBigInt(v)
	z := r.Shift(-1)
	t.Logf("v: %s, r: %d, c: %d, z: %s", v.String(), r.Round(5), c, z.String())

	//112414124253161461614613451
	//11241412425316146161461345100
	//type args struct {
	//	v decimal.Decimal
	//}
	//tests := []struct {
	//	name  string
	//	args  args
	//	want  decimal.Decimal
	//	want1 int
	//}{
	//	// TODO: Add test cases.
	//}
	//for _, tt := range tests {
	//	t.Run(tt.name, func(t *testing.T) {
	//		got, got1 := shrinkBigInt(tt.args.v)
	//		if !reflect.DeepEqual(got, tt.want) {
	//			t.Errorf("shrinkBigInt() got = %v, want %v", got, tt.want)
	//		}
	//		if got1 != tt.want1 {
	//			t.Errorf("shrinkBigInt() got1 = %v, want %v", got1, tt.want1)
	//		}
	//	})
	//}
}

func TestDynNum_ShrinkInt(t *testing.T) {

	tests := []struct {
		name   string
		fields string
		want   interface{}
		want1  int
	}{
		{
			name:   "1",
			fields: "1234567",
			want:   int64(1234567),
			want1:  0,
		},
		{
			name:   "2",
			fields: "1234567.3214",
			want:   int64(12345673214),
			want1:  -4,
		},
		{
			name:   "3",
			fields: "-1234567.3214",
			want:   int64(-12345673214),
			want1:  -4,
		},
		{
			name:   "4",
			fields: "-9223372036854775808",
			want:   mustBigInt("-9223372036854775808"),
			want1:  0,
		},
		{
			name:   "4",
			fields: "-92233720368547758084444",
			want:   mustBigInt("-92233720368547758084444"),
			want1:  0,
		},

		{
			name:   "5",
			fields: "0.00000123",
			want:   int64(123),
			want1:  -8,
		},
		{
			name:   "6",
			fields: "922337203685477580831412415235143461432542341235632632",
			want:   mustBigInt("922337203685477580831412415235143461432542341235632632"),
			want1:  0,
		},
		{
			name:   "7",
			fields: "0.00000123922337203685477580831412415235143461432542341235632632",
			want:   mustBigInt("123922337203685477580831412415235143461432542341235632632"),
			want1:  -62,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			d, _ := NewDynNumFromString(tt.fields)
			got, got1 := d.ShrinkInt()
			switch got.(type) {
			case decimal.Decimal:
				gotStr := got.(decimal.Decimal).String()
				wantStr := got.(decimal.Decimal).String()
				if !reflect.DeepEqual(gotStr, wantStr) {
					t.Errorf("ShrinkInt() got = %v, want %v", gotStr, tt.want)
				}
			case int64:
				if !reflect.DeepEqual(got.(int64), tt.want.(int64)) {
					t.Errorf("ShrinkInt() got = %v, want %v", got, tt.want)
				}
			}
			//if !reflect.DeepEqual(got, tt.want) {
			//	t.Errorf("ShrinkInt() got = %v, want %v", got, tt.want)
			//}
			if got1 != tt.want1 {
				t.Errorf("ShrinkInt() got1 = %v, want %v", got1, tt.want1)
			}
		})
	}
}

func TestDynNum_StressAdd(t *testing.T) {
	cnt := 50000
	type args struct {
		dec  decimal.Decimal
		dyn  DynNum
		s    string
		dynS string
		decS string
	}
	newArg := func(op func() float64) args {
		t := decimal.NewFromFloat(op()).Mul(decimal.NewFromFloat(op()))
		s := t.String()
		v, _ := NewDynNumFromString(s)
		return args{
			dec: t,
			s:   s,
			dyn: v,
		}
	}
	testV1 := make([]args, cnt)
	testV2 := make([]args, cnt)

	for i := 0; i < cnt; i++ {
		op := func() float64 {
			return (math.Round(rand.Float64()*100) + float64(rand.Int63n(1000000)+2000000)) / 100
		}
		op2 := func() float64 {
			return math.Round(rand.Float64()*100) / 100
		}
		//op := rand.NormFloat64
		//if rand.Float64() < 0.25 {
		//	op = rand.ExpFloat64
		//} else if rand.Float64() < 0.33 {
		//	op = rand.Float64
		//} else if rand.Float64() < 0.5 {
		//	op = func() float64 {
		//		return math.Pow(rand.Float64()*2, 100)
		//	}
		//}
		testV1[i] = newArg(op)
		testV2[i] = newArg(op2)
	}
	for j := 0; j < 100; j++ {
		for i := 0; i < cnt; i++ {
			testV1[i].decS = testV1[i].dec.Add(testV2[i].dec).String()
			testV1[i].dynS = testV1[i].dyn.Add(testV2[i].dyn).String()
			if j == 0 {
				assert.Equal(t, testV1[i].decS, testV1[i].dynS)
			}
			//if i%10 == 0 {
			//	t.Logf("%s+%s=%s", testV1[i].s, testV2[i].s, testV1[i].dynS)
			//}
		}
	}
}

func TestDynNum_Add(t *testing.T) {
	type args struct {
		n Number
	}
	tests := []struct {
		name   string
		fields string
		args   args
		want   Number
	}{
		{
			name:   "1",
			fields: "1.321351235",
			args: args{
				n: mustDynNum("3.2143252"),
			},
		},
		{
			name:   "2",
			fields: "1.321351235",
			args: args{
				n: Float64(3.21),
			},
		},
		{
			name:   "3",
			fields: "0.000000000012144124",
			args: args{
				n: mustDynNum("1441251643631164.999931231341535"),
			},
		},
		{
			name:   "4",
			fields: "4223372030",
			args: args{
				n: Float64(4223372030),
			},
		},
		{
			name:   "4",
			fields: "-4223372030",
			args: args{
				n: Float64(4223372030),
			},
		},
		{
			name:   "4",
			fields: "-4223372030",
			args: args{
				n: Float64(-49223372030),
			},
		},
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			want1, _ := decimal.NewFromString(tt.fields)
			add := mustBigInt(tt.args.n.String())
			want := want1.Add(add)
			t.Run(tt.name, func(t *testing.T) {
				d, _ := NewDynNumFromString(tt.fields)
				if got := d.Add(tt.args.n); !reflect.DeepEqual(got.String(), want.String()) {
					t.Errorf("Add(%s + %s) = %v, want %v", tt.fields, tt.args.n.String(), got, want.String())
				}
			})
		})
	}
}

func TestDynNum_String(t *testing.T) {

	tests := []struct {
		name   string
		fields string
		want   string
	}{
		// TODO: Add test cases.
		{
			name:   "1",
			fields: "1368559650",
			want:   "1368559650",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			d := mustDynNum(tt.fields)
			assert.Equalf(t, tt.want, d.String(), "String()")
		})
	}
}

func Test_tryPow(t *testing.T) {
	type args struct {
		v float64
		p float64
	}
	tests := []struct {
		name    string
		args    args
		wantR   float64
		wantMov int
	}{
		{
			name: "1",
			args: args{
				v: 10,
				p: 100000,
			},
			wantR:   1,
			wantMov: 1000000,
		},
		{
			name: "1",
			args: args{
				v: -10,
				p: 100000,
			},
			wantR:   math.NaN(),
			wantMov: 0,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotR, gotMov := tryPow(tt.args.v, tt.args.p)
			assert.Equalf(t, tt.wantR, gotR, "tryPow(%v, %v)", tt.args.v, tt.args.p)
			assert.Equalf(t, tt.wantMov, gotMov, "tryPow(%v, %v)", tt.args.v, tt.args.p)
		})
	}
}

func TestDynNum_Pow(t *testing.T) {

	var a Number
	a = Float64(0)
	b := NewDynNumFromFloat64(1)
	a = a.Add(b)
	fmt.Printf("%s", a.String())
	//type args struct {
	//	n Number
	//}
	//tests := []struct {
	//	name   string
	//	fields string
	//	args   args
	//	want   Number
	//}{
	//	{
	//		name:   "1",
	//		fields: "2",
	//		args: args{
	//			n: Float64(0.5),
	//		},
	//		want: nil,
	//	},
	//}
	//for _, tt := range tests {
	//	t.Run(tt.name, func(t *testing.T) {
	//		d := mustDynNum(tt.fields)
	//
	//		r := d.Pow(tt.args.n)
	//
	//		assert.Equalf(t, tt.want, r, "Pow(%v)", tt.args.n)
	//	})
	//}
}
