package numerical

import (
	"math"
	//"fmt"
)

//TODO: implement these methods for a more general `Function` interface

func BisectionSolve(p Poly, a, b float64) float64 {
	// Its upto user to ensure that the polynomial acually has a root in (a, b)
	var inc bool
	if p.Eval(a) < p.Eval(b) {
		inc = true
	}
	n := 0
	c := (a+b)/2
	for math.Abs(p.Eval(c)) > EPS {
		if p.Eval(c) > 0 {
			if inc {
				b = c
			} else {
				a = c
			}
		} else {
			if inc {
				a =c
			} else {
				b = c
			}
		}
		c = (a + b)/2
		n++
	}
	return c
}

func NewtonRaphson(p Poly, a, b float64) float64 {
	x := (a+b)/2
	n := 0
	for math.Abs(p.Eval(x)) > EPS {
		x = x - p.Eval(x)/p.Derivative().Eval(x)
		n++
	}
	return x
}
