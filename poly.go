package numerical

import (
	"strconv"
)

type Poly struct {
	Coeffs []float64
}

func NewPoly (deg int) Poly {
	return Poly{make([]float64, deg+1)}
}

func (p Poly) Str() string {
	d := p.Deg()

	if d == -1 {
		return "0"
	}

	s := ""
	if p.Coeffs[d] == -1 {
		if d == 0 {
			s += "-1"
		} else {
			s += "-"
		}
	} else if p.Coeffs[d] != 1 || d == 0 {
		s += strconv.FormatFloat(p.Coeffs[d], 'g', -1, 64)
	}
	if d == 1 {
		s += "x"
	} else if d > 1 {
		s += "x^" + strconv.Itoa(d)
	}
	for i := d-1; i >= 0; i-- {
		if p.Coeffs[i] > 0 {
			s +=  " + "
			if p.Coeffs[i] != 1 || i == 0 {
				s += strconv.FormatFloat(p.Coeffs[i], 'g', -1, 64)
			}
		} else if p.Coeffs[i] < 0 {
			s += " - "
			if p.Coeffs[i] != -1 || i == 0 {
				s += strconv.FormatFloat(-p.Coeffs[i], 'g', -1, 64)
			}
		}
		if p.Coeffs[i] != 0 {
			if i == 1 {
				s += "x"
			} else if i > 1 {
				s += "x^" + strconv.Itoa(i)
			}
		}
	}

	return s
}

func (p *Poly) Deg() int {
	i := len(p.Coeffs) - 1
	for ; i >= 0 && p.Coeffs[i] == 0; i-=1 {
	}
	return i
}

func (p Poly) Eval(x float64) float64 {
	val := p.Coeffs[0]
	pow := 1.0
	for i:=1; i <= p.Deg(); i+=1 {
		pow *= x
		val += p.Coeffs[i]*pow
	}
	return val
}

func PolyAdd(p, q Poly) Poly {
	degp := p.Deg()
	degq := q.Deg()
	var degr int
	if p.Deg() > q.Deg() {
		degr = p.Deg()
	} else {
		degr = q.Deg()
	}
	r := NewPoly(degr)
	for i := 0; i <= degr; i+=1 {
		if i <= degp {
			r.Coeffs[i] = p.Coeffs[i]
		}
		if i <= degq {
			r.Coeffs[i] += q.Coeffs[i]
		}
	}
	return r
}

func PolySub(p, q Poly) Poly {
	degp := p.Deg()
	degq := q.Deg()
	big := degp
	if degp < degq {
		big = degq
	}
	r := Poly{make([]float64, big+1)}
	for i := 0; i <= big; i++ {
		if i <= degp {
			r.Coeffs[i] = p.Coeffs[i]
		}
		if i <= degq {
			r.Coeffs[i] -= q.Coeffs[i]
		}
	}
	return r
}

func PolyMul(p, q Poly) Poly {
	m, n := p.Deg(), q.Deg()
	r := Poly{make([]float64, m+n+1)}
	for i := 0; i <= m; i++ {
		for j := 0; j <= n; j++ {
			r.Coeffs[i+j] += p.Coeffs[i]*q.Coeffs[j]
		}
	}
	return r
}

func PolyDiv(a, b Poly) (Poly, Poly) {
	var q Poly
	for b.Deg() <= a.Deg() {
		dega := a.Deg()
		degb := b.Deg()
		qc := a.Coeffs[dega]/b.Coeffs[degb]
		qp := dega - degb
		qi := Poly{make([]float64, qp+1)}
		qi.Coeffs[qp] = qc
		a = PolySub(a, PolyMul(qi, b))
		q = PolyAdd(q, qi)
	}
	return q, a
}

func (p Poly) Derivative() Poly {
	n := p.Deg()
	res := NewPoly(n-1)
	for i := 1; i <= n; i++ {
		res.Coeffs[i-1] = float64(i)*p.Coeffs[i]
	}
	return res
}
