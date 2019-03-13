package numerical

/*
	This file implements the following family of orthogonal polynomials:
	 - Legendre
	 - Chebyshev (First kind)
	 - Chebyshev (Second kind)
	 - Laguerre
	 - Hermite
	The coefficients of these polynomials are float64's (like any other polynomial in this package)
	so don't expect the .Str() of these polynomials to be very readable and thus,
	this implementation will not be very useful if you want to know what these polynomials are just for the heck of it.
	But, these are useful in doing acutal computation like least square approximation and numerical integration. 
*/
func Legendre(n int) Poly {
	var p0, p1, p2 Poly
	p0 = Poly{[]float64{1}}
	p1 = Poly{[]float64{0, 1}}
	if n == 0 {
		return p0
	} else if n == 1{
		return p1
	}
	for i:=1; i<n; i++ {
		A := PolyMul(Poly{[]float64{0, float64(2*i+1)/float64(i+1)}}, p1)
		B := PolyMul(Poly{[]float64{float64(-i)/float64(i+1)}}, p0)
		p2 = PolyAdd(A, B)
		p0 = p1
		p1 = p2
	}
	return p2
}

func Chebyshev1(n int) Poly {
	var p0, p1, p2 Poly
	p0 = Poly{[]float64{1}}
	p1 = Poly{[]float64{0, 1}}
	if n == 0 {
		return p0
	} else if n == 1{
		return p1
	}
	for i:=1; i<n; i++ {
		A := PolyMul(Poly{[]float64{0, 2}}, p1)
		B := PolyMul(Poly{[]float64{-1}}, p0)
		p2 = PolyAdd(A, B)
		p0 = p1
		p1 = p2
	}
	return p2
}
