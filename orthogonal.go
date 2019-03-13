package numerical

/*
	This file implements the following families of orthogonal polynomials (via their triple-recursion relations):
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
func Legendre(n int) []Poly {
	P := make([]Poly, n+1)
	P[0] = Poly{[]float64{1}}
	P[1] = Poly{[]float64{0, 1}}
	for i:=1; i<n; i++ {
		j := float64(i)
		A := PolyMul(Poly{[]float64{0, (2*j+1)/(j+1)}}, P[i])
		B := PolyMul(Poly{[]float64{j/(j+1)}}, P[i-1])
		P[i+1] = PolySub(A, B)
	}
	return P
}

func Chebyshev1(n int) []Poly {
	T := make([]Poly, n+1)
	T[0] = Poly{[]float64{1}}
	T[1] = Poly{[]float64{0, 1}}
	for i:=1; i<n; i++ {
		A := PolyMul(Poly{[]float64{0, 2}}, T[i])
		T[i+1] = PolySub(A, T[i-1])
	}
	return T
}

func Chebyshev2(n int) []Poly {
	U := make([]Poly, n+1)
	U[0] = Poly{[]float64{1}}
	U[1] = Poly{[]float64{0, 2}}
	for i:=1; i<n; i++ {
		A := PolyMul(Poly{[]float64{0, 2}}, U[i])
		U[i+1] = PolySub(A, U[i-1])
	}
	return U
}

func Laguerre(n int) []Poly {
	L := make([]Poly, n+1)
	L[0] = Poly{[]float64{1}}
	L[1] = Poly{[]float64{1, -1}}
	for i:=1; i<n; i++ {
		j := float64(i)
		A := PolyMul(Poly{[]float64{(2*j+1)/(j+1), -1/(j+1)}}, L[i])
		B := PolyMul(Poly{[]float64{j/(j+1)}}, L[i-1])
		L[i+1] = PolySub(A, B)
	}
	return L
}

func Hermite(n int) []Poly {
	H := make([]Poly, n+1)
	H[0] = Poly{[]float64{1}}
	H[1] = Poly{[]float64{0, 2}}
	for i:=1; i<n; i++ {
		j := float64(i)
		A := PolyMul(Poly{[]float64{0, 2}}, H[i])
		B := PolyMul(Poly{[]float64{2*j}}, H[i-1])
		H[i+1] = PolySub(A, B)
	}
	return H
}
