package numerical

//import "fmt"

func DefiniteIntegral(fn Integrable, a, b float64) float64 {
	return fn.AntiDerivative().Eval(b) - fn.AntiDerivative().Eval(a)
}

func NewtonCotes(fn Function, a, b float64, n int) float64 {
	h := (b - a)/float64(n)
	mu := Poly{[]float64{0, 1}}
	res := 0.0
	for i := 0; i <= n; i++ {
		integrand := Poly{[]float64{1}}
		f1, f2 := 1, 1
		for j := 0; j <= n; j++ {
			if j > 0 && j <= i {
				f1 *= j
			}
			if j < n && j >= i {
				f2 *= (n-j)
			}
			if j != i {
				integrand = PolyMul(integrand, PolySub(mu, Poly{[]float64{float64(j)}}))
			}
		}
		I := DefiniteIntegral(integrand, 0, float64(n))
		w := h/(float64(f1)*float64(f2)) * I
		if (n-i)%2 == 1 {
			w = -w
		}
		res += w*fn.Eval(a + float64(i)*h)
	}
	return res
}
