package numerical

import (
	//"math"
)

func LagrangeInterpolPoly(nodes, vals []float64) Poly {
	n := len(nodes)-1
	result := Poly{[]float64{0}}
	lagrangePolys := make([]Poly, n+1)
	for i := 0; i <= n; i++ {
		num := Poly{[]float64{1}}
		den := 1.0
		for j := 0; j <= n; j++ {
			if j != i {
				num = PolyMul(num, Poly{[]float64{-nodes[j], 1}})
				den *= (nodes[i] - nodes[j])
			}
		}
		lagrangePolys[i] = PolyMul(num, Poly{[]float64{1/den}})
	}
	for i := 0; i <= n; i++ {
		result = PolyAdd(result, PolyMul(Poly{[]float64{vals[i]}}, lagrangePolys[i]))
	}
	return result
}

func NewtonDivDiffTab(nodes, vals []float64) [][]float64 {
	n := len(nodes)-1
	divDiffTab := make([][]float64, n+1)
	for i := 0; i < n+1; i++ {
		divDiffTab[i] = make([]float64, n+1)
	}
	for i := 0; i <= n; i++ {
		divDiffTab[i][0] = vals[i]
	}
	for j := 1; j <= n; j++ {
		for i := 0; i <= n-j; i++ {
			divDiffTab[i][j] = (divDiffTab[i+1][j-1] - divDiffTab[i][j-1])/(nodes[i+j]-nodes[i])
		}
	}
	return divDiffTab
}

func NewtonDivDiff(nodes, vals []float64) Poly {
	n := len(nodes)-1
	divDiffTab := NewtonDivDiffTab(nodes, vals)
	multipliers := divDiffTab[0]
	result := Poly{[]float64{multipliers[0]}}
	nodePoly := Poly{[]float64{1}}
	for i := 1; i <= n; i++ {
		nodePoly = PolyMul(nodePoly, Poly{[]float64{-nodes[i-1], 1}})
		constMultiplier := Poly{[]float64{multipliers[i]}}
		result = PolyAdd(result, PolyMul(constMultiplier, nodePoly))
	}
	return result
}

func ForwardDiffTab(vals []float64) [][]float64 {
	n := len(vals)-1
	forwDiffTab := make([][]float64, n+1)
	for i := 0; i < n+1; i++ {
		forwDiffTab[i] = make([]float64, n+1)
	}
	for i := 0; i <= n; i++ {
		forwDiffTab[i][0] = vals[i]
	}
	for j := 1; j <= n; j++ {
		for i := 0; i <= n-j; i++ {
			forwDiffTab[i][j] = forwDiffTab[i+1][j-1] - forwDiffTab[i][j-1]
		}
	}
	return forwDiffTab
}

func ForwardDiffInterpolation(start, space float64, vals []float64) Poly {
	n := len(vals) - 1
	delta := ForwardDiffTab(vals)
	res := Poly{[]float64{delta[0][0]}}
	cm := Poly{[]float64{1}}
	pm := Poly{[]float64{1}}
	mu := Poly{[]float64{-start, 1}}
	fact := 1
	pow := 1.0
	for j := 1; j <= n; j++ {
		pow *= space
		fact *= j
		cm = PolyMul(Poly{[]float64{delta[0][j]}}, Poly{[]float64{1.0/(float64(fact)*pow)}})
		mu = Poly{[]float64{-start - space*float64(j-1), 1}}
		pm = PolyMul(pm, mu)
		res = PolyAdd(res, PolyMul(cm, pm))
	}
	return res
}
