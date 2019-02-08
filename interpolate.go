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
