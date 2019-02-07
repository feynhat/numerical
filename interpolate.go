package numerical

import (
	//"math"
)

func NewtonDivDiffTab(nodes, vals []float64) [][]float64 {
	n := len(vals)-1
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
		nodePoly = Multiply(nodePoly, Poly{[]float64{-nodes[i-1], 1}})
		constMultiplier := Poly{[]float64{multipliers[i]}}
		result = Add(result, Multiply(constMultiplier, nodePoly))
	}
	return result
}
