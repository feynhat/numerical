package numerical

type Function interface {
	Eval(x float64) float64
}

type Differentiable interface {
	Function
	Derivative() Function
}

type Integrable interface {
	Function
	AntiDerivative() Function
}
