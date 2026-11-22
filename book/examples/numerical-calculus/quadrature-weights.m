# Gauss-Legendre nodes and weights, and the error constant
GaussianQuadratureWeights[3, -1, 1]
GaussianQuadratureWeights[2, 0, 1, 20]
rule = GaussianQuadratureWeights[5, 0, 1];
Total[(#[[2]] Exp[#[[1]]]) & /@ rule] - (E - 1)
GaussianQuadratureError[3, f, -1, 1]
GaussianQuadratureError[5, Exp, 0, 1]
