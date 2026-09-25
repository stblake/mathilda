# 4.2 Sparse and multivariate views of a polynomial
MonomialList[(x + y)^3, {x, y}]
CoefficientRules[(x + y)^3, {x, y}]
FromCoefficientRules[{{3, 0} -> 1, {2, 1} -> 3, {1, 2} -> 3, {0, 3} -> 1}, {x, y}]
