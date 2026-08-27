# Section 4.5.3: NIntegrate basics -- infinite ranges, singularities, several dimensions.
NIntegrate[Exp[-x^2], {x, 0, Infinity}]
NIntegrate[Sin[x]/x, {x, 0, Infinity}]
NIntegrate[1/Sqrt[x], {x, 0, 1}]
NIntegrate[Exp[-x^2 - y^2], {x, -Infinity, Infinity}, {y, -Infinity, Infinity}]
