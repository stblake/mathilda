# FunctionInterpolation: an InterpolatingFunction from an expression
fi = FunctionInterpolation[NIntegrate[Exp[-t^2], {t, 0, x}], {x, 0, 3}]
{fi[1.5], N[Sqrt[Pi]/2 Erf[1.5]]}
fi'[1.0]
g = FunctionInterpolation[Sin[x] Exp[-x/4], {x, 0, 10}];
g[2.3] - Sin[2.3] Exp[-2.3/4]
