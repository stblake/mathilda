# MeijerGReduce: rewrite an expression as inactive Meijer G functions of x
MeijerGReduce[Exp[-x], x]
MeijerGReduce[Sin[x], x]
MeijerGReduce[x^2 BesselJ[0, x] + Log[1 + a x], x]
g = MeijerGReduce[ArcTan[x], x]
{Activate[g] /. x -> 0.7, ArcTan[0.7]}
