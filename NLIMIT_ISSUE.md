The following is a bug in NLimit: 

In[6]:= NLimit[(Cos[x^2]/(x^2) - Cos[(x + 1)^2]/((x + 1)^2))/(1/(x^3)), x -> Infinity] // Timing
Out[6]= {0.00055, -5.28256}

We should numerically detect diverging oscillatory limits and report an informative error message. We should have extensive unit tests for this feature and not introduce any memory leaks. 

Timing[NLimit[x Sin[x], x -> Infinity]]
NLimit::noise: Cannot recognize a limiting value.  This may be due to  noise resulting from roundoff errors in which case higher WorkingPrecision,  fewer Terms, or a different Scale might help.
{0.0026, NLimit[x Sin[x], x -> \[Infinity]]}
Timing[NLimit[Sin[x] Sin[x^2], x -> Infinity]]
NLimit::noise: Cannot recognize a limiting value.  This may be due to  noise resulting from roundoff errors in which case higher WorkingPrecision,  fewer Terms, or a different Scale might help.
{0.00252, NLimit[Sin[x] Sin[x^2], x -> \[Infinity]]}
Timing[NLimit[Sin[1/x]/x, x -> 0]]
NLimit::noise: Cannot recognize a limiting value.  This may be due to  noise resulting from roundoff errors in which case higher WorkingPrecision,  fewer Terms, or a different Scale might help.
{0.001772, NLimit[Sin[1/x]/x, x -> 0]}
