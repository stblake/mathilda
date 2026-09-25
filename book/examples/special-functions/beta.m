# 4.7.1 The Euler beta function B(a,b) = Gamma[a] Gamma[b] / Gamma[a+b]
# A positive integer argument collapses the gamma ratio to a rational
Beta[2, 3]
Beta[3, 1/3]
# Half-integer arguments fold to rational multiples of Pi
Beta[1/2, 1/2]
Beta[5/2, 7/2]
# Otherwise the residual gamma form survives
Beta[1/3, 1/3]
# A fully numeric call evaluates
N[Beta[2.3, 3.2], 10]
