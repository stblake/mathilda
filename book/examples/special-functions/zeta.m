# 4.7.2 The Riemann zeta function and its relatives
# Even positive integers close into rational multiples of a power of Pi
Zeta[2]
Zeta[4]
# ... which is exactly the Basel sum, recognized directly by Sum
Sum[1/n^2, {n, 1, Infinity}]
# Special values off the critical strip, from the functional equation
Zeta[0]
Zeta[-1]
# The trivial zeros sit at the negative even integers
Zeta[-2]
# The odd values have no known closed form and stay symbolic
Zeta[3]
N[Zeta[3], 30]
# On the critical line, at the height of the first nontrivial zero, Zeta nearly vanishes
N[Zeta[1/2 + 14.134725 I], 10]
# The Hurwitz zeta generalizes the sum to shifted denominators
HurwitzZeta[2, 3]
# The Stieltjes constants are the Laurent coefficients of Zeta about s = 1
StieltjesGamma[0]
