# 4.4.6 Exact reduction and the ZeroTest option
# The trap: this identity is TRUE but does not auto-simplify to 1.
Cos[x]^2 + Sin[x]^2
# Two rows that are mathematically equal (top-left simplifies to 1).
mat = {{Cos[x]^2 + Sin[x]^2, 1}, {1, 1}}
NullSpace[mat]
NullSpace[mat, ZeroTest -> (Simplify[#] === 0 &)]
# The numeric analogue: tiny pivots kept exactly, lost at machine precision.
big = {{1, 1, 1}, {0, 10^-10, 0}, {0, 0, 10^-20}};
MatrixRank[big]
MatrixRank[N[big]]
MatrixRank[N[big], Tolerance -> 0]
