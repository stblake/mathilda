# One integer type, several machine representations, chosen automatically.
Head[3]
Head[3.0]
Head[1/3]
Head[2 + 3 I]
# Exact stays exact; one machine real makes the whole result inexact:
1/2 + 1/3
1/2 + 0.5
# N gives machine reals by default, or MPFR arithmetic at higher precision:
N[Pi]
N[Pi, 40]
Precision[N[Pi]]
Precision[N[Pi, 40]]
