# 4.4.9 The machine-exact divide: exact rational, machine LAPACK, MPFR, packed
Det[HilbertMatrix[5]]
Det[N[HilbertMatrix[5]]]
Det[N[HilbertMatrix[7], 40]]
Inverse[HilbertMatrix[3]]
a = NDArray[{{1., 2.}, {3., 4.}}]
a . a
Det[a]
