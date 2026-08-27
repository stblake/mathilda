# Section 4.5.4: NSum -- convergence acceleration, and a method mismatch.
NSum[1/n^2, {n, 1, Infinity}]
NSum[(-1)^n/(2 n + 1), {n, 0, Infinity}, Method -> "AlternatingSigns"]
NSum[1/n^(11/10), {n, 1, Infinity}, WorkingPrecision -> 40]
NSum[1/n^2, {n, 1, Infinity}, Method -> "WynnEpsilon"]
