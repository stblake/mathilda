# Section 4.6.4 (classical thread): perfect numbers and Mersenne primes.
Select[Range[10000], DivisorSigma[1, #] == 2 # &]
Table[2^(p - 1) (2^p - 1), {p, {2, 3, 5, 7}}]
PrimeQ[2^13 - 1]
PrimeQ[2^11 - 1]
FactorInteger[2^11 - 1]
