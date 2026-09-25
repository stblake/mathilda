# 4.7.1 The gamma function: exact values, poles, the incomplete form
# Gamma interpolates the factorial: Gamma[n] = (n-1)!
Table[Gamma[n], {n, 1, 6}]
# Half-integers close into rational multiples of Sqrt[Pi]
Gamma[1/2]
Gamma[7/2]
Gamma[-1/2]
# The non-positive integers are poles
Gamma[0]
# The upper incomplete gamma at a positive integer first argument is elementary
Gamma[3, x]
Gamma[5, 2]
# Complex numeric evaluation to 30 digits (Spouge's approximation)
N[Gamma[I], 30]
# The derivative closes through the digamma function
D[Gamma[z], z]
