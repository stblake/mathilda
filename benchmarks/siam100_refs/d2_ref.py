# D2 reference: S = sum_{n>=1} 1/f(n), f(1)=1, f(2)=2, f(n)=n f(bitlength(n)).
# Group by bit length k: block k (n in [2^(k-1), 2^k-1]) contributes
#   B(k)/f(k),  B(k) = psi(2^k) - psi(2^(k-1)) = H(2^k-1) - H(2^(k-1)-1),
# for k >= 3 (f(n) = n f(k)); k = 1, 2 are the special values 1, 1/2 + 1/6.
# Since sum_{k>=3} 1/f(k) = S - 3/2, writing B(k) = log 2 + (B(k) - log 2):
#   S = 5/3 + log2 (S - 3/2) + E,   E = sum_{k>=3} (B(k) - log 2)/f(k),
#   S = (5/3 - 3/2 log 2 + E)/(1 - log 2),
# and B(k) - log 2 = O(2^-k), so E converges geometrically.
import mpmath as mp
mp.mp.dps = 90
def f(n):
    if n <= 2: return n
    return n * f(n.bit_length())
E = mp.mpf(0)
for k in range(3, 330):
    B = mp.digamma(mp.mpf(2)**k) - mp.digamma(mp.mpf(2)**(k-1))
    E += (B - mp.log(2)) / f(k)
S = (mp.mpf(5)/3 - mp.mpf(3)/2*mp.log(2) + E) / (1 - mp.log(2))
print(mp.nstr(S, 75))
