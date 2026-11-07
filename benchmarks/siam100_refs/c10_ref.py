"""C10 reference: P(hit ends) from the centre of a 10x1 rectangle.
Method 1: the book's closed form (2/pi) arcsin(k), k from the singular moduli.
Method 2: the separated series (4/pi) sum (-1)^k/(2k+1) sech(5(2k+1)pi)."""
from mpmath import mp, mpf, sqrt, asin, pi, sech, nsum, inf, nstr, root
mp.dps = 80
k = (3 - 2*sqrt(2))**2 * (2 + sqrt(5))**2 * (sqrt(10) - 3)**2 * (root(5, 4) - sqrt(2))**4
a = 2/pi*asin(k)
b = 4/pi*nsum(lambda j: (-1)**int(j)/(2*j+1)*sech(5*(2*j+1)*pi), [0, inf])
print(nstr(a, 70)); print(nstr(b, 70)); print("rel diff", nstr(abs(a-b)/a, 5))
