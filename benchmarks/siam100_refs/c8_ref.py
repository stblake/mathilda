"""C8 reference: centre of [-1,1]^2 reaches u=1 (one side held at 5).
u(0,0,t) = (5/4)(1 - w(t)^2), w = centre value of the 1-D problem on [-1,1]
(zero ends, initial 1).  Solve w(t) = 1/sqrt(5).
Method 1: eigen series  w = (4/pi) sum (-1)^k/(2k+1) exp(-(2k+1)^2 pi^2 t/4).
Method 2: image (theta-transformed) series w = 1 - 2 sum (-1)^n erfc((2n+1)/(2 sqrt t))."""
from mpmath import mp, mpf, sqrt, pi, exp, erfc, nsum, inf, findroot, nstr
mp.dps = 80
w1 = lambda t: 4/pi*nsum(lambda k: (-1)**int(k)/(2*k+1)*exp(-(2*k+1)**2*pi**2*t/4), [0, inf])
w2 = lambda t: 1 - 2*nsum(lambda n: (-1)**int(n)*erfc((2*n+1)/(2*sqrt(t))), [0, inf])
t1 = findroot(lambda t: w1(t) - 1/sqrt(5), mpf('0.42'))
t2 = findroot(lambda t: w2(t) - 1/sqrt(5), mpf('0.42'))
print(nstr(t1, 70)); print(nstr(t2, 70)); print("rel diff", nstr(abs(t1-t2)/t1, 5))
