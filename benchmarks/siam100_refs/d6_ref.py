# D6 reference: Int_{-1}^{1} exp(x + sin(exp(exp(x+1/3)))) dx
# Method A: mpmath tanh-sinh, interval split at the points where
# exp(exp(x+1/3)) = k*pi (the zeros of the sine), dps 80.
import mpmath as mp
mp.mp.dps = 80
f = lambda x: mp.exp(x + mp.sin(mp.exp(mp.exp(x + mp.mpf(1)/3))))
lo = mp.exp(mp.exp(mp.mpf(-2)/3)); hi = mp.exp(mp.exp(mp.mpf(4)/3))
k0 = int(mp.ceil(lo/mp.pi)); k1 = int(mp.floor(hi/mp.pi))
pts = [mp.mpf(-1)] + [mp.log(mp.log(k*mp.pi)) - mp.mpf(1)/3 for k in range(k0, k1+1)] + [mp.mpf(1)]
A = mp.quad(f, pts)
# Method B: Gauss-Legendre on the same split, dps 80
B = mp.quad(f, pts, method='gauss-legendre')
print(mp.nstr(A, 70)); print(mp.nstr(B, 70)); print(mp.nstr(A-B, 5))
