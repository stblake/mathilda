"""c6 reference (method 1): expected visits to the origin
   G(eps) = (2 pi)^-2 Int_{[-pi,pi]^2} dk / (1 - phi(k)),
   phi = (1/2)cos k1 + 2 i eps sin k1 + (1/2) cos k2,
with the k2 integral done in closed form (Int_0^pi dk/(A - cos(k)/2) =
pi/sqrt((A-1/2)(A+1/2))), leaving a 1-D integral; return probability
1 - 1/G = 1/2  <=>  G = 2, solved by mpmath.findroot.  The integrand has a
k^(-1/2) endpoint singularity at k = 0, removed by k = u^2.
Also (check) the book's AGM closed form (6.1) in mpmath."""
import mpmath as mp
def G(e):
    def f(u):
        k = u*u
        A = 1 - mp.cos(k)/2 - 2j*e*mp.sin(k)
        return 2*u*mp.re(1/(mp.sqrt(A - 0.5)*mp.sqrt(A + 0.5)))
    return mp.quad(f, mp.linspace(0, mp.sqrt(mp.pi), 9))/mp.pi
for d in (80, 110):
    mp.mp.dps = d
    e = mp.findroot(lambda e: G(e) - 2, mp.mpf('0.0619'))
    print(d, mp.nstr(e, 70))
mp.mp.dps = 110
eta = mp.findroot(lambda t: mp.agm(mp.sqrt(4-(t-1)**2), mp.sqrt(4-(t+1)**2)) - 1, mp.mpf('0.99'))
print("agm", mp.nstr(mp.sqrt(1-eta**2)/4, 70))
