"""bc6.5 second method (the I0 integral of bc6_ref.py converges only like
t^(-3/2) and reaches ~40 digits): Watson's integral as a 1-D integral of the
2-D square-lattice Green's function, G = (1/pi) Int_0^pi (2/(pi A)) K(m = 4/(9 A^2)) dk,
A = 1 - cos(k)/3 (mpmath.ellipk; k = u^2), against the Glasser-Zucker Gamma form."""
import mpmath as mp
from mpmath import mpf, pi, sqrt, cos, gamma, quad, nstr
mp.mp.dps = 80
def f(u):
    k = u*u
    with mp.extradps(60):
        A = 1 - cos(k)/3
        return +(2*u*2/(pi*A)*mp.ellipk(4/(9*A*A)))
G = quad(f, mp.linspace(0, sqrt(pi), 17))/pi
Gz = sqrt(6)/(32*pi**3)*gamma(mpf(1)/24)*gamma(mpf(5)/24)*gamma(mpf(7)/24)*gamma(mpf(11)/24)
print('bc6.5', nstr(1 - 1/Gz, 70)); print('     ', nstr(1 - 1/G, 70), ' agree to', nstr(-mp.log10(abs(G-Gz)/G), 3), 'digits')
