"""Reference for holdout-B family d15 (experiment 250): blow-up (or quenching)
time of a reaction-diffusion PDE or an ODE, i.e. where the solution's domain
ends.

Every problem is built to have an exact reduction:
  * spatially constant data with Neumann conditions -> the ODE u' = f(u),
    t* = Int du/f(u);
  * the substitution w = e^-u (resp. w = 1/u) turns u_t = Lap u - |grad u|^2 + e^u
    (resp. u_t = u_xx - 2 u_x^2/u + u^2) into the LINEAR w_t = Lap w - 1, and
    blow-up of u is the first zero of w (eigenfunction series / Laplace transform);
  * Riccati / Hamiltonian ODEs with Bessel / Beta-function blow-up times.
Two independent methods per problem (A, B).
Usage: bd15_ref.py k [dps]
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import (mp, mpf, log, exp, quad, inf, lambertw, findroot, nstr, pi, sqrt, besselj, besseljzero,
                    besseli, beta, hyp2f1, sech, odefun, invertlaplace, cos, nsum)

k = int(sys.argv[1]); mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 70
out = {}
t0 = time.time()
if k == 1:   # u_t = u_xx + u^2 + u, Neumann, u0 = 1
    out = {"A log 2": log(2), "B quad Int_1^oo du/(u^2+u)": quad(lambda u: 1 / (u * u + u), [1, inf])}
elif k == 2:  # u_t = u_xx - u_x^2 + e^u on [0,pi], u = 0 on the boundary and at t = 0
    def phiA(t):   # centre of phi_t = phi_xx + 1 (sine series)
        return pi ** 2 / 8 - nsum(lambda j: (-1) ** j * 4 / (pi * (2 * j + 1) ** 3) * exp(-(2 * j + 1) ** 2 * t), [0, inf])
    A = findroot(lambda t: phiA(t) - 1, mpf("1.0"))
    phiB = lambda t: invertlaplace(lambda s: (1 - sech(sqrt(s) * pi / 2)) / s ** 2, t, method="talbot")
    B = findroot(lambda t: phiB(t) - 1, A + mpf(10) ** -6)
    out = {"A sine series of w = e^-u (w_t = w_xx - 1), w(pi/2, t) = 0": A,
           "B Talbot inversion of (1 - sech(pi sqrt(s)/2))/s^2": B}
elif k == 3:  # 2-D [0,pi]^2 Neumann, u0 = -log(2 + cos x cos y)
    A = 2 + lambertw(-2 * exp(-4)).real / 2
    B = findroot(lambda t: 2 - t - exp(-2 * t), mpf("1.98"))
    out = {"A 2 + W0(-2 e^-4)/2": A, "B findroot of min w = 2 - t - e^(-2t)": B}
elif k == 4:  # u_t = u_xx - 2 u_x^2/u + u^2, Neumann [0,pi], u0 = 1/(2 + cos x)
    A = 2 + lambertw(-exp(-2)).real
    B = findroot(lambda t: 2 - t - exp(-t), mpf("1.84"))
    out = {"A 2 + W0(-e^-2)": A, "B findroot of min w = 2 - t - e^-t": B}
elif k == 5:  # quenching u_t = u_xx + e^u/(1-u), Neumann, u0 = 0
    out = {"A 1/e": 1 / exp(1), "B quad Int_0^1 (1-u) e^-u du": quad(lambda u: (1 - u) * exp(-u), [0, 1])}
elif k == 6:  # y' = y^2 + t^2, y(0) = 0
    z = findroot(lambda z: besselj(-mpf(1) / 4, z), mpf(2.0))
    A = sqrt(2 * z)
    f = odefun(lambda t, w: [w[1], -t * t * w[0]], 0, [mpf(1), mpf(0)])
    B = findroot(lambda t: f(t)[0], mpf(2.0))
    out = {"A sqrt(2 j), j = first zero of J_{-1/4}": A, "B odefun w'' + t^2 w = 0, first zero of w": B}
elif k == 7:  # y'' = 6 y^2, y(0) = 1, y'(0) = 0
    out = {"A B(1/6, 1/2)/6": beta(mpf(1) / 6, mpf(1) / 2) / 6,
           "B quad Int_1^oo dy/(2 sqrt(y^3-1)), y = 1 + s^2": quad(lambda s: 1 / sqrt(3 + 3 * s * s + s ** 4), [0, 1, inf])}
elif k == 8:  # u_t = u_xx + v^2, v_t = v_xx + u^2, Neumann, u0 = 1, v0 = 1/2
    c = mpf(7) / 8
    out = {"A quad Int_1^oo (u^3 - 7/8)^(-2/3) du": quad(lambda u: (u ** 3 - c) ** (-mpf(2) / 3), [1, 2, inf]),
           "B 2F1(1/3, 2/3; 4/3; 7/8)": hyp2f1(mpf(1) / 3, mpf(2) / 3, mpf(4) / 3, c)}
elif k == 9:  # disk radius 3, u_t = Lap u - |grad u|^2 + e^u, u = 0 on r = 3 and at t = 0
    R = 3
    def phiA(t):
        return mpf(R * R) / 4 - nsum(lambda n: 2 * R * R / (besseljzero(0, int(n)) ** 3 * besselj(1, besseljzero(0, int(n))))
                                     * exp(-besseljzero(0, int(n)) ** 2 * t / R ** 2), [1, inf])
    A = findroot(lambda t: phiA(t) - 1, mpf("0.5"))
    phiB = lambda t: invertlaplace(lambda s: (1 - 1 / besseli(0, R * sqrt(s))) / s ** 2, t, method="talbot")
    B = findroot(lambda t: phiB(t) - 1, A + mpf(10) ** -6)
    out = {"A Fourier-Bessel series of w = e^-u, w(0, t) = 0": A, "B Talbot inversion of (1 - 1/I0(3 sqrt s))/s^2": B}
elif k == 10:  # square [0,4]^2, u_t = Lap u - |grad u|^2 + e^u, u = 0 on boundary and at t = 0
    L = 4
    sg = lambda m: (-1) ** ((m - 1) // 2)
    S = mpf(L * L) / 8 - nsum(lambda j: (-1) ** j * 4 * L * L / (pi * (2 * j + 1)) ** 3 * sech((2 * j + 1) * pi / 2), [0, inf])
    def phiA(t):   # steady (sech series) minus the transient double eigen-series
        tr = mpf(0)
        for m in range(1, 80, 2):
            for n in range(1, 80, 2):
                lam = (pi / L) ** 2 * (m * m + n * n)
                term = 16 / (pi ** 2 * m * n * lam) * exp(-lam * t)
                if term < mpf(10) ** (-mp.dps - 5):
                    break
                tr += sg(m) * sg(n) * term
        return S - tr
    A = findroot(lambda t: phiA(t) - 1, mpf("1.5"))
    def Fhat(s):
        rs = sqrt(s)
        base = (1 - sech(rs * L / 2)) / s
        corr = mpf(0)
        for m in range(1, 400, 2):
            q = sqrt(s + (m * pi / L) ** 2)
            term = sg(m) * 4 / (m * pi) / (s + (m * pi / L) ** 2) * sech(q * L / 2)
            corr += term
            if abs(term) < mpf(10) ** (-mp.dps - 5):
                break
        return (base - corr) / s
    B = findroot(lambda t: invertlaplace(Fhat, t, method="talbot") - 1, A + mpf(10) ** -6)
    out = {"A steady sech series - transient double sine series": A,
           "B Talbot inversion of the x-sine / y-cosh Laplace-domain series": B}
for kk, vv in out.items():
    print(kk); print("   ", nstr(vv, 70))
print("secs", round(time.time() - t0))
