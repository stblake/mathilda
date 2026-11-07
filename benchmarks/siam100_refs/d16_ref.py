"""Reference A for SIAM100 App. D problem 16 (exp 191).

Heat equation u_t = u_xx on [0,3], u(0)=u(3)=0, u(x,0) = phi_2(x) (Daubechies
D4 scaling function, support [0,3]).  Find t with
    L(t) = Int_0^3 sqrt(1 + u_x(x,t)^2) dx = 5.4.
(L(t) -> infinity as t -> 0: the graph of phi_2 has infinite length; t* ~ 5.2e-7.)

Method A (this file; gmpy2 MPFR, own radix-2 FFT, no Wolfram code):
  * phihat(w) = Int phi_2(x) e^{-iwx} dx = prod_{j>=1} m0(w/2^j) (exact infinite
    product from the two-scale relation), m0(xi) = sum_k p_k e^{-ik xi},
    p = ((1+r3)/8, (3+r3)/8, (3-r3)/8, (1-r3)/8); truncated at J factors, tail
    replaced by exp(-i mu1 w/2^J) (mu1 = sum k p_k = first moment), error
    O((w/2^J)^2) < 1e-70.
  * sine coefficients b_k = (2/3) Int phi sin(k pi x/3) = -(2/3) Im phihat(k pi/3);
    u_x(x,t) = sum_k b_k w_k e^{-w_k^2 t} cos(w_k x), w_k = k pi/3, truncated
    where w_k e^{-w_k^2 t} < 1e-75.
  * u_x is an even 6-periodic cosine series, so the integrand is periodic and
    analytic: the trapezoidal rule with M panels on [0,3] converges
    geometrically; u_x on the grid by one length-2M FFT.  M is doubled until
    L(t) is stable to the working precision.
  * t from the secant method on log t.
Usage: python3.11 d16_ref.py BITS LOG2M [TARGET TGUESS]   (defaults 5.4, 5.2127825e-7)
       python3.11 d16_ref.py BITS LOG2M conv T     (trapezoid convergence check)
Results: L=5.4 -> t* = 5.21278254742981720768011478822745542811125648357811816e-7
         L=4.5 -> see truth.tsv of experiment 191 (d16.machine-reduced)
"""
import sys
import numpy as np
import gmpy2
from gmpy2 import mpfr, mpc

BITS = int(sys.argv[1]) if len(sys.argv) > 1 else 260
LOG2M = int(sys.argv[2]) if len(sys.argv) > 2 else 19
gmpy2.get_context().precision = BITS
PI = gmpy2.const_pi()
r3 = gmpy2.sqrt(mpfr(3))
p = [(1 + r3) / 8, (3 + r3) / 8, (3 - r3) / 8, (1 - r3) / 8]
mu1 = p[1] + 2 * p[2] + 3 * p[3]
CONV = len(sys.argv) > 3 and sys.argv[3] == "conv"
TARGET = "5.4" if CONV or len(sys.argv) <= 3 else sys.argv[3]
TGUESS = "5.2127825e-7" if CONV or len(sys.argv) <= 4 else sys.argv[4]
TSTAR_APPROX = float(sys.argv[4]) if CONV else float(TGUESS)
# number of modes: w^2 t > (BITS*ln2 + 20) at t = 0.9 t*
K = int(np.sqrt((BITS * 0.6931 + 25) / (0.9 * TSTAR_APPROX)) * 3 / np.pi) + 1
J = int(np.log2(K * 3.2)) + BITS // 2 + 10


def coefficients():
    ks = np.arange(1, K + 1)
    prod = np.array([mpc(1)] * K, dtype=object)
    for j in range(1, J + 1):
        base = gmpy2.exp(mpc(0, -PI / (3 * mpfr(2) ** j)))
        z = np.multiply.accumulate(np.array([base] * K, dtype=object))  # base^k
        z2 = z * z
        prod = prod * (p[0] + p[1] * z + p[2] * z2 + p[3] * z2 * z)
    tail = np.array([gmpy2.exp(mpc(0, -mu1 * int(k) * PI / (3 * mpfr(2) ** J))) for k in ks], dtype=object)
    prod = prod * tail
    b = np.array([-(mpfr(2) / 3) * c.imag for c in prod], dtype=object)
    w = np.array([int(k) * PI / 3 for k in ks], dtype=object)
    return b, w


def fft(x, roots):
    n = len(x)
    bits = n.bit_length() - 1
    rev = np.zeros(n, dtype=np.int64)
    for i in range(bits):
        rev |= ((np.arange(n) >> i) & 1) << (bits - 1 - i)
    x = x[rev]
    m = 2
    while m <= n:
        half = m // 2
        wv = roots[:: n // m][:half]
        x = x.reshape(-1, m)
        u = x[:, :half]
        t = x[:, half:] * wv
        x = np.concatenate([u + t, u - t], axis=1)
        m *= 2
    return x.reshape(-1)


def make_L(b, w, logM):
    M = 2 ** logM
    N = 2 * M
    ang = 2 * PI / N
    roots = np.array([gmpy2.exp(mpc(0, -ang * k)) for k in range(N // 2)], dtype=object)
    def L(t):
        c = b * w * np.array([gmpy2.exp(-wk * wk * t) for wk in w], dtype=object)
        a = np.array([mpc(0)] * N, dtype=object)
        for k in range(len(c)):
            a[(k + 1) % N] += c[k]
        ux = np.array([z.real for z in fft(a, roots)[: M + 1]], dtype=object)
        f = np.array([gmpy2.sqrt(1 + v * v) for v in ux], dtype=object)
        s = sum(f) - (f[0] + f[-1]) / 2
        return 3 * s / M
    return L


if __name__ == "__main__":
    b, w = coefficients()
    print("K =", K, "J =", J, "b1 =", gmpy2.mpfr(b[0], 80), flush=True)
    if CONV:
        t = mpfr(sys.argv[4])
        for lm in range(LOG2M - 2, LOG2M + 1):
            print(lm, make_L(b, w, lm)(t), flush=True)
        sys.exit()
    L = make_L(b, w, LOG2M)
    target = mpfr(TARGET)
    x0 = gmpy2.log(mpfr(TGUESS)); x1 = x0 + mpfr("1e-8")
    f0, f1 = L(gmpy2.exp(x0)) - target, L(gmpy2.exp(x1)) - target
    for it in range(12):
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0 = x1, f1
        x1 = x2
        f1 = L(gmpy2.exp(x1)) - target
        print(it, gmpy2.exp(x1), float(f1), flush=True)
        if abs(f1) < mpfr(2) ** (-BITS + 12):
            break
    print("t* =", gmpy2.exp(x1))
