#!/usr/bin/env python3
"""Independent check of every hand-written reference value in corpus.m.

Each entry recomputes the series / product directly with mpmath (nsum /
nprod or an explicit partial sum) and compares it with the reference
expression the corpus states, both at 30 digits.  Nothing is taken from
Mathilda.  Run:  python3 -I tests/sum_product_oracle/verify_refs.py
"""
from mpmath import mp, mpf, nsum, nprod, inf, pi, zeta, log, e, coth, tanh, \
    sqrt, sinh, sin, cos, exp, besseli, polylog, binomial, factorial, harmonic

mp.dps = 40
X = mpf(1) / 3

# label -> (direct computation, reference value as stated in corpus.m)
REFS = {
    "inf-basel":         (lambda: nsum(lambda k: 1/k**2, [1, inf]), pi**2/6),
    "inf-zeta3":         (lambda: nsum(lambda k: 1/k**3, [1, inf]), zeta(3)),
    "inf-alt-log2":      (lambda: nsum(lambda k: (-1)**k/k, [1, inf]), -log(2)),
    "inf-k2p1":          (lambda: nsum(lambda k: 1/(k**2+1), [1, inf]), (pi*coth(pi)-1)/2),
    "inf-k2p1-from0":    (lambda: nsum(lambda k: 1/(k**2+1), [0, inf]), (pi*coth(pi)+1)/2),
    "inf-k2m1":          (lambda: nsum(lambda k: 1/(k**2-1), [2, inf]), mpf(3)/4),
    "inf-k2m4":          (lambda: nsum(lambda k: 1/(k**2-4), [3, inf]), mpf(25)/48),
    "inf-exp":           (lambda: nsum(lambda k: X**k/factorial(k), [0, inf]), exp(X)),
    "inf-log":           (lambda: nsum(lambda k: X**k/k, [1, inf]), log(mpf(3)/2)),
    "inf-k2-xk":         (lambda: nsum(lambda k: k**2*X**k, [0, inf]), mpf(3)/2),
    "inf-geom-sym":      (lambda: nsum(lambda k: X**k, [0, inf]), mpf(3)/2),
    "inf-odd-sq":        (lambda: nsum(lambda k: 1/(2*k+1)**2, [0, inf]), pi**2/8),
    "inf-alt-cube":      (lambda: nsum(lambda k: (-1)**k/k**3, [1, inf]), -3*zeta(3)/4),
    "inf-tele3":         (lambda: nsum(lambda k: 1/(k*(k+1)*(k+2)), [1, inf]), mpf(1)/4),
    "inf-central-binom": (lambda: nsum(lambda k: 1/binomial(2*k, k), [1, inf]),
                          mpf(1)/3 + 2*sqrt(3)*pi/27),
    # Slow (log) tail: nsum's default acceleration stalls, use Euler-Maclaurin.
    "inf-euler-h-k2":    (lambda: nsum(lambda k: harmonic(k)/k**2, [1, inf],
                                       method="euler-maclaurin"), 2*zeta(3)),
    "inf-sin-k-over-k":  (lambda: nsum(lambda k: sin(k)/k, [1, inf]), (pi-1)/2),
    "inf-cos-k-over-k2": (lambda: nsum(lambda k: cos(k)/k**2, [1, inf]),
                          pi**2/6 - pi/2 + mpf(1)/4),
    "inf-k4k2p1":        (lambda: nsum(lambda k: 1/(k**4+k**2+1), [1, inf]),
                          -mpf(1)/2 + pi*tanh(sqrt(3)*pi/2)/(2*sqrt(3))),
    "inf-k-over-k4p4":   (lambda: nsum(lambda k: k/(k**4+4), [1, inf]), mpf(3)/8),
    "inf-bessel":        (lambda: nsum(lambda k: 1/factorial(k)**2, [0, inf]), besseli(0, 2)),
    "inf-central-gf":    (lambda: nsum(lambda k: binomial(2*k, k)/mpf(8)**k, [0, inf]), sqrt(2)),
    "inf-4k2m1":         (lambda: nsum(lambda k: 1/(4*k**2-1), [1, inf]), mpf(1)/2),
    # Interior zeros of the term fool nsum into stopping early (the same trap
    # as the Mathilda bug this record guards), so sum the factorial series
    # explicitly: 80 terms leave a tail below 1e-100.
    "inf-interior-zero": (lambda: sum((k-3)**2/factorial(k) for k in range(80)), 5*e),
    "inf-interior-zero2": (lambda: sum((k-2)*(k-5)/factorial(k) for k in range(80)), 5*e),
    "inf-zeta-s":        (lambda: nsum(lambda k: 1/k**3, [1, inf]), zeta(3)),
    "inf-zeta-3half":    (lambda: nsum(lambda k: 1/k**mpf(1.5), [1, inf],
                                       method="euler-maclaurin"), zeta(mpf(3)/2)),
    "inf-eta-half":      (lambda: nsum(lambda k: (-1)**k/sqrt(k), [1, inf]),
                          (sqrt(2)-1)*zeta(mpf(1)/2)),
    "inf-eta-s":         (lambda: nsum(lambda k: (-1)**k/k**3, [1, inf]), -3*zeta(3)/4),
    "inf-polylog2":      (lambda: nsum(lambda k: X**k/k**2, [1, inf]), polylog(2, X)),
    "inf-log-sinh":      (lambda: nsum(lambda k: log(1+1/k**2), [1, inf]), log(sinh(pi)/pi)),
    "inf-double-ij3":    (lambda: nsum(lambda m: (m-1)/m**3, [2, inf]), zeta(2)-zeta(3)),
    "inf-step2-basel":   (lambda: nsum(lambda j: 1/(2*j+1)**2, [0, inf]), pi**2/8),
    "inf-step2-geom":    (lambda: nsum(lambda j: mpf(1)/4**j, [0, inf]), mpf(4)/3),
    "inf-deg6":          (lambda: nsum(lambda k: 1/(k**6+1), [0, inf]), None),
    "inf-deg7":          (lambda: nsum(lambda k: 1/(k**7+1), [0, inf]), None),
    "inf-deg7-3":        (lambda: nsum(lambda k: 1/(k**7+3), [0, inf]), None),
    "rng-neg-basel":     (lambda: nsum(lambda k: 1/k**2, [-inf, -1]), pi**2/6),
    "rng-2side-k2p1":    (lambda: nsum(lambda k: 1/(k**2+1), [-inf, inf]), pi*coth(pi)),
    "rng-2side-half":    (lambda: nsum(lambda k: 1/(k+mpf(1)/2)**2, [-inf, inf]), pi**2),
    "rng-2side-alt-half": (lambda: nsum(lambda k: (-1)**k/(k+mpf(1)/2), [0, inf])
                                  + nsum(lambda k: (-1)**k/(k+mpf(1)/2), [-inf, -1]), pi),
    "rng-2side-abs":     (lambda: 1 + 2*nsum(lambda k: mpf(1)/2**k, [1, inf]), mpf(3)),
    "rng-neg-geom":      (lambda: nsum(lambda k: mpf(3)**k, [-inf, 0]), mpf(3)/2),
    "rng-neg-geom2":     (lambda: nsum(lambda k: mpf(2)**k, [-inf, 0]), mpf(2)),
    "iprod-sinh":        (lambda: nprod(lambda k: 1+1/k**2, [1, inf]), sinh(pi)/pi),
    "iprod-sin":         (lambda: nprod(lambda k: 1-(mpf(1)/2)**2/(k**2*pi**2), [1, inf]),
                          2*sin(mpf(1)/2)),
    "iprod-exp":         (lambda: exp(nsum(lambda k: 1/k**2, [1, inf])), exp(pi**2/6)),
    "iprod-cantor":      (lambda: nprod(lambda k: 1+X**(2**int(k)), [0, 12]), mpf(3)/2),
    "iprod-1m1k2":       (lambda: nprod(lambda k: 1-1/k**2, [2, inf]), mpf(1)/2),
    "iprod-wallis":      (lambda: nprod(lambda k: 1-1/(4*k**2), [1, inf]), 2/pi),
    "iprod-zero-skip":   (lambda: nprod(lambda k: 1-4/k**2, [3, inf]), mpf(1)/6),
}


def primes_below(n):
    """Sieve of Eratosthenes (no sympy dependency)."""
    sieve = bytearray([1]) * n
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(range(i*i, n, i)))
    return [i for i in range(n) if sieve[i]]


def main():
    bad = 0
    for label, (direct, ref) in REFS.items():
        d = direct()
        if ref is None:
            print(f"{label:22s} value {mp.nstr(d, 32)}")
            continue
        err = abs(d - ref) / (1 + abs(ref))
        ok = err < mpf(10)**-20
        bad += not ok
        print(f"{label:22s} {'ok ' if ok else 'BAD'} {mp.nstr(d, 25)}  err={mp.nstr(err, 3)}")
    # Euler products: checked against the elementary partial product.
    for label, skip in (("iprod-euler-zeta2", False), ("iprod-euler-from2", True)):
        p = mpf(1)
        for q in primes_below(2_000_000):
            if skip and q == 2:
                continue
            p /= (1 - mpf(q)**-2)
        ref = pi**2/8 if skip else pi**2/6
        err = abs(p - ref)
        ok = err < mpf(10)**-5      # partial product to 2e6: ~1e-7 tail
        bad += not ok
        print(f"{label:22s} {'ok ' if ok else 'BAD'} {mp.nstr(p, 12)}  err={mp.nstr(err, 3)}")
    print("BAD:", bad)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
