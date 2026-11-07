# Holdout-A D7 family (experiment 212) references, method A: mpmath tanh-sinh
# of the explicit 1-D form (x-interval length w(y), or the y-extent as a
# function of x) with breakpoints at every island end and every cube-root kink,
# island ends located by a scan + findroot, dps 90.  Prints island ends / kinks
# (starting values for the cases).  Method B: hd7_ref.wl.
import mpmath as mp
mp.mp.dps = 90
pi = mp.pi
rc = lambda t, n: mp.sign(t) * abs(t) ** (mp.mpf(1) / n)

def islands(g, lo=-12, hi=12, n=24000):
    """Intervals of y where g(y) < 1 (g(0) = 1 is an exact end)."""
    ys = [lo + (hi - lo) * mp.mpf(i) / n for i in range(n + 1)]
    h = lambda y: g(y) - 1
    ends = []
    for a, b in zip(ys, ys[1:]):
        if b == 0 or a == 0:
            continue
        if h(a) * h(b) < 0:
            ends.append(mp.findroot(h, (a, b), solver='anderson'))
    ends.append(mp.mpf(0)); ends = sorted(ends)
    out = []
    for a, b in zip(ends, ends[1:]):
        if g((a + b) / 2) < 1:
            out.append((a, b))
    return out

def kinks(g, sc, a, b, n=2000):
    """Zeros of -y - sqrt(-sc*log g(y)) strictly inside (a, b)."""
    k = lambda y: -y - mp.sqrt(-sc * mp.log(g(y)))
    out = []
    ys = [a + (b - a) * mp.mpf(i) / n for i in range(1, n)]
    for u, v in zip(ys, ys[1:]):
        if k(u) * k(v) < 0:
            out.append(mp.findroot(k, (u, v), solver='anderson'))
    return out

def area_fg(g, inv, sc=1):
    """Area of {exp(-(y + p(x))^2/sc) > g(y)}, p^-1 = inv."""
    tot = 0
    for a, b in islands(g):
        ks = kinks(g, sc, a, b)
        print("   island", mp.nstr(a, 20), mp.nstr(b, 20), "kinks", [mp.nstr(k, 20) for k in ks])
        def w(y):
            s = mp.sqrt(-sc * mp.log(g(y)))
            return inv(-y + s) - inv(-y - s)
        tot += mp.quad(w, [a] + ks + [b])
    return tot

g32 = lambda y: y**2 / 32 + mp.exp(mp.sin(y))
g16 = lambda y: y**2 / 16 + mp.exp(mp.sin(y))
g100 = lambda y: y**2 / 100 + mp.exp(mp.sin(y))
R = {}
R[1] = mp.quad(lambda x: 2 * (1 - x**4) ** mp.mpf('0.25'), [-1, 0, 1])
a2 = mp.findroot(lambda x: mp.cos(x) - x**2, 0.8)
R[2] = mp.quad(lambda x: mp.cos(x) - x**2, [-a2, 0, a2])
print("hd7.2 end", mp.nstr(a2, 30))
print("hd7.3"); R[3] = area_fg(g16, lambda t: rc(t, 3))
print("hd7.4"); R[4] = area_fg(g32, mp.asinh)
print("hd7.5"); R[5] = area_fg(g32, lambda t: rc(t, 5))
print("hd7.6"); R[6] = area_fg(g32, lambda t: rc(t, 3), sc=2)
R[7] = mp.quad(lambda x: pi - 2 * mp.asin(1 - mp.sin(x)), [0, pi / 2, pi])
L = mp.log(2 + mp.sqrt(3))
R[8] = mp.quad(lambda x: mp.log(4 - mp.exp(x)) + x, [-L, 0, L])
a9 = mp.findroot(lambda x: 3 - x**2 - mp.exp(x), -1.7); b9 = mp.findroot(lambda x: 3 - x**2 - mp.exp(x), 0.8)
print("hd7.9 ends", mp.nstr(a9, 30), mp.nstr(b9, 30))
R[9] = mp.quad(lambda x: 3 - x**2 - mp.exp(x), [a9, 0, b9])
print("hd7.10"); R[10] = area_fg(g100, lambda t: rc(t, 3))
closed = {1: 4 * mp.gamma(mp.mpf(5) / 4) ** 2 / mp.gamma(mp.mpf(3) / 2),
          2: 2 * mp.sin(a2) - 2 * a2**3 / 3,
          9: 3 * (b9 - a9) - (b9**3 - a9**3) / 3 - (mp.exp(b9) - mp.exp(a9))}
for k in range(1, 11):
    s = "hd7.%d  %s" % (k, mp.nstr(R[k], 75))
    if k in closed:
        s += "  closed-diff %s" % mp.nstr(R[k] - closed[k], 3)
    print(s, flush=True)
