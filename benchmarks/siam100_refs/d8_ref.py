"""Reference for d8: c_10(R^10), the Sobolev constant.  CONSENSUS digits only.

c^2 = (2 pi)^-N Int_{R^N} dx / w(x),  w = Sum_{|k|<=m} x^{2k} = h_m(1, x_1^2, ..., x_N^2)
(h_m = complete homogeneous symmetric polynomial), N = m = 10.

Reduction (exact): x = sqrt(rho) * omega, omega uniform on S^{N-1}, so
u = omega^2 ~ Dirichlet(1/2, ..., 1/2) and

  c^2 = (2 pi)^-N * Gamma(1/2)^N / Gamma(N/2) * E_u[ F(u) ],
  F(u) = Int_0^inf rho^(N/2-1) / P_u(rho) d rho,   P_u(rho) = Sum_{j<=m} h_j(u) rho^j,

with F done in closed form by residues over the roots rho_k of P_u:
  F(u) = - Sum_k rho_k^(N/2-1) log(-rho_k) / P_u'(rho_k)       (N even).
(Checked: N=1, m=3 gives c = 0.5 exactly; N=2, m=2 agrees with a direct 2-D
NIntegrate to 7 digits.)

E_u[F] is a 9-dimensional integral of a smooth symmetric function.  It is
estimated with a control variate: all products of the power sums
p_2..p_10 of u up to degree 4 (990 terms), whose Dirichlet expectations are
EXACT rationals (moments of the Dirichlet law), fitted by least squares on an
independent pilot sample; the residual has ~1/1000 of F's spread.

Two estimators (run both; they share F and the control variate, not the sampler):
  qmc : scrambled Sobol points -> Gaussian -> omega = z/|z|, R independent
        scramblings for the error bar;
  mc  : numpy PCG64 pseudo-random stick-breaking, v_i ~ Beta(1/2, (N-i)/2),
        u_i = v_i prod_{j<i}(1 - v_j), R independent streams.
Usage: python d8_ref.py {qmc|mc} log2(points per rep) reps
"""
import sys, math, itertools
from fractions import Fraction
import numpy as np
from scipy.stats import qmc, norm

N = m = 10
DEG = 4

def F_of(u):
    n = u.shape[0]
    p = np.stack([(u ** j).sum(1) for j in range(1, m + 1)], 1)
    h = np.zeros((n, m + 1)); h[:, 0] = 1
    for j in range(1, m + 1):
        h[:, j] = sum(p[:, i - 1] * h[:, j - i] for i in range(1, j + 1)) / j
    Cm = np.zeros((n, m, m)); Cm[:, 1:, :-1] = np.eye(m - 1)
    Cm[:, :, -1] = -h[:, :m] / h[:, m:m + 1]
    rts = np.linalg.eigvals(Cm)
    dP = sum(j * h[:, j:j + 1] * rts ** (j - 1) for j in range(1, m + 1))
    return (-(rts ** (N // 2 - 1) * np.log(-rts.astype(complex)) / dP).sum(1)).real

# ---- control variates: products of power sums p_2..p_10, exact means ----
combos = [c for d in range(1, DEG + 1) for c in itertools.combinations_with_replacement(range(2, 11), d)]

def poch(a, s):
    r = Fraction(1)
    for i in range(s):
        r *= a + i
    return r

def dir_moment(blocks):
    """E[u_1^{s_1} ... u_r^{s_r}] (distinct indices), Dirichlet(1/2 x N)."""
    S = sum(blocks)
    num = Fraction(1)
    for s in blocks:
        num *= poch(Fraction(1, 2), s)
    return num / poch(Fraction(N, 2), S)

def set_partitions(lst):
    if not lst:
        yield []
        return
    first, rest = lst[0], lst[1:]
    for part in set_partitions(rest):
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]
        yield [[first]] + part

def mean_of(c):
    """E[p_{c1} p_{c2} ...] = Sum over set partitions of the factors."""
    tot = Fraction(0)
    for part in set_partitions(list(c)):
        r = len(part)
        ways = math.perm(N, r)
        tot += ways * dir_moment([sum(b) for b in part])
    return tot

EX = np.array([float(mean_of(c)) for c in combos])

def design(u):
    p = {j: (u ** j).sum(1) for j in range(2, 11)}
    return np.stack([np.prod([p[j] for j in c], 0) for c in combos], 1)

def sample(kind, rng_or_sob, n):
    if kind == "qmc":
        z = norm.ppf(rng_or_sob)
        u = z * z
        return u / u.sum(1, keepdims=True)
    rng = rng_or_sob
    u = np.empty((n, N)); rem = np.ones(n)
    for i in range(1, N):
        v = rng.beta(0.5, (N - i) / 2.0, n)
        u[:, i - 1] = v * rem; rem *= 1 - v
    u[:, N - 1] = rem
    return u

if __name__ == "__main__":
    kind, L, R = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    # pilot fit (independent of the estimation samples)
    pr = np.random.default_rng(12345)
    up = sample("mc", pr, 300000)
    fp = F_of(up); Xp = design(up)
    beta, *_ = np.linalg.lstsq(Xp - EX, fp - fp.mean(), rcond=None)
    print("pilot: F std %.3g, residual std %.3g" % (fp.std(), (fp - fp.mean() - (Xp - EX) @ beta).std()), flush=True)
    const = (2 * math.pi) ** -N * math.gamma(0.5) ** N / math.gamma(N / 2)
    ests = []
    CH = 2 ** 15
    for r in range(R):
        tot, cnt = 0.0, 0
        if kind == "qmc":
            pts = qmc.Sobol(d=N, scramble=True, seed=2000 + r).random_base2(L)
            chunks = (pts[b:b + CH] for b in range(0, pts.shape[0], CH))
            for ch in chunks:
                u = sample("qmc", ch, ch.shape[0])
                tot += (F_of(u) - (design(u) - EX) @ beta).sum(); cnt += u.shape[0]
        else:
            rng = np.random.default_rng(7000 + r)
            for b in range(0, 2 ** L, CH):
                u = sample("mc", rng, CH)
                tot += (F_of(u) - (design(u) - EX) @ beta).sum(); cnt += CH
        c = math.sqrt(const * tot / cnt)
        ests.append(c)
        print(kind, "rep", r, repr(c), flush=True)
    e = np.array(ests)
    print("%s mean %.12g  std-err %.3g" % (kind, e.mean(), e.std(ddof=1) / math.sqrt(len(e))), flush=True)
