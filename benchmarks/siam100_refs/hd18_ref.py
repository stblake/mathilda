"""Reference A for holdout families hd18 (exp 223, Dirichlet eigenvalues) and
hd17 (exp 222, nodal-line lengths): closed forms evaluated in mpmath, plus
the eigenvalue ORDERING checks the ranks depend on.

hd18 (lambda_k of -Laplacian, zero Dirichlet data):
  1 equilateral triangle side 1, k=4: (16 Pi^2/9)(m^2+mn+n^2) at (2,2) = 64 Pi^2/3
  2 unit square, k=4: Pi^2 (m^2+n^2) at (2,2) = 8 Pi^2
  3 rectangle 2 x 1, k=1: Pi^2 (1/4 + 1) = 5 Pi^2/4
  4 right isosceles triangle, legs 1, k=1: Pi^2 (m^2+n^2), m>n>=1 -> 5 Pi^2
  5 30-60-90 triangle (half of the side-1 equilateral), k=1: 112 Pi^2/9
  6 unit disk, k=1: j_{0,1}^2
  7 half disk, k=1: j_{1,1}^2
  8 quarter disk, k=1: j_{2,1}^2
  9 sector of angle 2Pi/3, radius 1, k=1: j_{3/2,1}^2 (root of tan x = x)
 10 annulus 1/2 < r < 1, k=1: k^2, J0(k/2) Y0(k) = J0(k) Y0(k/2)
hd17 (length of the nodal set, in the open domain, of the k-th eigenfunction):
  1 unit square k=4: 2       2 rectangle 2x1 k=2: 1     3 equilateral k=4: 3/2
  4 unit disk k=6: 2 Pi j01/j02          5 half disk k=4: Pi j11/j12
  6 half disk k=2: 1         7 unit disk k=2: 2 (any diameter; lambda_2 double)
  8 unit square k=11: 4      9 unit disk k=15: 2 Pi (j01 + j02)/j03
 10 right isosceles triangle (legs 1) k=2: Sqrt2/2 (segment of x + y = 1)
Usage: python3.11 hd18_ref.py
"""
import mpmath as mp
mp.mp.dps = 80
pi = mp.pi
j = lambda nu, k: mp.besseljzero(nu, k)

# ---- ordering checks -------------------------------------------------------
def ranks(vals):
    vals = sorted(vals)
    return vals
sq = sorted(m * m + n * n for m in range(1, 10) for n in range(1, 10))
print("square m^2+n^2:", sq[:12])
tri = sorted(m * m + m * n + n * n for m in range(1, 10) for n in range(1, 10))
print("equilateral m^2+mn+n^2 (with multiplicity):", tri[:8])
rect = sorted((m * m / 4.0 + n * n, m, n) for m in range(1, 10) for n in range(1, 10))
print("rect 2x1:", rect[:4])
iso = sorted(m * m + n * n for m in range(1, 10) for n in range(1, m))
print("right isosceles:", iso[:5])
disk = sorted([(float(j(0, k)), 0, k) for k in range(1, 5)] +
              [(float(j(m, k)), m, k) for m in range(1, 9) for k in range(1, 5)] * 1)
seq = []
for v, m, k in disk:
    seq += [(v, m, k)] * (1 if m == 0 else 2)
print("disk ranks:", [(i + 1, round(v, 4), m, k) for i, (v, m, k) in enumerate(seq[:16])])
half = sorted((float(j(m, k)), m, k) for m in range(1, 9) for k in range(1, 4))
print("half disk:", half[:5])

# ---- values -----------------------------------------------------------------
x32 = mp.findroot(lambda x: mp.tan(x) - x, 4.49)
ann = mp.findroot(lambda k: mp.besselj(0, k / 2) * mp.bessely(0, k) - mp.besselj(0, k) * mp.bessely(0, k / 2), 6.2)
hd18 = [64 * pi ** 2 / 3, 8 * pi ** 2, 5 * pi ** 2 / 4, 5 * pi ** 2, 112 * pi ** 2 / 9,
        j(0, 1) ** 2, j(1, 1) ** 2, j(2, 1) ** 2, x32 ** 2, ann ** 2]
print("j_{3/2,1} check (besseljzero):", mp.nstr(j(mp.mpf(3) / 2, 1) - x32, 3))
for i, v in enumerate(hd18):
    print("hd18.%d %s" % (i + 1, mp.nstr(v, 72)))
hd17 = [2, 1, mp.mpf(3) / 2, 2 * pi * j(0, 1) / j(0, 2), pi * j(1, 1) / j(1, 2), 1, 2, 4,
        2 * pi * (j(0, 1) + j(0, 2)) / j(0, 3), mp.sqrt(2) / 2]
for i, v in enumerate(hd17):
    print("hd17.%d %s" % (i + 1, mp.nstr(v, 72)))
