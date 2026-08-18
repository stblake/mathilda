#!/usr/bin/env python3
"""benchmarks/87-diophantine-integers/cases.py -- the single source of truth.

One row per well-known Diophantine family.  Each case carries:

  * label     -- join key between the Mathilda and sympy halves
  * family    -- the mathematical family (for the report grouping)
  * category  -- A: sympy `diophantine` returns a finite set we can count
                    directly and compare;
                 B: sympy `diophantine` solves the *family* but returns an
                    unbounded PARAMETRIC set -- the bounded answer the user
                    actually asked for needs a hand-written search on top;
                 C: sympy `diophantine` raises NotImplementedError -- there is
                    no solver, so a naive Python search is the only fallback.
  * math      -- the exact Solve[..., Integers] Mathilda evaluates
  * reps      -- timed repetitions on the Mathilda side (min is reported)
  * expect    -- the known solution count (ground-truth anchor from the
                 literature / cross-checks); None where only agreement matters
  * sympy     -- () -> (status, count|None, seconds) using sympy.diophantine
  * brute     -- () -> count, a reasonable Python search over the SAME box
                 (what a sympy user writes when diophantine can't answer); may
                 exceed the wall-clock budget, in which case run.py records
                 TIMEOUT
  * box       -- human description of the search box (for the report)
  * note      -- one-line remark

The Mathilda expressions and the Python searches are deliberately over the
*same* box so the timings are head-to-head.
"""

from math import gcd, isqrt, prod
from time import perf_counter

from sympy import symbols
from sympy.solvers.diophantine import diophantine

# ----------------------------------------------------------------------------
#  sympy driver: classify what diophantine returns for one polynomial.
# ----------------------------------------------------------------------------

def run_sympy(poly, syms, permute=False):
    """Return (status, count|None, seconds).

    status is one of:
      'finite'      -- a finite set of concrete integer tuples (count is exact)
      'parametric'  -- a family whose tuples still contain free parameters
      'empty'       -- provably no solution (e.g. gcd rules it out)
      'unsupported' -- NotImplementedError: no solver written
    """
    t = perf_counter()
    try:
        s = diophantine(poly, permute=permute)
    except NotImplementedError:
        return ('unsupported', None, perf_counter() - t)
    except Exception:                                  # keep the sweep going
        return ('unsupported', None, perf_counter() - t)
    secs = perf_counter() - t
    if len(s) == 0:
        return ('empty', 0, secs)
    params = set()
    for tup in s:
        for comp in tup:
            params |= comp.free_symbols
    params -= set(syms)
    if params:
        return ('parametric', None, secs)
    return ('finite', len(s), secs)


# ----------------------------------------------------------------------------
#  Python "search" fallbacks -- the same box the Mathilda case uses.
#  Each returns the number of solutions.  They are written the way a sympy
#  user would when `diophantine` cannot answer: solve the innermost variable
#  in closed form wherever possible instead of looping over it.
# ----------------------------------------------------------------------------

def _is_sq(n):
    return n >= 0 and isqrt(n) ** 2 == n

def brute_two_squares(N):
    B = isqrt(N)
    c = 0
    for x in range(-B, B + 1):
        r = N - x * x
        if _is_sq(r):
            c += 2 if r > 0 else 1          # y = +/-sqrt(r); r==0 -> one y
    return c

def brute_three_squares(N):
    B = isqrt(N)
    c = 0
    for x in range(-B, B + 1):
        for y in range(-B, B + 1):
            r = N - x * x - y * y
            if _is_sq(r):
                c += 2 if r > 0 else 1
    return c

def brute_linear_pos(a, rhs):
    # sum a_i x_i == rhs, all x_i > 0 ; solve the last coefficient in closed form
    a0, a1, a2 = a
    c = 0
    x = 1
    while a0 * x < rhs:
        rem1 = rhs - a0 * x
        y = 1
        while a1 * y < rem1:
            rem2 = rem1 - a1 * y
            if rem2 > 0 and rem2 % a2 == 0:
                c += 1                       # z = rem2/a2 > 0
            y += 1
        x += 1
    return c

def brute_linear_unsolvable(a, rhs, _bound):
    g = 0
    for ai in a:
        g = gcd(g, ai)
    return 0 if rhs % g != 0 else None       # None => would need real search

def brute_pell(D, N, xbound):
    # x^2 - D y^2 == N, 0 < x < xbound, y > 0
    c = 0
    y = 1
    while True:
        val = D * y * y + N
        if val <= 0:
            y += 1
            continue
        x = isqrt(val)
        if x * x == val:
            if x <= 0 or x >= xbound:
                break
            if x > 0:
                c += 1
        # advance; stop once the smallest possible x exceeds the bound
        if isqrt(D * y * y + N) >= xbound:
            break
        y += 1
    return c

def brute_pythagorean_perimeter(P):
    c = 0
    for x in range(1, P):
        for y in range(x + 1, P):
            z = P - x - y
            if z <= y:
                continue
            if x * x + y * y == z * z:
                c += 1
    return c

def brute_motivating():
    c = 0
    y = 1
    while 2 * y ** 3 < 3681:
        r = 3681 - 2 * y ** 3
        if _is_sq(r) and r > 0:
            c += 1                           # x = +sqrt(r) > 0
        y += 1
    return c

def brute_two_quadratics():
    c = 0
    B = isqrt(1000)
    for x in range(1, B + 1):
        for y in range(1, B + 1):
            zz = 1000 - x * x - y * y
            if zz <= 0 or not _is_sq(zz):
                continue
            z = isqrt(zz)
            if z > 0 and x * x + 2 * y * y - z * z == 200:
                c += 1
    return c

def brute_egyptian():
    # 4/2027 == 1/x + 1/y + 1/z, 0 < x <= y <= z
    from math import ceil
    R = (4, 2027)                            # rational 4/2027
    c = 0
    # x in [ceil(1/R)... floor(3/R)] with R = 4/2027 -> [ceil(2027/4)..floor(3*2027/4)]
    xlo = ceil(2027 / 4)
    xhi = (3 * 2027) // 4
    for x in range(max(1, xlo), xhi + 1):
        # remaining = 4/2027 - 1/x = (4x - 2027) / (2027 x)
        rn = 4 * x - 2027
        rd = 2027 * x
        if rn <= 0:
            continue
        # 1/y + 1/z == rn/rd, y <= z, y >= x, y in [ceil(rd/rn) .. floor(2 rd/rn)]
        ylo = max(x, -(-rd // rn))            # ceil(rd/rn)
        yhi = (2 * rd) // rn
        for y in range(ylo, yhi + 1):
            # 1/z == rn/rd - 1/y == (rn*y - rd) / (rd*y)
            zn = rn * y - rd
            zd = rd * y
            if zn <= 0:
                continue
            if zd % zn == 0:
                z = zd // zn
                if z >= y:
                    c += 1
    return c

def brute_markov(bound):
    # x^2+y^2+z^2 == 3xyz, 0 < x <= y <= z <= bound ; naive triple loop
    c = 0
    for x in range(1, bound + 1):
        for y in range(x, bound + 1):
            s = x * x + y * y
            p = 3 * x * y
            # z^2 - p z + s == 0 -> integer z
            disc = p * p - 4 * s
            if disc < 0 or not _is_sq(disc):
                continue
            r = isqrt(disc)
            for z in ((p + r) // 2, (p - r) // 2):
                if (2 * z == p + r or 2 * z == p - r) and y <= z <= bound and z * z - p * z + s == 0:
                    c += 1
    return c

def brute_sum3cubes_ordered(target, hi):
    # x^3+y^3+z^3 == target, 0 < x <= y <= z < hi
    c = 0
    for x in range(1, hi):
        x3 = x ** 3
        if 3 * x3 > target:
            break
        for y in range(x, hi):
            s = x3 + y ** 3
            if s + y ** 3 > target:          # z >= y so s + y^3 is the min
                break
            zz = target - s
            if zz < y ** 3:
                continue
            z = round(zz ** (1 / 3))
            for zc in (z - 1, z, z + 1):
                if y <= zc < hi and zc ** 3 == zz:
                    c += 1
                    break
    return c

def brute_sum3cubes_box(target, B):
    # x^3+y^3+z^3 == target, |x|,|y|,|z| < B  -- deliberately the full signed box
    c = 0
    for x in range(-B + 1, B):
        for y in range(-B + 1, B):
            zz = target - x ** 3 - y ** 3
            z = round(abs(zz) ** (1 / 3)) * (1 if zz >= 0 else -1)
            for zc in (z - 1, z, z + 1):
                if -B < zc < B and zc ** 3 == zz:
                    c += 1
                    break
    return c

def brute_mordell(k, xbound):
    # y^2 == x^3 + k, 0 < x < xbound, y > 0
    c = 0
    for x in range(1, xbound):
        r = x ** 3 + k
        if r > 0 and _is_sq(r):
            c += 1
    return c

def brute_taxicab(limit):
    # x^3+y^3 == z^3+w^3, 0<x<y, 0<z<w, x<z, x^3+y^3 < limit  (naive 4-index)
    c = 0
    xs = []
    x = 1
    while x ** 3 < limit:
        xs.append(x)
        x += 1
    for x in xs:
        for y in xs:
            if y <= x:
                continue
            s1 = x ** 3 + y ** 3
            if s1 >= limit:
                break
            for z in xs:
                if z <= x:
                    continue
                for w in xs:
                    if w <= z:
                        continue
                    if x ** 3 + y ** 3 == z ** 3 + w ** 3 and z ** 3 + w ** 3 < limit:
                        c += 1
    return c

def brute_euler_fifth():
    # a^5+b^5+c^5+d^5 == e^5, 0<a<=b<=c<=d<150, e<200 ; solve e by fifth-root
    c = 0
    for a in range(1, 150):
        a5 = a ** 5
        for b in range(a, 150):
            b5 = b ** 5
            for cc in range(b, 150):
                c5 = cc ** 5
                for d in range(cc, 150):
                    s = a5 + b5 + c5 + d ** 5
                    e = round(s ** 0.2)
                    for ec in (e - 1, e, e + 1):
                        if 0 < ec < 200 and ec ** 5 == s:
                            c += 1
                            break
    return c

def brute_catalan():
    # x^p - y^q == 1, 1<x<100, 1<y<100, 1<p<10, 1<q<10
    c = 0
    for x in range(2, 100):
        for p in range(2, 10):
            xp = x ** p
            for y in range(2, 100):
                for q in range(2, 10):
                    if xp - y ** q == 1:
                        c += 1
    return c

def brute_boxed_power():
    # x^2 + y^3 == z^5, 0<x<1000, 0<y<1000, 0<z<100 ; solve x by square-root
    c = 0
    for y in range(1, 1000):
        y3 = y ** 3
        for z in range(1, 100):
            r = z ** 5 - y3
            if r > 0 and _is_sq(r):
                x = isqrt(r)
                if 0 < x < 1000:
                    c += 1
    return c


# ----------------------------------------------------------------------------
#  The case table.
# ----------------------------------------------------------------------------

_x, _y, _z, _a, _b, _c, _d, _e = symbols('x y z a b c d e', integer=True)

CASES = [
    # ---- Category A: sympy diophantine returns a finite, directly comparable set
    dict(label="two-squares-25", family="Sum of two squares", category="A",
         math="Solve[x^2 + y^2 == 25, {x, y}, Integers]",
         reps=5, expect=12, box="|x|,|y| <= 5",
         sympy=lambda: run_sympy(_x**2 + _y**2 - 25, (_x, _y), permute=True),
         brute=lambda: brute_two_squares(25),
         note="x^2+y^2=N: fixed by even-only symmetric bound (was {} before)."),

    dict(label="two-squares-1e6", family="Sum of two squares", category="A",
         math="Solve[x^2 + y^2 == 1000000, {x, y}, Integers]",
         reps=5, expect=28, box="|x|,|y| <= 1000",
         sympy=lambda: run_sympy(_x**2 + _y**2 - 1000000, (_x, _y), permute=True),
         brute=lambda: brute_two_squares(1000000),
         note="Gaussian-integer count; Mathilda enumerates the exact leaf."),

    dict(label="three-squares-29", family="Sum of three squares", category="A",
         math="Solve[x^2 + y^2 + z^2 == 29, {x, y, z}, Integers]",
         reps=5, expect=72, box="|x|,|y|,|z| <= 5",
         sympy=lambda: run_sympy(_x**2 + _y**2 + _z**2 - 29, (_x, _y, _z), permute=True),
         brute=lambda: brute_three_squares(29),
         note="sympy needs permute=True to expand base reps to all 72 tuples."),

    # ---- Category B: sympy solves the family but only PARAMETRICALLY
    dict(label="linear-3var-box", family="Linear (bounded box)", category="B",
         math="Solve[15 x + 21 y + 35 z == 4000 && x > 0 && y > 0 && z > 0, {x, y, z}, Integers]",
         reps=3, expect=703, box="x,y,z > 0",
         sympy=lambda: run_sympy(15*_x + 21*_y + 35*_z - 4000, (_x, _y, _z)),
         brute=lambda: brute_linear_pos((15, 21, 35), 4000),
         note="diophantine returns a 2-parameter family; the positive count needs a search on top."),

    dict(label="linear-huge-coeff", family="Linear (unsolvable)", category="B",
         math="Solve[314159265 x + 271828182 y + 161803398 z == 1 "
              "&& Abs[x] < 10^6 && Abs[y] < 10^6 && Abs[z] < 10^6, {x, y, z}, Integers]",
         reps=3, expect=0, box="|x|,|y|,|z| < 10^6",
         sympy=lambda: run_sympy(314159265*_x + 271828182*_y + 161803398*_z - 1, (_x, _y, _z)),
         brute=lambda: 0,     # gcd(coeffs)=3 does not divide 1 -> provably empty
         note="gcd(coeffs)=3 does not divide 1: no solution. sympy sees it via gcd instantly."),

    dict(label="pell-61", family="Pell x^2 - D y^2 = 1", category="B",
         math="Solve[x^2 - 61 y^2 == 1 && x > 0 && y > 0 && x < 10^10, {x, y}, Integers]",
         reps=5, expect=1, box="0 < x < 10^10",
         sympy=lambda: run_sympy(_x**2 - 61*_y**2 - 1, (_x, _y)),
         brute=lambda: brute_pell(61, 1, 10**10),
         note="D=61 has a famously large fundamental solution (1766319049, 226153980)."),

    dict(label="pell-negative", family="Negative Pell x^2 - D y^2 = -1", category="B",
         math="Solve[x^2 - 2 y^2 == -1 && 0 < x < 1000 && y > 0, {x, y}, Integers]",
         reps=5, expect=4, box="0 < x < 1000",
         sympy=lambda: run_sympy(_x**2 - 2*_y**2 + 1, (_x, _y)),
         brute=lambda: brute_pell(2, -1, 1000),
         note="Negative Pell is solvable iff the sqrt(D) CF period is odd; D=2 qualifies."),

    dict(label="pythagorean-perimeter", family="Pythagorean + linear system", category="B",
         math="Solve[x^2 + y^2 == z^2 && x + y + z == 3000 && 0 < x < y && z > 0, {x, y, z}, Integers]",
         reps=3, expect=3, box="perimeter 3000, 0 < x < y",
         sympy=lambda: run_sympy(_x**2 + _y**2 - (3000 - _x - _y)**2, (_x, _y)),
         sympy_caveat="raw conic after manual z-elimination; only 3 of the 252 satisfy 0<x<y<z",
         brute=lambda: brute_pythagorean_perimeter(3000),
         note="A 2-equation system: diophantine takes one equation, so z is eliminated by hand, "
              "then the 252 conic points are filtered down to the 3 triangles."),

    # ---- Category C: sympy diophantine raises NotImplementedError
    dict(label="motivating", family="Mixed quadratic-cubic", category="C",
         math="Solve[x^2 + 2 y^3 == 3681 && x > 0 && y > 0, {x, y}, Integers]",
         reps=5, expect=3, box="x,y > 0",
         sympy=lambda: run_sympy(_x**2 + 2*_y**3 - 3681, (_x, _y)),
         brute=lambda: brute_motivating(),
         note="The motivating case: {15,12},{41,10},{57,6}."),

    dict(label="two-quadratics", family="Two-quadratic system", category="C",
         math="Solve[x^2 + y^2 + z^2 == 1000 && x^2 + 2 y^2 - z^2 == 200 "
              "&& x > 0 && y > 0 && z > 0, {x, y, z}, Integers]",
         reps=5, expect=0, box="x,y,z > 0",
         sympy=lambda: ('unsupported', None, 0.0),   # diophantine cannot take a system
         sympy_caveat="system: diophantine takes a single equation only (SympifyError on a list)",
         brute=lambda: brute_two_quadratics(),
         note="Positivity of the first equation bounds the whole box; no solution. "
              "sympy cannot express the coupled system at all."),

    dict(label="egyptian-4/2027", family="Egyptian fractions", category="C",
         math="Solve[4/2027 == 1/x + 1/y + 1/z && 0 < x <= y <= z, {x, y, z}, Integers]",
         reps=1, expect=73, box="0 < x <= y <= z",
         sympy=lambda: run_sympy(2027*(_y*_z + _x*_z + _x*_y) - 4*_x*_y*_z, (_x, _y, _z)),
         brute=lambda: brute_egyptian(),
         note="Unit-fraction expansion; 73 ordered decompositions."),

    dict(label="markov-1000", family="Markov triples", category="C",
         math="Solve[x^2 + y^2 + z^2 == 3 x y z && 0 < x <= y <= z <= 1000, {x, y, z}, Integers]",
         reps=1, expect=13, box="0 < x <= y <= z <= 1000",
         sympy=lambda: run_sympy(_x**2 + _y**2 + _z**2 - 3*_x*_y*_z, (_x, _y, _z)),
         brute=lambda: brute_markov(1000),
         note="A naive triple loop is 10^9; both Mathilda and the Python column instead solve the "
              "quadratic in z by its discriminant (O(bound^2)). One of the few cases the tuned "
              "Python search edges Mathilda -- sympy still has no Markov solver at all."),

    dict(label="sum3cubes-3", family="Sum of three cubes (small box)", category="C",
         math="Solve[x^3 + y^3 + z^3 == 3 && 0 < x <= y <= z && z < 100, {x, y, z}, Integers]",
         reps=5, expect=1, box="0 < x <= y <= z < 100",
         sympy=lambda: run_sympy(_x**3 + _y**3 + _z**3 - 3, (_x, _y, _z)),
         brute=lambda: brute_sum3cubes_ordered(3, 100),
         note="Only {1,1,1} in the positive box (the other, {4,4,-5}, is negative)."),

    dict(label="sum3cubes-42-bigbox", family="Sum of three cubes (10^5 box)", category="C",
         math="Solve[x^3 + y^3 + z^3 == 42 && Abs[x] < 10^5 && Abs[y] < 10^5 && Abs[z] < 10^5, {x, y, z}, Integers]",
         reps=1, expect=0, box="|x|,|y|,|z| < 10^5",
         sympy=lambda: run_sympy(_x**3 + _y**3 + _z**3 - 42, (_x, _y, _z)),
         brute=lambda: brute_sum3cubes_box(42, 10**5),
         note="8x10^15-point box: infeasible to brute; Mathilda uses the s=x+y | m divisor method."),

    dict(label="mordell-10000", family="Mordell y^2 = x^3 + k", category="C",
         math="Solve[y^2 == x^3 - 10000 && x > 0 && y > 0 && x < 10^4, {x, y}, Integers]",
         reps=3, expect=1, box="0 < x < 10^4",
         sympy=lambda: run_sympy(_y**2 - (_x**3 - 10000), (_x, _y)),
         brute=lambda: brute_mordell(-10000, 10**4),
         note="Elliptic curve; sympy reports cubic_thue NotImplemented."),

    dict(label="taxicab-1729", family="Taxicab (two cubes two ways)", category="C",
         math="Solve[x^3 + y^3 == z^3 + w^3 && 0 < x < y && 0 < z < w && x < z "
              "&& x^3 + y^3 < 100000, {x, y, z, w}, Integers]",
         reps=3, expect=10, box="x^3+y^3 < 10^5",
         sympy=lambda: run_sympy(_x**3 + _y**3 - _z**3 - _d**3, (_x, _y, _z, _d)),
         brute=lambda: brute_taxicab(100000),
         note="Naive 4-index search; Mathilda uses meet-in-the-middle."),

    dict(label="euler-5th-powers", family="Euler sum of like powers", category="C",
         math="Solve[a^5 + b^5 + c^5 + d^5 == e^5 && 0 < a <= b <= c <= d < 150 && e < 200, {a, b, c, d, e}, Integers]",
         reps=1, expect=1, box="ordered, d < 150, e < 200",
         sympy=lambda: run_sympy(_a**5 + _b**5 + _c**5 + _d**5 - _e**5, (_a, _b, _c, _d, _e)),
         brute=lambda: brute_euler_fifth(),
         note="The Lander-Parkin counterexample 27^5+84^5+110^5+133^5 = 144^5."),

    dict(label="catalan", family="Catalan / exponential", category="C",
         math="Solve[x^p - y^q == 1 && 1 < x < 100 && 1 < y < 100 && 1 < p < 10 && 1 < q < 10, {x, y, p, q}, Integers]",
         reps=5, expect=1, box="bases,exponents in (1,100)/(1,10)",
         sympy=lambda: ('unsupported', None, 0.0),   # variable exponent: not a polynomial
         brute=lambda: brute_catalan(),
         note="Variable exponents: not a polynomial at all. Only 3^2 - 2^3 = 1 (Mihailescu)."),

    dict(label="boxed-power", family="Bounded x^2 + y^3 = z^5", category="C",
         math="Solve[x^2 + y^3 == z^5 && 0 < x < 1000 && 0 < y < 1000 && 0 < z < 100, {x, y, z}, Integers]",
         reps=5, expect=2, box="0<x<1000, 0<y<1000, 0<z<100",
         sympy=lambda: run_sympy(_x**2 + _y**3 - _z**5, (_x, _y, _z)),
         brute=lambda: brute_boxed_power(),
         note="Explicit box + perfect-square leaf."),
]
