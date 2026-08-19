#!/usr/bin/env python3
"""benchmarks/88-thue-equations/cases.py -- the single source of truth.

An adversarial stress corpus for the Thue-equation solver
(`Solve[F(x,y) == m && Element[{x,y}, Integers], {x,y}, Integers]`,
src/solvethue.c).  The point is to BREAK it: find wrong answers, crashes,
timeouts, and bottlenecks -- not to look good.

Every case is cross-checked against PARI/GP's `thue()`, the gold-standard
reference (it solves arbitrary Thue equations: any degree, any m, non-
monogenic fields).  A case where Mathilda's finite set differs from PARI's is
a COMPLETENESS BUG.

Each case is defined by the DEFINING POLYNOMIAL  P(x) = F(x, 1)  as a coeff
vector ascending in x (`pcoef[i]` = coeff of x^i, degree n = len-1) and the
right-hand side `m`.  The homogeneous form is  F(x,y) = sum_i pcoef[i] x^i y^(n-i).
Both the Mathilda `Solve[...]` and the PARI `thue(thueinit(P,0), m)` are built
from these, so the two halves solve the SAME equation.

`family` groups the report; `note` is a one-line remark; `expect` (optional) is
a literature/known count used only as an extra anchor -- PARI is the oracle.
"""

# ----------------------------------------------------------------------------
#  builders: coeff-vector -> Mathilda form string / PARI polynomial string
# ----------------------------------------------------------------------------

def mathilda_form(pcoef):
    """F(x,y) = sum_i pcoef[i] x^i y^(n-i)  as a Mathilda/WL string."""
    n = len(pcoef) - 1
    terms = []
    for i, c in enumerate(pcoef):
        if c == 0:
            continue
        xp = "" if i == 0 else ("x" if i == 1 else f"x^{i}")
        yp = "" if (n - i) == 0 else ("y" if (n - i) == 1 else f"y^{n-i}")
        mon = "*".join(t for t in (xp, yp) if t) or "1"
        terms.append((c, mon))
    s = ""
    for k, (c, mon) in enumerate(terms):
        sign = "-" if c < 0 else "+"
        mag = abs(c)
        piece = mon if mag == 1 and mon != "1" else (f"{mag}" if mon == "1" else f"{mag}*{mon}")
        if k == 0:
            s = ("-" + piece) if c < 0 else piece
        else:
            s += f" {sign} {piece}"
    return s

def pari_poly(pcoef):
    """P(x) = sum_i pcoef[i] x^i  as a PARI string."""
    n = len(pcoef) - 1
    terms = []
    for i, c in enumerate(pcoef):
        if c == 0:
            continue
        xp = "1" if i == 0 else ("x" if i == 1 else f"x^{i}")
        terms.append(f"({c})*{xp}" if xp != "1" else f"({c})")
    return " + ".join(terms) if terms else "0"


CASES = []
def add(label, family, pcoef, m, note="", expect=None, reps=1):
    CASES.append(dict(label=label, family=family, pcoef=list(pcoef), m=m,
                      note=note, expect=expect, reps=reps,
                      form=mathilda_form(pcoef), P=pari_poly(pcoef)))


# === A. Binomial cubics  x^3 - d y^3 = m  (Delone-Nagell) ====================
# monogenic iff d squarefree and d^2 !~ 1 (mod 9); m = +-1 in scope.
for d in [2, 3, 5, 6, 7, 10, 11, 12, 13, 15, 20, 30, 41, 42, 43]:
    add(f"binom-cubic-d{d}-p1", "binomial-cubic", [-d, 0, 0, 1], 1,
        note=f"x^3 - {d} y^3 = 1")
    add(f"binom-cubic-d{d}-m1", "binomial-cubic", [-d, 0, 0, 1], -1,
        note=f"x^3 - {d} y^3 = -1")
# large fundamental solution (stress the bound reduction reaching it)
add("binom-cubic-d17-p1", "binomial-cubic-large", [-17, 0, 0, 1], 1,
    note="x^3-17y^3=1 -> (18,7); Q(cbrt17) non-monogenic (index 3)")
add("binom-cubic-d20-p1", "binomial-cubic-large", [-20, 0, 0, 1], 1,
    note="x^3-20y^3=1")
add("binom-cubic-plus-d2", "binomial-cubic", [2, 0, 0, 1], 1, note="x^3+2y^3=1")
add("binom-cubic-plus-d7", "binomial-cubic", [7, 0, 0, 1], 1, note="x^3+7y^3=1")

# === B. |m| > 1  (out of the |m|=1 scope -> Mathilda should DECLINE) =========
for m in [2, 3, 4, 5, 9, 10, 100]:
    add(f"binom-cubic-d2-m{m}", "high-m", [-2, 0, 0, 1], m, note=f"x^3-2y^3={m}")
add("binom-cubic-d2-m73", "high-m", [-2, 0, 0, 1], 73, note="x^3-2y^3=73 (has solutions)")
add("cyclic-cubic-m2", "high-m", [1, -3, 0, 1], 2, note="x^3-3xy^2+y^3=2")

# === C. Simplest cubics / cyclic cubic fields (Thomas family) ================
add("cyclic-cond9-p1", "simplest-cubic", [1, -3, 0, 1], 1, note="cond 9, x^3-3xy^2+y^3=1")
add("cyclic-cond9-m1", "simplest-cubic", [1, -3, 0, 1], -1, note="cond 9, =-1")
add("thomas-cond7-p1", "simplest-cubic", [-1, -2, 1, 1], 1,
    note="x^3+x^2y-2xy^2-y^3=1 (cond 7, 9 solutions)")
add("thomas-cond7-m1", "simplest-cubic", [-1, -2, 1, 1], -1, note="cond 7, =-1")
# Thomas family x^3 - t x^2 y - (t+2) x y^2 - y^3 = 1, various t
for t in [1, 2, 3, 4, 5, 10]:
    add(f"thomas-t{t}", "simplest-cubic-family", [-1, -(t+2), -t, 1], 1,
        note=f"x^3-{t}x^2y-{t+2}xy^2-y^3=1")
add("shanks-simplest-n3", "simplest-cubic", [1, -3, -6, 1], 1,
    note="Shanks x^3-3x^2y-6xy^2-... variant")

# === D. General cubics (assorted discriminants, signatures) =================
add("cubic-irr-1", "general-cubic", [1, 1, -2, 1], 1, note="x^3-2x^2y+xy^2+y^3")
add("cubic-irr-2", "general-cubic", [-1, 3, -3, 1], 1, note="totally real-ish")
add("cubic-irr-3", "general-cubic", [1, 0, -1, 1], 1, note="x^3-xy^2+y^3? small disc")
add("cubic-neg-disc", "general-cubic", [1, 1, 1, 1], 1, note="x^3+x^2y+xy^2+y^3 (reducible? check)")
add("cubic-big-coef", "general-cubic", [-5, 11, -7, 1], 1, note="bigger coeffs")

# === E. No-solution cubics (must return {} == PARI) =========================
add("nosol-d2-m4", "no-solution", [-2, 0, 0, 1], 4, note="x^3-2y^3=4 (|m|!=1 -> decline; PARI {})")
add("nosol-cyclic", "no-solution", [1, -3, 0, 1], 4, note="cyclic, =4")

# === F. Binomial quartics  x^4 - d y^4 = m  (subfield -> degenerate lattice) =
for d in [2, 3, 5, 6, 7, 8, 10, 17]:
    add(f"binom-quartic-d{d}-p1", "binomial-quartic", [-d, 0, 0, 0, 1], 1,
        note=f"x^4 - {d} y^4 = 1")
    add(f"binom-quartic-d{d}-m1", "binomial-quartic", [-d, 0, 0, 0, 1], -1,
        note=f"x^4 - {d} y^4 = -1 (subfield unit -> case iii)")
add("ljunggren-x4-2y4", "binomial-quartic", [-2, 0, 0, 0, 1], -1,
    note="Ljunggren-type x^4-2y^4=-1")

# === G. General quartics (signatures (4,0),(2,1),(0,2); biquadratic subfield)=
add("quartic-biquad-1", "general-quartic", [1, 0, -4, 0, 1], 1,
    note="x^4-4x^2y^2+y^4 (biquadratic, quadratic subfield)")
add("quartic-tr-1", "general-quartic", [1, -4, 6, -4, 1], 1, note="(x-y)^4+... totally real?")
add("quartic-mixed-1", "general-quartic", [-1, 0, 0, 2, 1], 1, note="x^4+2x^3y-y^4")
add("quartic-cyclotomic", "general-quartic", [1, 1, 1, 1, 1], 1,
    note="x^4+x^3y+x^2y^2+xy^3+y^4 (cyclotomic, reducible over cyclo)")
add("quartic-x4-x-1", "general-quartic", [-1, -1, 0, 0, 1], 1, note="x^4-xy^3-y^4")
add("quartic-totallyimag", "general-quartic", [1, 0, 1, 0, 1], 1,
    note="x^4+x^2y^2+y^4 (may be reducible)")

# === H. Quintics / sextics (rank up to 4/5 -> unit search / Baker stress) ====
for d in [2, 3, 5]:
    add(f"binom-quintic-d{d}", "quintic", [-d, 0, 0, 0, 0, 1], 1, note=f"x^5-{d}y^5=1")
add("quintic-general", "quintic", [1, -1, 0, 0, -1, 1], 1, note="general quintic form")
for d in [2, 3]:
    add(f"binom-sextic-d{d}", "sextic", [-d, 0, 0, 0, 0, 0, 1], 1, note=f"x^6-{d}y^6=1")
add("septic-x7-2", "septic", [-2, 0, 0, 0, 0, 0, 0, 1], 1, note="x^7-2y^7=1")

# === I. Non-monogenic fields (Gate 1 decline; PARI still solves) =============
add("nonmono-d17", "non-monogenic", [-17, 0, 0, 1], 1, note="Q(cbrt17) index 3 at p=3")
add("nonmono-d19", "non-monogenic", [-19, 0, 0, 1], 1, note="x^3-19y^3=1")
add("nonmono-cubic-idx", "non-monogenic", [-8, -2, -1, 1], 1,
    note="Dedekind cubic t^3-t^2-2t-8 (index 2 at 2)")

# === J. Reducible forms (should DECLINE; PARI errors or handles factors) =====
add("reducible-diff-cubes", "reducible", [0, 0, 0, 1], None, note="placeholder-skip")  # replaced below
CASES.pop()
add("reducible-x3-y3", "reducible", [-1, 0, 0, 1], 7,
    note="x^3 - y^3 = 7 -> (x-y)(x^2+xy+y^2); reducible")
add("reducible-x4-y4", "reducible", [-1, 0, 0, 0, 1], 15, note="x^4 - y^4 = 15 reducible")
add("reducible-biquad", "reducible", [4, 0, -5, 0, 1], 3,
    note="x^4-5x^2y^2+4y^4=(x^2-y^2)(x^2-4y^2)")

# === K. Adversarial precision (roots close / big disc / big coeffs) =========
add("adv-close-roots-1", "adversarial-precision", [1, -7, 14, -7, 1], 1,
    note="near-symmetric quartic, clustered roots")
add("adv-big-coef-cubic", "adversarial-precision", [-97, 0, 0, 1], 1, note="x^3-97y^3=1")
add("adv-big-coef-2", "adversarial-precision", [131, -211, 97, 1], 1, note="large mixed coeffs")
add("adv-scaled-1", "adversarial-precision", [-999, 0, 0, 1], 1, note="x^3-999y^3=1 (999=3^3*37)")
add("adv-large-reg", "adversarial-precision", [-1, -1, -1, 1], 1,
    note="x^3-x^2y-xy^2-y^3=1 (tribonacci field, larger regulator)")

# === L. Many-solution / large-solution ======================================
add("many-sol-cond7", "many-solutions", [-1, -2, 1, 1], 1, note="Thomas cond7 (9 solutions)")
add("large-sol-1", "large-solutions", [-3, 0, 0, 1], 1, note="x^3-3y^3=1")


if __name__ == "__main__":
    print(f"{len(CASES)} cases")
    for c in CASES[:6]:
        print(f"  {c['label']:28s} F = {c['form']}  == {c['m']}   | PARI P = {c['P']}")
