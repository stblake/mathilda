#!/usr/bin/env python3
"""
osp_famous.py -- the 222 "famous" problems of the math.osp.no benchmark
(suite 2026.1, CC BY 4.0), solved by Mathilda and graded locally.

    python3 benchmarks/osp_famous/osp_famous.py                 # solve + grade all
    python3 benchmarks/osp_famous/osp_famous.py --only basel_problem,ahmed_integral -v
    python3 benchmarks/osp_famous/osp_famous.py --runner run.jsonl   # grade a math-run file
    python3 benchmarks/osp_famous/osp_famous.py --regrade            # re-grade saved answers

TWO INPUTS.  The site drives Mathilda through its runner, math-run 0.1.0, whose
Mathilda adapter shares another adapter's translation.  That translation has no
form for `harmonic`, `integrate`, `sum` or `product` inside an expression and
passes them through as unknown lower-case heads (`harmonic[k]`,
`integrate[f, y, 0, 1]`), so nine problems are unanswerable as the runner poses
them.  This script therefore translates every problem itself, into the Mathilda
input the runner means (same heads, same `Assumptions -> ...`, same `mr`
prefix on the single capitals C D E I N O K), and grades that: the MATH score.
`--runner FILE` grades the answers in a math-run upload file instead: the
RUNNER score, which is what the site would report.

GRADING (local, by Mathilda itself).  The
reference answer is translated the same way and the check is op-specific:

  value ops (integrate, sum, product, limit, evaluate, simplify, laplace,
  determinant, ...)  the answer must be free of the task's own head, and
  N[answer - reference, 60] must vanish to 1e-30 (relative) at sample values of
  the free parameters chosen inside the assumptions.  simplify/evaluate also ask
  what the site asks: no head of the input that the reference has worked out
  (Cos, PolyLog, EllipticK, ArcTan, Root, ...) may survive, and the answer may
  not be much larger than the reference (LeafCount <= 1.5 ref + 4).
  solve_ode / solve_pde  the answer is substituted into the equation and the
  conditions and the residuals evaluated at sample points with the constants
  C[i] sampled; an ODE without conditions must carry `order` constants.  An
  implicit first-order answer F(x, y) == 0 is checked by implicit
  differentiation.
  solve_recurrence  the answer must satisfy the recurrence and the conditions
  at n = n0 .. n0 + 10, and agree with the reference there.
  solve  every answer must satisfy the equations; finite solution sets must
  equal the reference's numerically; a parametric family must contain the
  reference's k = 0 member.
  minimal_polynomial  the answer over the reference is a nonzero constant.
  holds  the answer is True.
A problem without a reference answer (open, or proven impossible: Painleve,
Blasius, MinimalPolynomial[Pi], ...) is graded "expected-unevaluated": solved
when Mathilda leaves it unevaluated, wrong when it claims an answer.

Grades: solved | wrong | unsolved (left unevaluated) | unclear (the checker
could not decide; manual review) | timeout | crash | error (a message or a
malformed answer).
"""
import argparse, json, os, re, subprocess, sys, tempfile, time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

# ---------------------------------------------------------------- translator
# samara's neutral syntax -> Mathilda input, as math-run's adapter writes it.
CAP = """exp log sqrt sin cos tan cot sec csc sinh cosh tanh coth sech csch abs sign
floor ceiling conjugate factorial gamma erf erfc erfi beta binomial pochhammer""".split()
FUNCS = {
    "re": "Re", "im": "Im", "arg": "Arg", "log_gamma": "LogGamma", "polygamma": "PolyGamma",
    "fresnel_s": "FresnelS", "fresnel_c": "FresnelC", "ei": "ExpIntegralEi",
    "li": "LogIntegral", "si": "SinIntegral", "ci": "CosIntegral", "shi": "SinhIntegral",
    "chi": "CoshIntegral", "expint": "ExpIntegralE", "polylog": "PolyLog",
    "lerch_phi": "LerchPhi", "lambert_w": "ProductLog", "bessel_j": "BesselJ",
    "bessel_y": "BesselY", "bessel_i": "BesselI", "bessel_k": "BesselK",
    "hankel_1": "HankelH1", "hankel_2": "HankelH2", "airy_ai": "AiryAi", "airy_bi": "AiryBi",
    "airy_ai_prime": "AiryAiPrime", "airy_bi_prime": "AiryBiPrime",
    "elliptic_k": "EllipticK", "elliptic_e": "EllipticE", "elliptic_f": "EllipticF",
    "elliptic_pi": "EllipticPi", "hypergeometric": "HypergeometricPFQ",
    "legendre_p": "LegendreP", "legendre_q": "LegendreQ", "hermite_h": "HermiteH",
    "laguerre_l": "LaguerreL", "chebyshev_t": "ChebyshevT", "chebyshev_u": "ChebyshevU",
    "gegenbauer_c": "GegenbauerC", "jacobi_p": "JacobiP", "heaviside": "HeavisideTheta",
    "dirac_delta": "DiracDelta", "kronecker_delta": "KroneckerDelta", "struve_h": "StruveH",
    "struve_l": "StruveL", "kelvin_ber": "KelvinBer", "kelvin_bei": "KelvinBei",
    # the four the runner lacks
    "harmonic": "HarmonicNumber",
}
for f in CAP:
    FUNCS[f] = f[0].upper() + f[1:]
for f in "sin cos tan cot sec csc sinh cosh tanh coth sech csch".split():
    FUNCS["a" + f] = "Arc" + f[0].upper() + f[1:]
SYMS = {"pi": "Pi", "e": "E", "i": "I", "infinity": "Infinity", "euler_gamma": "EulerGamma",
        "catalan": "Catalan", "undefined": "Indeterminate"}
SETS = {"integers": "Integers", "reals": "Reals", "rationals": "Rationals",
        "complexes": "Complexes", "naturals": "Integers"}


def ident(s):
    s = s.replace("_", "$")
    return "mr" + s if s in ("C", "D", "E", "I", "N", "O", "K") else s


TOK = re.compile(r"\s*(?:(\d+\.\d*|\d*\.\d+|\d+)|([A-Za-z_][A-Za-z_0-9]*)|(>=|<=|!=|[-+*/^(),=<>\[\]]))")


def tokenize(s):
    pos, out = 0, []
    s = s.strip()
    while pos < len(s):
        m = TOK.match(s, pos)
        if not m or m.end() == pos:
            raise ValueError("cannot tokenize %r at %d" % (s, pos))
        pos = m.end()
        num, name, op = m.groups()
        out.append(("num", num) if num else ("name", name) if name else ("op", op))
    out.append(("end", None))
    return out


class P:
    """Precedence parser of the neutral syntax, emitting Mathilda input text."""
    BIN = {"=": 1, "!=": 1, "<": 1, ">": 1, "<=": 1, ">=": 1, "+": 2, "-": 2, "*": 3, "/": 3, "^": 5}

    def __init__(self, s):
        self.t, self.i = tokenize(s), 0

    def peek(self):
        return self.t[self.i]

    def take(self):
        self.i += 1
        return self.t[self.i - 1]

    def expect(self, v):
        k, x = self.take()
        if x != v:
            raise ValueError("expected %s, got %s" % (v, x))

    def parse(self):
        e = self.expr(0)
        if self.peek()[0] != "end":
            raise ValueError("trailing input")
        return e

    def expr(self, minp):
        lhs = self.unary()
        while True:
            k, op = self.peek()
            if k != "op" or op not in self.BIN or self.BIN[op] < minp:
                return lhs
            p = self.BIN[op]
            self.take()
            rhs = self.expr(p if op == "^" else p + 1)
            wl = {"=": "==", "!=": "!="}.get(op, op)
            lhs = "(%s)%s(%s)" % (lhs, wl, rhs) if op == "^" else "%s %s %s" % (lhs, wl, rhs)
            if op in ("+", "-", "*", "/") or p == 1:
                lhs = "(" + lhs + ")"

    def unary(self):
        k, v = self.peek()
        if k == "op" and v == "-":
            self.take()
            return "(-(%s))" % self.expr(4)
        if k == "op" and v == "+":
            self.take()
            return self.expr(4)
        return self.atom()

    def args(self, close):
        out = []
        if self.peek() == ("op", close):
            self.take()
            return out
        while True:
            out.append(self.expr(0))
            k, v = self.take()
            if v == close:
                return out
            if v != ",":
                raise ValueError("bad argument list")

    def atom(self):
        k, v = self.take()
        if k == "num":
            return v if "." not in v else "%s`60" % v
        if k == "op" and v == "(":
            e = self.expr(0)
            self.expect(")")
            return "(" + e + ")"
        if k == "op" and v == "[":
            return "{" + ", ".join(self.args("]")) + "}"
        if k != "name":
            raise ValueError("unexpected %s" % v)
        if self.peek() != ("op", "("):
            return SYMS.get(v, SETS.get(v, ident(v)))
        self.take()
        a = self.args(")")
        return call(v, a)


def call(name, a):
    if name == "diff":
        f, rest = a[0], a[1:]
        vs, j = [], 0
        while j < len(rest):
            if j + 1 < len(rest) and re.fullmatch(r"\d+", rest[j + 1]):
                vs.append("{%s, %s}" % (rest[j], rest[j + 1])); j += 2
            else:
                vs.append(rest[j]); j += 1
        return "D[%s, %s]" % (f, ", ".join(vs))
    if name == "subs":        # subs(e, v = a): a == has been emitted
        m = re.fullmatch(r"\((.*) == (.*)\)", a[1])
        return "((%s) /. %s -> %s)" % (a[0], m.group(1), m.group(2))
    if name in ("integrate", "sum", "product") and len(a) == 4:
        h = {"integrate": "Integrate", "sum": "Sum", "product": "Product"}[name]
        return "%s[%s, {%s, %s, %s}]" % (h, a[0], a[1], a[2], a[3])
    if name == "element":
        return "Element[%s, %s]" % (a[0], a[1])
    if name == "root":
        return "((%s)^(1/(%s)))" % (a[0], a[1])
    if name == "atan2":
        return "ArcTan[%s, %s]" % (a[1], a[0])
    if name == "zeta":
        return ("Zeta[%s]" if len(a) == 1 else "HurwitzZeta[%s]") % ", ".join(a)
    if name == "digamma":
        return "PolyGamma[0, %s]" % a[0]
    if name == "root_of":     # root_of(p(_z), k): the k-th root, Root's numbering
        return "Root[(%s) &, %s]" % (re.sub(r"\$z\b", "#1", a[0]), a[1])
    return "%s[%s]" % (FUNCS.get(name, ident(name)), ", ".join(a))


def wl(s):
    return P(s).parse()


# ---------------------------------------------------------------- problems
def unknown_functions(prob):
    """The unknown function's head and variables, as `y(x)` / `u(x, t)`."""
    f = prob["args"]["function"]
    m = re.fullmatch(r"\s*([A-Za-z_]\w*)\((.*)\)\s*", f)
    return ident(m.group(1)), [ident(v.strip()) for v in m.group(2).split(",")]


def assumptions(prob):
    a = prob["args"]
    out = []
    for s in a.get("assumptions", []) or []:
        m = re.fullmatch(r"element\((\w+),\s*(\w+)\)", s.strip())
        if m and m.group(2) == "naturals":
            out.append("Element[%s, Integers] && %s >= 0" % (ident(m.group(1)), ident(m.group(1))))
        else:
            out.append(wl(s))
    if prob["op"] in ("integrate", "limit"):
        out.append("Element[%s, Reals]" % ident(a["var"]))
    if prob["op"] in ("sum", "product"):
        out.append("Element[%s, Integers]" % ident(a["var"]))
    return " && ".join(out) if out else "True"


def task(prob):
    op, a = prob["op"], prob["args"]
    asm = assumptions(prob)
    opt = "Assumptions -> " + asm
    if op == "integrate":
        f, x = wl(a["integrand"]), ident(a["var"])
        if "from" in a:
            return "Integrate[%s, {%s, %s, %s}, %s]" % (f, x, wl(a["from"]), wl(a["to"]), opt)
        return "Integrate[%s, %s, %s]" % (f, x, opt)
    if op in ("sum", "product"):
        return "%s[%s, {%s, %s, %s}, %s]" % ("Sum" if op == "sum" else "Product", wl(a["term"]),
                                             ident(a["var"]), wl(a["from"]), wl(a["to"]), opt)
    if op in ("solve_ode", "solve_pde", "solve_recurrence"):
        eqs = [wl(e) for e in (a["equations"] if "equations" in a else [a["equation"]])]
        eqs += [wl(c) for c in a.get("conditions", []) or []]
        f, vs = unknown_functions(prob)
        fx = "%s[%s]" % (f, ", ".join(vs))
        if op == "solve_recurrence":
            return "RSolve[{%s}, {%s}, %s]" % (", ".join(eqs), fx, vs[0])
        v = vs[0] if len(vs) == 1 else "{" + ", ".join(vs) + "}"
        return "DSolve[{%s}, {%s}, %s]" % (", ".join(eqs), fx, v)
    if op == "limit":
        pt = a["point"]
        d = ""
        if pt.endswith("+") or pt.endswith("-"):
            d = ', Direction -> "FromAbove"' if pt.endswith("+") else ', Direction -> "FromBelow"'
            pt = pt[:-1]
        return "Limit[%s, %s -> %s%s, %s]" % (wl(a["expr"]), ident(a["var"]), wl(pt), d, opt)
    if op == "solve":
        vs = ", ".join(ident(v) for v in a["vars"])
        eqs = [wl(e) for e in a["equations"]]
        ineq = any(re.search(r"[<>]|!=", e) for e in a["equations"])
        dom = "Reals" if a.get("domain") == "real" or ineq else "Complexes"
        f = "Reduce" if ineq else "Solve"
        if a.get("assumptions"):
            eqs.append(asm)
        return "%s[{%s}, {%s}, %s]" % (f, ", ".join(eqs), vs, dom)
    if op == "simplify":
        return "FullSimplify[%s, %s]" % (wl(a["expr"]), asm)
    if op == "evaluate":
        return "FunctionExpand[%s, %s]" % (wl(a["expr"]), asm)
    if op == "minimal_polynomial":
        return "MinimalPolynomial[%s, %s]" % (wl(a["expr"]), ident(a["var"]))
    if op == "factor":
        g = ", GaussianIntegers -> True" if a.get("domain") == "gaussian" else ""
        return "Factor[%s%s]" % (wl(a["expr"]), g)
    if op == "holds":
        vs = ", ".join(ident(v) for v in a["vars"])
        return "Resolve[ForAll[{%s}, %s, %s], Reals]" % (vs, asm, wl(a["relation"]))
    if op == "determinant" and "power" in a:
        return "Simplify[Det[MatrixPower[%s, %s]], %s]" % (wl(a["matrix"]), wl(a["power"]), opt)
    if op in ("determinant", "inverse", "eigenvalues"):
        h = {"determinant": "Det", "inverse": "Inverse", "eigenvalues": "Eigenvalues"}[op]
        return "%s[%s]" % (h, wl(a["matrix"]))
    if op == "laplace":
        return "LaplaceTransform[%s, %s, %s, %s]" % (wl(a["expr"]), ident(a["var"]),
                                                     ident(a["transform_var"]), opt)
    raise ValueError("op " + op)


TASK_HEAD = {"integrate": "Integrate", "sum": "Sum", "product": "Product", "limit": "Limit",
             "solve_ode": "DSolve", "solve_pde": "DSolve", "solve_recurrence": "RSolve",
             "solve": "Solve", "minimal_polynomial": "MinimalPolynomial", "holds": "Resolve",
             "laplace": "LaplaceTransform", "evaluate": None, "simplify": None, "factor": None,
             "determinant": "Det", "inverse": "Inverse", "eigenvalues": "Eigenvalues"}


# ---------------------------------------------------------------- sampling
def sample_values(prob, names):
    """Exact sample values for the free parameters, inside the assumptions."""
    lo, hi, integer = {}, {}, set()
    for s in prob["args"].get("assumptions", []) or []:
        s = s.strip()
        m = re.fullmatch(r"element\((\w+),\s*(\w+)\)", s)
        if m:
            if m.group(2) in ("integers", "naturals"):
                integer.add(m.group(1))
            if m.group(2) == "naturals":
                lo[m.group(1)] = max(lo.get(m.group(1), -1e9), 0)
            continue
        m = re.fullmatch(r"(\w+)\s*(>=|>|<=|<)\s*(.+)", s)
        if m:
            try:
                v = float(Fraction(m.group(3).replace(" ", "")))
            except (ValueError, ZeroDivisionError):
                continue
            if m.group(2)[0] == ">":
                lo[m.group(1)] = max(lo.get(m.group(1), -1e9), v)
            else:
                hi[m.group(1)] = min(hi.get(m.group(1), 1e9), v)
    pool = [Fraction(37, 100), Fraction(71, 100), Fraction(13, 10), Fraction(23, 100),
            Fraction(57, 100), Fraction(11, 10), Fraction(19, 10), Fraction(43, 100)]
    ints = [5, 4, 3, 6, 7, 2]
    out, k, j = {}, 0, 0
    for n in names:
        l, h = lo.get(n, -1e9), hi.get(n, 1e9)
        if n in integer:
            cand = [i for i in ints if l <= i <= h] or [int(l) + 1]
            out[n] = str(cand[j % len(cand)]); j += 1
            continue
        cand = [p for p in pool if l < p < h] or ([Fraction((Fraction(l) + Fraction(h)) / 2)]
                                                   if h < 1e8 and l > -1e8 else
                                                   [Fraction(l) + Fraction(1, 3)] if l > -1e8 else
                                                   [Fraction(h) - Fraction(1, 3)])
        p = cand[k % len(cand)]; k += 1
        out[n] = "%d/%d" % (p.numerator, p.denominator)
    return out


# ---------------------------------------------------------------- mathilda
def run_mathilda(binary, code, timeout):
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as f:
        f.write(code)
        path = f.name
    t0 = time.time()
    try:
        p = subprocess.run([binary, "-file", path], capture_output=True, text=True,
                           timeout=timeout)
        return p.returncode, p.stdout, p.stderr, time.time() - t0
    except subprocess.TimeoutExpired:
        return "timeout", "", "", timeout
    finally:
        os.unlink(path)


def solve_code(prob):
    return ('mrR = %s;\nWriteString["stdout", "\\n@@ANSWER ", ToString[mrR, InputForm, '
            'PageWidth -> Infinity], "\\n"];\n' % task(prob))


def wlstr(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def check_code(prob, answer):
    a = prob["args"]
    ref = prob["answer"] or ""
    head = TASK_HEAD.get(prob["op"])
    names = set()
    blob = json.dumps(a) + " " + ref
    for n in re.findall(r"\b([A-Za-z_]\w*)\b(?!\s*\()", blob):
        names.add(n)
    skip = set(SYMS) | set(SETS) | {"C1", "C2", "C3", "F1", "F2", "integers", "reals", "_z"}
    if prob["op"] == "solve":      # k numbers a reference solution family
        skip.add("k")
    if "var" in a:
        skip.add(a["var"])
    if "vars" in a:
        skip |= set(a["vars"])
    if "function" in a:
        skip |= {v.strip() for v in re.findall(r"\((.*)\)", a["function"])[0].split(",")}
        skip.add(re.match(r"\s*(\w+)", a["function"]).group(1))
    if "transform_var" in a:
        skip.add(a["transform_var"])
    params = sample_values(prob, sorted(n for n in names if n not in skip))
    prules = "{" + ", ".join("%s -> %s" % (ident(k), v) for k, v in params.items()) + "}"
    lines = ['Get[%s];' % wlstr(os.path.join(HERE, "check_famous.m")),
             'mrHeld = ToExpression[%s, InputForm, Hold];' % wlstr(answer),
             'mrAns = If[ospUneval[mrHeld, %s], mrHeld, ReleaseHold[mrHeld]];' % (head or "None"),
             'mrP = %s;' % prules,
             'mrHead = %s;' % (head or "None"),
             'mrOp = %s;' % wlstr(prob["op"]),
             'mrIndef = %s;' % (ident(a["var"]) if prob["op"] == "integrate" and "from" not in a else "False")]
    if not ref:
        lines.append('mrG = ospExpectUneval[mrAns, mrHead];')
    elif prob["op"] in ("solve_ode", "solve_pde", "solve_recurrence"):
        f, vs = unknown_functions(prob)
        eqs = [wl(e) for e in (a["equations"] if "equations" in a else [a["equation"]])]
        conds = [wl(c) for c in a.get("conditions", []) or []]
        lines += ['mrF = %s; mrVars = {%s};' % (f, ", ".join(vs)),
                  'mrEqs = Hold[{%s}]; mrConds = Hold[{%s}];' % (", ".join(eqs), ", ".join(conds)),
                  'mrRef = Hold[%s];' % wl(ref),
                  'mrG = ospCheckDE[mrOp, mrAns, mrF, mrVars, mrEqs, mrConds, mrRef, mrP, %d];'
                  % ode_order(prob)]
    elif prob["op"] == "solve":
        eqs = [wl(e) for e in a["equations"]]
        lines += ['mrVars = {%s};' % ", ".join(ident(v) for v in a["vars"]),
                  'mrG = ospCheckSolve[mrAns, Hold[{%s}], %s, mrVars, mrP];'
                  % (", ".join(eqs), wl(ref))]
    elif prob["op"] == "minimal_polynomial":
        lines.append('mrG = ospCheckMinPoly[mrAns, %s, %s];' % (wl(ref), ident(a["var"])))
    elif prob["op"] == "holds":
        lines.append('mrG = If[mrAns === True, "solved", If[mrAns === False, "wrong", "unsolved"]];')
    elif prob["op"] == "eigenvalues":
        lines.append('mrG = ospCheckMultiset[mrAns, %s, mrP];' % wl(ref))
    else:
        inp = wl(a["expr"]) if prob["op"] in ("simplify", "evaluate") else "None"
        lines.append('mrG = ospCheckValue[mrAns, %s, mrHead, mrP, %s, %s];'
                     % (wl(ref), "True" if prob["op"] in ("simplify", "evaluate") else "False",
                        "Hold[" + inp + "]"))
    lines.append('WriteString["stdout", "\\n@@GRADE ", mrG, "\\n"];')
    return "\n".join(lines) + "\n"


def ode_order(prob):
    a = prob["args"]
    eqs = a["equations"] if "equations" in a else [a["equation"]]
    best = 0
    for e in eqs:
        for m in re.finditer(r"diff\([^,]+\([^)]*\),((?:\s*\w+\s*,?)+?)\)", e):
            parts = [p.strip() for p in m.group(1).split(",") if p.strip()]
            n = sum(int(p) for p in parts if p.isdigit()) + sum(
                1 for i, p in enumerate(parts)
                if not p.isdigit() and not (i + 1 < len(parts) and parts[i + 1].isdigit()))
            best = max(best, n)
    if prob["op"] == "solve_recurrence":
        best = max([int(x) for e in eqs for x in re.findall(r"n \+ (\d+)", e)] + [1])
    return best


# ---------------------------------------------------------------- driver
def grade(binary, prob, answer, timeout, verbose):
    rc, out, err, _ = run_mathilda(binary, check_code(prob, answer), timeout)
    if rc == "timeout":
        return "unclear", "checker timed out"
    m = re.search(r"@@GRADE (.*)", out)
    if not m:
        return "unclear", "checker printed no grade (rc=%s) %s" % (rc, (err or out)[-200:])
    g = m.group(1).strip()
    return (g.split(":", 1)[0], g.split(":", 1)[1].strip() if ":" in g else "")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--binary", default=os.path.join(ROOT, "Mathilda"))
    ap.add_argument("--suite", default=os.path.join(HERE, "famous.jsonl"))
    ap.add_argument("--only", help="comma-separated ids (with or without famous-)")
    ap.add_argument("--timeout", type=float, default=60)
    ap.add_argument("--runner", help="grade the answers of a math-run file instead")
    ap.add_argument("--regrade", action="store_true",
                    help="re-grade the answers saved in OUT/math.jsonl (after a checker change)")
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    probs = [json.loads(l) for l in open(args.suite)]
    if args.only:
        want = {w if w.startswith("famous-") else "famous-" + w for w in args.only.split(",")}
        probs = [p for p in probs if p["id"] in want]
    runner = {}
    if args.regrade:            # saved answers, re-graded; solve failures kept as they were
        for l in open(os.path.join(args.out, "math.jsonl")):
            r = json.loads(l)
            st = r["grade"] if r["grade"] in ("timeout", "crash", "error") else "ok"
            runner[r["id"]] = {"status": st, "answer": r["answer"], "wall": r["seconds"]}
        args.runner = "regrade"
    if args.runner and not args.regrade:
        for l in open(args.runner):
            r = json.loads(l)
            if "id" in r:
                runner[r["id"]] = r
    base = {}
    for l in open(os.path.join(HERE, "baseline_v0.6222.tsv")):
        if not l.startswith("#"):
            k, v = l.rstrip("\n").split("\t")
            base[k] = v

    os.makedirs(args.out, exist_ok=True)
    tag = "math" if args.regrade else "runner" if args.runner else "math"
    rows = []
    for i, r in enumerate(probs, 1):
        pid, prob = r["id"], r["problem"]
        prob = dict(prob, answer=r["answer"])
        t = 0.0
        if args.runner:
            rr = runner.get(pid, {"status": "missing"})
            status, answer = rr["status"], rr.get("answer") or ""
            t = rr.get("wall", 0.0)
        else:
            try:
                code = solve_code(prob)
            except Exception as e:      # translation gap: our bug, report it
                rows.append((pid, "error", "translation: %s" % e, "", 0.0)); continue
            rc, out, err, t = run_mathilda(args.binary, code, args.timeout)
            m = re.search(r"@@ANSWER (.*)", out)
            status = ("timeout" if rc == "timeout" else "ok" if m else
                      "crash" if isinstance(rc, int) and rc < 0 or rc in (139, 134, 138) else "error")
            answer = m.group(1).strip() if m else (err or out)[-300:].strip()
        if status != "ok":
            g, why = ("timeout" if status == "timeout" else "crash" if status == "crash" else "error",
                      answer[:200])
        else:
            g, why = grade(args.binary, prob, answer, max(60, args.timeout), args.verbose)
        rows.append((pid, g, why, answer, t))
        if args.verbose or g != "solved":
            print("%3d/%d %-9s %6.1fs %s  %s" % (i, len(probs), g, t, pid[7:], why[:110]),
                  flush=True)
        if args.verbose:
            print("        answer:", answer[:300])

    with open(os.path.join(args.out, "%s.jsonl" % tag), "w") as f:
        for pid, g, why, ans, t in rows:
            f.write(json.dumps({"id": pid, "grade": g, "why": why, "answer": ans,
                                "seconds": round(t, 3), "baseline": base.get(pid)}) + "\n")
    from collections import Counter
    c = Counter(g for _, g, _, _, _ in rows)
    gained = [p for p, g, *_ in rows if g == "solved" and base.get(p) != "solved"]
    lost = [p for p, g, *_ in rows if g != "solved" and base.get(p) == "solved"]
    with open(os.path.join(args.out, "%s.md" % tag), "w") as f:
        f.write("# OSP famous problems, %s score\n\n" % tag)
        f.write("%d problems: %s\n\n" % (len(rows), ", ".join("%s %d" % kv for kv in c.most_common())))
        f.write("Against the official v0.6222 run: +%d solved, -%d lost.\n\n" % (len(gained), len(lost)))
        if lost:
            f.write("Lost: %s\n\n" % ", ".join(p[7:] for p in lost))
        f.write("| problem | grade | v0.6222 | s | note |\n|---|---|---|---|---|\n")
        for pid, g, why, ans, t in rows:
            f.write("| %s | %s | %s | %.1f | %s |\n" % (pid[7:], g, base.get(pid, ""), t,
                                                    why.replace("|", "/")[:120]))
    print("%s: %s; +%d / -%d vs v0.6222" % (tag, dict(c), len(gained), len(lost)))
    if lost:
        print("LOST:", ", ".join(p[7:] for p in lost))


if __name__ == "__main__":
    main()
