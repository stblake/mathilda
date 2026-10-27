#!/usr/bin/env python3.11
"""Polynomial Factor shootout harness.

For every (engine, case) pair: emit the engine's script, run it in its own
process with a timeout, parse the timing and the factor list it prints, and
verify the answer independently of every engine (in sympy):

  * reconstruction -- the product of the returned factors^multiplicity equals
    the input exactly (modulo p for Modulus; modulo the minimal polynomials of
    the adjoined generators for GaussianIntegers / Extension; by
    cross-multiplication for a rational function);
  * the count of non-constant factors (sum of |multiplicity|) is recorded.
    report.py compares it against the largest count any verified engine found
    (PARI and the reference engine included): every verified factorisation's
    count is a lower bound on the true irreducible count, so a lower count is
    an under-factored answer.  Over-factoring cannot pass reconstruction.

Timing is taken inside each engine: one untimed warm-up, then the best of 3
(only the warm-up time when it exceeds 5 s).  The warm-up is bounded by an
in-engine time limit where the engine has one (TimeConstrained, alarm); the
process as a whole is killed after a wall limit.

Per-case status: OK, CHECK-FAIL (reconstruction failed), TIMEOUT,
UNSUPPORTED (Mathilda/reference returned Factor[...] unevaluated), SKIP (the
engine has no such operation), ERROR (crash / no parseable answer).

    python3 run.py [--engines mathilda,reference,pari,maxima,sympy]
                   [--only family[,family]] [--cases id,id] [--tag baseline]
                   [--timeout 60]

Results go to results/<tag>-<engine>.json; report.py turns them into tables.
"""
import os
import sys

try:
    import sympy  # noqa: F401
except ImportError:                      # re-exec under the Python that has sympy
    os.execv("/usr/local/bin/python3.11", ["/usr/local/bin/python3.11"] + sys.argv)

import argparse
import json
import re
import signal
import subprocess
import tempfile
import time

from fractions import Fraction

from sympy import Poly, QQ, Symbol, reduced, sympify

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MATHILDA = os.environ.get("MATHILDA", os.path.join(ROOT, "Mathilda"))
REFERENCE = "/Applications/Mathematica.app/Contents/MacOS/wolframscript"
PY = "/usr/local/bin/python3.11"
SCRATCH = os.environ.get("TMPDIR") or tempfile.gettempdir()

# Minimal polynomial of the field generator t that PARI (nffactor) and Maxima
# (factor(p, minpoly)) express their GaussianIntegers / Extension answers in.
MINPOLY_T = {(2,): "t^2-2", (2, 3): "t^4-10*t^2+1", "gauss": "t^2+1"}


def field_minpoly(opts):
    if opts.get("gaussian"):
        return MINPOLY_T["gauss"]
    if "extension" in opts:
        return MINPOLY_T[tuple(opts["extension"])]
    return None


# --------------------------------------------------------------------------
# Emitters.  Each returns (suffix, script).  Every engine prints one line
#   RES <TAB> <status> <TAB> <seconds> <TAB> <f1 @@ m1 ;; f2 @@ m2 ...>
# (status OK / TIMEOUT / UNSUPPORTED).  The engines sharing one script write it to a
# file instead of stdout, so no page-width line breaking can touch it.
# --------------------------------------------------------------------------

def wl_opts(opts):
    if "modulus" in opts:
        return f", Modulus -> {opts['modulus']}"
    if opts.get("gaussian"):
        return ", GaussianIntegers -> True"
    if "extension" in opts:
        rs = [f"Sqrt[{r}]" for r in opts["extension"]]
        return ", Extension -> " + (rs[0] if len(rs) == 1 else "{" + ", ".join(rs) + "}")
    return ""


def wl_script(case, tlimit, outpath, engine):
    # Globals only (a Module-local Do/Table count is a known Mathilda quirk).
    # ClearSystemCache[] and the assignment r = ... defeat result memos, so
    # each timed run really factors.
    return ".m", f"""in0 = {case['input']};
tl = {tlimit};
run[] := (ClearSystemCache[]; First[AbsoluteTiming[r = TimeConstrained[Factor[in0{wl_opts(case['options'])}], tl, $Aborted];]]);
t1 = run[];
st = "OK"; best = t1;
If[r === $Aborted, st = "TIMEOUT"];
If[st === "OK" && Head[r] === Factor, st = "UNSUPPORTED"];
If[st === "OK" && t1 < 5, best = 10.^9; Do[t2 = run[]; If[t2 < best, best = t2], {{3}}]];
fl[e_] := Map[If[Head[#] === Power && IntegerQ[#[[2]]], {{#[[1]], #[[2]]}}, {{#, 1}}] &, If[Head[e] === Times, List @@ e, {{e}}]];
fs = If[st === "OK", StringRiffle[Map[ToString[FullForm[#[[1]]]] <> " @@ " <> ToString[#[[2]]] &, fl[r]], " ;; "], ""];
h = OpenWrite["{outpath}"];
WriteString[h, "RES\\t" <> st <> "\\t" <> ToString[Round[best*10^9]] <> "e-9\\t" <> fs <> "\\n"];
Close[h];
"""


def gp_script(case, tlimit, outpath, engine):
    o = case["options"]
    mp = field_minpoly(o)
    if mp and len(case["vars"]) > 1:
        return None                       # PARI's nffactor is univariate only
    setup, call, cnum = "", "factor(q)", "p"
    if "modulus" in o:
        call, cnum = f"factor(Mod(1, {o['modulus']})*q)", f"Mod(1, {o['modulus']})*p"
    elif mp:
        setup, call = f"nf = nfinit({mp});", "nffactor(nf, q)"
    # gp reads its input line by line: keep each statement on one line.
    return ".gp", f"""default(parisizemax, 2^31);
p = {case['input']};
{setup}
F(q) = {call};
fl(f) = my(v = List()); for(i = 1, #f~, listput(v, Str(lift(f[i,1]), " @@ ", f[i,2]))); strjoin(Vec(v), " ;; ");
t0 = getabstime(); f = alarm({int(tlimit)}, F(p)); t1 = (getabstime() - t0) / 1000.;
if(type(f) == "t_ERROR", print("RES\\tTIMEOUT\\t", t1, "\\t"); quit);
best = t1; reps = if(t1 < 0.02, 50, 1);
if(t1 < 5, best = 10^9; for(k = 1, 3, t0 = getabstime(); for(j = 1, reps, g = F(p)); dt = (getabstime() - t0) / 1000. / reps; if(dt < best, best = dt)));
c = lift(simplify(({cnum}) / factorback(f)));
s = fl(f); if(c != 1, s = Str(c, " @@ 1 ;; ", s));
print("RES\\tOK\\t", best, "\\t", s);
quit;
"""


MAXIMA_FL = r"""fl(e) := if mapatom(e) then [[e,1]] else if op(e) = "*" then apply(append, map(fl1, args(e))) else fl1(e)$
fl1(e) := if mapatom(e) then [[e,1]] else if op(e) = "^" and integerp(part(e,2)) then [[part(e,1), part(e,2)]] else if op(e) = "-" and length(args(e)) = 1 then append([[-1,1]], fl(part(e,1))) else if op(e) = "/" then append(fl(num(e)), map(lambda([q], [q[1], -q[2]]), fl(denom(e)))) else [[e,1]]$
fls(e) := simplode(map(lambda([q], concat(string(q[1]), " @@ ", string(q[2]))), fl(e)), " ;; ")$"""


def maxima_script(case, tlimit, outpath, engine):
    o = case["options"]
    if "modulus" in o and o["modulus"] > 2 ** 31:
        return None                       # Maxima's modular factor hangs at 60-bit p
    mp = field_minpoly(o)
    if o.get("gaussian"):
        call = "gfactor(q)"
    elif mp:
        call = f"factor(q, {mp.replace('t', 'a')})"
    else:
        call = "factor(q)"
    mod_on = f"modulus:{o['modulus']}$" if "modulus" in o else ""
    return ".mac", f"""linel:1000000$
display2d:false$
factor_max_degree:1000000$
{MAXIMA_FL}
p : {case['input']}$
F(q) := {call}$
{mod_on}
t0 : elapsed_real_time()$
r : F(p)$
t1 : elapsed_real_time() - t0$
best : t1$
if t1 < 5 then (best : 10^9, for k : 1 thru 3 do (t0 : elapsed_real_time(), r2 : F(p), dt : elapsed_real_time() - t0, if dt < best then best : dt))$
modulus:false$
printf(true, "RES~aOK~a~a~a~a~%", ascii(9), ascii(9), float(best), ascii(9), fls(r))$
"""


def sympy_script(case, tlimit, outpath, engine):
    o = case["options"]
    if "modulus" in o and len(case["vars"]) > 1:
        return None                       # sympy: multivariate over GF(p) not implemented
    kw = ""
    if "modulus" in o:
        kw = f", modulus={o['modulus']}"
    elif o.get("gaussian"):
        kw = ", gaussian=True"
    elif "extension" in o:
        kw = ", extension=[" + ", ".join(f"sqrt({r})" for r in o["extension"]) + "]"
    syms = " ".join(case["vars"])
    if "/" in case["input"]:              # rational function: cancel, then factor num/den
        call = f"factor_list(cancel(q), *gens{kw}, frac=True)"
        unpack = "c, nf, df = r; fs = list(nf) + [(f, -m) for f, m in df]"
    else:
        call = f"factor_list(q, *gens{kw})"
        unpack = "c, fs = r"
    return ".py", f"""
import time
from sympy import Add, symbols, sympify, factor_list, sqrt, cancel
from sympy.core.cache import clear_cache
gens = symbols({syms!r}, seq=True)
loc = {{str(g): g for g in gens}}
def poly_text(t):   # term by term: one big sympify of a 5000-term sum overflows the compiler
    return Add(*[sympify(u, locals=loc) for u in t.replace(' - ', ' + -').split(' + ')])
src = {case['input'].replace('^', '**')!r}
if src.startswith('('):
    n, d = src[1:-1].split(')/(')
    p = poly_text(n) / poly_text(d)
else:
    p = poly_text(src)
F = lambda q: {call}
clear_cache(); t0 = time.perf_counter(); r = F(p); t1 = time.perf_counter() - t0
best = t1
if t1 < 5:
    best = 1e9
    for _ in range(3):
        clear_cache(); t0 = time.perf_counter(); F(p); best = min(best, time.perf_counter() - t0)
{unpack}
s = " ;; ".join([f"{{c}} @@ 1"] + [f"{{f}} @@ {{m}}" for f, m in fs])
print("RES\\tOK\\t%.9f\\t%s" % (best, s))
"""


ENGINES = {
    "mathilda":  (wl_script,     lambda p: [MATHILDA, "-file", p]),
    "reference": (wl_script,     lambda p: [REFERENCE, "-file", p]),
    "pari":      (gp_script,     lambda p: ["gp", "-q", "-f", p]),
    "maxima":    (maxima_script, lambda p: ["maxima", "--very-quiet", "-b", p]),
    "sympy":     (sympy_script,  lambda p: [PY, p]),
}
SHARED_SCRIPT_ENGINES = {"mathilda", "reference"}


def run_case(engine, case, tlimit, wall):
    """Run one case in a fresh process; return (status, seconds, factor text)."""
    emit, cmd = ENGINES[engine]
    fd, outpath = tempfile.mkstemp(suffix=".out", dir=SCRATCH)
    os.close(fd)
    em = emit(case, tlimit, outpath, engine)
    if em is None:
        os.unlink(outpath)
        return "SKIP", None, ""
    suffix, text = em
    fd, path = tempfile.mkstemp(suffix=suffix, dir=SCRATCH)
    with os.fdopen(fd, "w") as f:
        f.write(text)
    status, out = None, ""
    try:
        # A new session so a timeout kills the whole tree (maxima's shell
        # wrapper leaves its Lisp child behind otherwise).
        pr = subprocess.Popen(cmd(path), stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL, text=True, start_new_session=True,
                              cwd=SCRATCH)
        try:
            out, _ = pr.communicate(timeout=wall)
        except subprocess.TimeoutExpired:
            os.killpg(pr.pid, signal.SIGKILL)
            pr.communicate()
            status = "TIMEOUT"
        if engine in SHARED_SCRIPT_ENGINES and os.path.exists(outpath):
            out = open(outpath).read()
    finally:
        os.unlink(path)
        if os.path.exists(outpath):
            os.unlink(outpath)
    if status == "TIMEOUT":
        return "TIMEOUT", None, ""
    for line in out.splitlines():
        if line.startswith("RES\t"):
            parts = line.split("\t", 3) + [""] * 4
            st, t, fs = parts[1], parts[2].strip(), parts[3].strip()
            try:
                secs = float(t.replace(" ", "").replace("`", "").replace("*^", "e"))
            except ValueError:
                secs = None
            return st, secs, fs
    return "ERROR", None, ""


# --------------------------------------------------------------------------
# Parsing engine output, then independent verification.
#
# Polynomials are sparse dicts {exponent tuple: Fraction} over the generators
# (case variables..., _t, _s2, _s3, _i): _t is the PARI/Maxima field generator,
# _s2/_s3 are Sqrt[2]/Sqrt[3] and _i is I.  Both parsers below build these
# dicts directly -- building a sympy Add term by term is quadratic and took
# minutes on a 6000-term input.  A value is a pair (num, den) so the
# rational-function inputs parse too; den is the constant 1 almost always.
# --------------------------------------------------------------------------
ALG = ["_t", "_s2", "_s3", "_i"]


class PolyRing:
    def __init__(self, names):
        self.names = list(names) + ALG
        self.idx = {n: i for i, n in enumerate(self.names)}
        self.n = len(self.names)
        self.zero_mon = (0,) * self.n
        self.one = {self.zero_mon: Fraction(1)}

    def const(self, c):
        return {self.zero_mon: Fraction(c)} if c else {}

    def var(self, name, power=1):
        m = [0] * self.n
        m[self.idx[name]] = power
        return {tuple(m): Fraction(1)}

    @staticmethod
    def add(a, b, sign=1):
        out = dict(a)
        for m, c in b.items():
            v = out.get(m, 0) + sign * c
            if v:
                out[m] = v
            else:
                out.pop(m, None)
        return out

    @staticmethod
    def mul(a, b):
        if len(a) < len(b):
            a, b = b, a
        out = {}
        for mb, cb in b.items():
            for ma, ca in a.items():
                m = tuple(x + y for x, y in zip(ma, mb))
                v = out.get(m, 0) + ca * cb
                if v:
                    out[m] = v
                else:
                    out.pop(m, None)
        return out

    def is_const(self, a):
        return all(m == self.zero_mon for m in a)

    def cval(self, a):
        return a.get(self.zero_mon, Fraction(0))

    def pow(self, a, e):
        out, base = self.one, a
        while e:
            if e & 1:
                out = self.mul(out, base)
            e >>= 1
            if e:
                base = self.mul(base, base)
        return out

    # rational values (num, den)
    def r_add(self, x, y, sign=1):
        (a, b), (c, d) = x, y
        if b == self.one and d == self.one:
            return self.add(a, c, sign), b
        return self.add(self.mul(a, d), self.mul(c, b), sign), self.mul(b, d)

    def r_mul(self, x, y):
        return self.mul(x[0], y[0]), (self.one if x[1] == self.one and y[1] == self.one
                                      else self.mul(x[1], y[1]))

    def r_inv(self, x):
        a, b = x
        if self.is_const(a):
            c = self.cval(a)
            if c == 0:
                raise ZeroDivisionError("division by zero")
            return {m: v / c for m, v in b.items()}, self.one
        return b, a

    def r_pow(self, x, e):
        if e < 0:
            x, e = self.r_inv(x), -e
        return self.pow(x[0], e), (self.one if x[1] == self.one else self.pow(x[1], e))

    def sqrt_const(self, n):
        """Sqrt[n] for n in {2, 3, 6} in the generator encoding."""
        if n == 2:
            return self.var("_s2")
        if n == 3:
            return self.var("_s3")
        if n == 6:
            return self.mul(self.var("_s2"), self.var("_s3"))
        raise ValueError(f"Sqrt[{n}] not supported")


INFIX = re.compile(r"\s*(?:(\d+)|(%?[A-Za-z_][A-Za-z0-9_]*)|(\*\*|[-+*/^(),]))")


def parse_infix(s, R):
    """Infix text (PARI / Maxima / sympy / corpus) -> (num, den)."""
    toks = []
    pos = 0
    s = s.strip()
    while pos < len(s):
        m = INFIX.match(s, pos)
        if not m or m.end() == pos:
            raise ValueError(f"bad token at {s[pos:pos + 30]!r}")
        pos = m.end()
        if m.group(1):
            toks.append(("n", int(m.group(1))))
        elif m.group(2):
            toks.append(("v", m.group(2)))
        else:
            toks.append(("o", "^" if m.group(3) == "**" else m.group(3)))
    toks.append(("o", "END"))
    k = 0

    def peek():
        return toks[k]

    def take():
        nonlocal k
        k += 1
        return toks[k - 1]

    def expr():
        neg = False
        if peek() == ("o", "-"):
            take()
            neg = True
        elif peek() == ("o", "+"):
            take()
        acc = term()
        if neg:
            acc = ({m: -c for m, c in acc[0].items()}, acc[1])
        while peek() in (("o", "+"), ("o", "-")):
            sign = 1 if take()[1] == "+" else -1
            acc = R.r_add(acc, term(), sign)
        return acc

    def term():
        acc = power()
        while peek() in (("o", "*"), ("o", "/")):
            op = take()[1]
            rhs = power()
            acc = R.r_mul(acc, rhs if op == "*" else R.r_inv(rhs))
        return acc

    def power():
        base = atom()
        if peek() == ("o", "^"):
            take()
            neg = False
            if peek() == ("o", "-"):
                take()
                neg = True
            if peek() == ("o", "("):           # x^(2) / x^(-1)
                take()
                e = expr()
                assert take() == ("o", ")")
                ev = R.cval(e[0]) / R.cval(e[1])
            else:
                ev = Fraction(take()[1])
            if neg:
                ev = -ev
            if ev.denominator == 2 and R.is_const(base[0]) and base[1] == R.one:
                n = R.cval(base[0])
                k2 = ev.numerator               # n^(k2/2)
                root = R.sqrt_const(int(n))
                val = R.mul(R.const(Fraction(n) ** ((k2 - 1) // 2)), root) if k2 % 2 \
                    else R.const(Fraction(n) ** (k2 // 2))
                return val, R.one
            if ev.denominator != 1:
                raise ValueError("non-integer exponent")
            return R.r_pow(base, int(ev))
        return base

    def atom():
        t = take()
        if t[0] == "n":
            return R.const(t[1]), R.one
        if t[0] == "v":
            name = t[1]
            if peek() == ("o", "("):               # a function call: sqrt(n)
                take()
                arg = expr()
                assert take() == ("o", ")")
                if name != "sqrt":
                    raise ValueError(f"function {name}")
                return R.sqrt_const(int(R.cval(arg[0]))), R.one
            if name in ("I", "%i"):
                return R.var("_i"), R.one
            if name in ("t", "a"):
                return R.var("_t"), R.one
            if name not in R.idx:
                raise ValueError(f"unknown symbol {name}")
            return R.var(name), R.one
        if t == ("o", "("):
            e = expr()
            assert take() == ("o", ")")
            return e
        if t == ("o", "-"):
            v = power()
            return {m: -c for m, c in v[0].items()}, v[1]
        raise ValueError(f"unexpected {t}")

    out = expr()
    if peek() != ("o", "END"):
        raise ValueError(f"trailing input at token {k}")
    return out


TOK = re.compile(r"\s*(?:(-?\d+)|([A-Za-z$][A-Za-z0-9$`]*)|(\[)|(\])|(,))")


def parse_fullform(s, R):
    """FullForm text (Plus/Times/Power/Rational/Complex/Integer/Symbol) ->
    (num, den).  Sqrt[2], Sqrt[3], Sqrt[6] and I become generators."""
    pos = 0

    def node():
        nonlocal pos
        m = TOK.match(s, pos)
        if not m:
            raise ValueError(f"FullForm parse error at {pos}: {s[pos:pos + 40]!r}")
        pos = m.end()
        if m.group(1):
            return ("num", Fraction(int(m.group(1))))
        name = m.group(2)
        m2 = TOK.match(s, pos)
        if not (m2 and m2.group(3)):
            if name not in R.idx:
                raise ValueError(f"unknown symbol {name}")
            return (R.var(name), R.one)
        pos = m2.end()
        args = []
        while True:
            m3 = TOK.match(s, pos)
            if m3 and m3.group(4):
                pos = m3.end()
                break
            args.append(node())
            m3 = TOK.match(s, pos)
            if m3 and m3.group(5):
                pos = m3.end()
        return build(name, args)

    def val(x):
        return (R.const(x[1]), R.one) if x[0] == "num" else x

    def build(h, a):
        if h == "Rational":
            return ("num", a[0][1] / a[1][1])
        if h == "Complex":
            return R.r_add(val(a[0]), R.r_mul(val(a[1]), (R.var("_i"), R.one)))
        if h == "Plus":
            acc = val(a[0])
            for t in a[1:]:
                acc = R.r_add(acc, val(t))
            return acc
        if h == "Times":
            acc = val(a[0])
            for t in a[1:]:
                acc = R.r_mul(acc, val(t))
            return acc
        if h == "Power":
            b, e = a
            if e[0] != "num":
                raise ValueError("symbolic exponent")
            ev = e[1]
            if ev.denominator == 2 and b[0] == "num":
                n, k2 = int(b[1]), ev.numerator
                root = R.sqrt_const(n)
                v = R.mul(R.const(Fraction(n) ** ((k2 - 1) // 2)), root) if k2 % 2 \
                    else R.const(Fraction(n) ** (k2 // 2))
                return (v, R.one)
            if ev.denominator != 1:
                raise ValueError("non-integer exponent")
            return R.r_pow(val(b), int(ev))
        raise ValueError(f"unexpected head {h}")

    return val(node())


def parse_factors(engine, text, R):
    """Factor text -> [((num, den), multiplicity)]."""
    out = []
    for item in text.split(" ;; "):
        if not item.strip():
            continue
        if " @@ " not in item:
            raise ValueError(f"malformed factor item {item[:60]!r}")
        f, m = item.rsplit(" @@ ", 1)
        if engine in SHARED_SCRIPT_ENGINES:
            v = parse_fullform(f.strip(), R)
        else:
            v = parse_infix(f.strip().replace("%i", "I"), R)
        out.append((v, int(m)))
    return out


def verify(case, text, engine):
    """(reconstruction ok, non-constant factor count)."""
    o = case["options"]
    R = PolyRing(case["vars"])
    nv = len(case["vars"])
    factors = parse_factors(engine, text, R)
    num, den = R.one, R.one
    net = {}                          # monic-normalised factor -> net multiplicity
    for (a, b), m in factors:
        if b != R.one:
            if not R.is_const(b):
                raise ValueError("factor with a non-constant denominator")
            a = {k: v / R.cval(b) for k, v in a.items()}
        if any(any(mon[:nv]) for mon in a):
            lc = a[max(a)]
            key = tuple(sorted((k, v / lc) for k, v in a.items()))
            net[key] = net.get(key, 0) + m
        if m >= 0:
            num = R.mul(num, R.pow(a, m)) if m else num
        else:
            den = R.mul(den, R.pow(a, -m))
    count = sum(abs(m) for m in net.values())
    n_in, d_in = parse_infix(case["input"], R)
    gens = [Symbol(n) for n in R.names]
    if "modulus" in o:
        p = o["modulus"]
        lhs, rhs = R.mul(num, d_in), R.mul(n_in, den)
        diff = R.add(lhs, rhs, -1)
        if any(any(mon[nv:]) for mon in diff):
            return False, count
        ok = all(c.denominator % p and c.numerator % p == 0 for c in diff.values())
        ok = ok and all(c.denominator % p for c in list(num.values()) + list(den.values()))
        return ok, count
    P = lambda d: Poly.from_dict(d or {R.zero_mon: 0}, *gens, domain=QQ)
    diff = P(R.mul(num, d_in)) - P(R.mul(n_in, den))
    if diff.is_zero:
        return True, count
    # Reduce by the minimal polynomials of the generators that occur.
    rels = []
    if "extension" in o or o.get("gaussian"):
        mp = field_minpoly(o).replace("t", "_t").replace("^", "**")
        rels.append(sympify(mp, locals={"_t": Symbol("_t")}))
    rels += [Symbol("_s2") ** 2 - 2, Symbol("_s3") ** 2 - 3, Symbol("_i") ** 2 + 1]
    _, rem = reduced(diff.as_expr(), rels, *gens, domain=QQ)
    return rem == 0, count


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engines", default="mathilda,reference,pari,maxima,sympy")
    ap.add_argument("--only", default="", help="comma-separated families")
    ap.add_argument("--cases", default="", help="comma-separated case ids")
    ap.add_argument("--tag", default="run")
    ap.add_argument("--timeout", type=float, default=60, help="per-case limit, seconds")
    a = ap.parse_args()
    corpus = json.load(open(os.path.join(HERE, "corpus.json")))
    fams = set(a.only.split(",")) if a.only else None
    ids = set(a.cases.split(",")) if a.cases else None
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    wall = 1.5 * a.timeout + 30           # warm-up + 3 runs (each < 5 s) + start-up
    for eng in a.engines.split(","):
        outpath = os.path.join(HERE, "results", f"{a.tag}-{eng}.json")
        results = json.load(open(outpath)) if os.path.exists(outpath) else {}
        for c in corpus["cases"]:
            if (fams and c["family"] not in fams) or (ids and c["id"] not in ids):
                continue
            t0 = time.time()
            st, secs, txt = run_case(eng, c, a.timeout, wall)
            if st == "OK" and secs is not None and secs > a.timeout:
                st = "TIMEOUT"                   # engines with no in-engine limit
            ok, count, err = None, None, ""
            if st == "OK":
                try:
                    ok, count = verify(c, txt, eng)
                    if not ok:
                        st = "CHECK-FAIL"
                except Exception as ex:          # an answer the parser cannot read
                    st, err = "ERROR", f"{type(ex).__name__}: {ex}"[:300]
            if st not in ("OK", "UNSUPPORTED"):
                secs = None
            rec = {"status": st, "t": secs, "count": count, "wall": round(time.time() - t0, 2)}
            if err:
                rec["error"] = err
            if st in ("CHECK-FAIL", "ERROR") and len(txt) < 4000:
                rec["answer"] = txt
            results[c["id"]] = rec
            print(f"{eng:10s} {c['id']:22s} {st:12s} t={secs} n={count} {err}", flush=True)
            json.dump(results, open(outpath, "w"), indent=1)


if __name__ == "__main__":
    main()
