#!/usr/bin/env python3.11
"""Together / Cancel / Apart shootout harness.

For every (engine, case) pair: emit the engine's script, run it in its own
process with a timeout, parse the timing and the answer it prints, and verify
the answer independently of every engine:

  * value -- the returned terms sum to the input as a rational function
    (exactly; modulo p for Modulus; with I and Sqrt[2] reduced by their
    minimal polynomials);
  * Together / Cancel -- the answer is ONE fraction N/D with gcd(N, D) a
    constant (over Q, Q(i), Q(Sqrt[2]) or GF(p) as the case requires);
  * Apart -- every term is either free of the variable (polynomial part) or
    proper in it with a denominator that is a power of ONE irreducible (over
    the coefficient field; parameters are coefficients), and no two terms
    share a (base, exponent).  The number of variable-dependent terms is
    recorded.

Timing is taken inside each engine: one untimed warm-up, then the best of 3
(only the warm-up time when it exceeds 5 s).  The process is killed after a
wall limit.

Status: OK, CHECK-FAIL, TIMEOUT, UNSUPPORTED (Mathilda/reference returned the
call unevaluated), SKIP (the engine has no such operation), ERROR.

    python3 run.py [--engines mathilda,reference,maxima,sympy,pari]
                   [--only family[,family]] [--cases id,id] [--tag baseline]
                   [--timeout 60]

Results go to results/<tag>-<engine>.json; report.py turns them into tables.
The polynomial arithmetic and the two parsers are adapted from
../99-factor-shootout/run.py (kernels LG = log(x), EX = exp(x) added).
"""
import os
import sys

try:
    import sympy  # noqa: F401
except ImportError:                      # re-exec under the Python that has sympy
    os.execv("/usr/local/bin/python3.11", ["/usr/local/bin/python3.11"] + sys.argv)

import argparse
import json
import random
import re
import signal
import subprocess
import tempfile
import time

from fractions import Fraction

from sympy import I as SI, Poly, Symbol, factor_list, sqrt as ssqrt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MATHILDA = os.environ.get("MATHILDA", os.path.join(ROOT, "Mathilda"))
REFERENCE = "/Applications/Mathematica.app/Contents/MacOS/wolframscript"
PY = "/usr/local/bin/python3.11"
SCRATCH = os.environ.get("TMPDIR") or tempfile.gettempdir()
OPNAME = {"together": "Together", "cancel": "Cancel", "apart": "Apart"}


def is_alg(case):
    return "I" in case["input"] or "sqrt(" in case["input"]


# --------------------------------------------------------------------------
# Emitters.  Each returns (suffix, script) or None (SKIP).  Every engine
# prints one line
#   RES <TAB> <status> <TAB> <seconds> <TAB> <n1 ## d1 ;; n2 ## d2 ...>
# Mathilda and the reference write it to a file so no page width can break it.
# --------------------------------------------------------------------------

def to_wl(s):
    s = s.replace("sqrt(2)", "Sqrt[2]").replace("LG", "Log[x]").replace("EX", "Exp[x]")
    return s


def wl_opts(case):
    o = case["options"]
    if "modulus" in o:
        return f", Modulus -> {o['modulus']}"
    if o.get("extension") == "auto":
        return ", Extension -> Automatic"
    return ""


def wl_script(case, tlimit, outpath, engine):
    op = OPNAME[case["op"]]
    call = f"{op}[in0, {case['var']}{wl_opts(case)}]" if case["op"] == "apart" \
        else f"{op}[in0{wl_opts(case)}]"
    split = "If[Head[r] === Plus, List @@ r, {r}]" if case["op"] == "apart" else "{r}"
    # Globals only; ClearSystemCache[] and the assignment r = ... defeat memos.
    return ".m", f"""in0 = {to_wl(case['input'])};
tl = {tlimit};
run[] := (ClearSystemCache[]; First[AbsoluteTiming[r = TimeConstrained[{call}, tl, $Aborted];]]);
t1 = run[];
st = "OK"; best = t1;
If[r === $Aborted, st = "TIMEOUT"];
If[st === "OK" && Head[r] === {op}, st = "UNSUPPORTED"];
If[st === "OK" && t1 < 5, best = 10.^9; Do[t2 = run[]; If[t2 < best, best = t2], {{3}}]];
fs = If[st === "OK", StringRiffle[Map[ToString[FullForm[Numerator[#]]] <> " ## " <> ToString[FullForm[Denominator[#]]] &, {split}], " ;; "], ""];
h = OpenWrite["{outpath}"];
WriteString[h, "RES\\t" <> st <> "\\t" <> ToString[Round[best*10^9]] <> "e-9\\t" <> fs <> "\\n"];
Close[h];
"""


def gp_script(case, tlimit, outpath, engine):
    if case["op"] == "apart" or is_alg(case) or "LG" in case["vars"] or case.get("nested"):
        return None                       # no general Apart; no I / sqrt / kernels here
    o = case["options"]
    m = f"Mod(1, {o['modulus']})*" if "modulus" in o else ""
    if case["op"] == "together":
        setup = f"v = {m}[{', '.join(case['terms'])}];"
        call = "vecsum(v)"
    else:
        setup = f"nn = {m}({case['num']});\ndd = {m}({case['den']});"
        call = "nn/dd"
    lift = "lift" if m else ""
    return ".gp", f"""default(parisizemax, 2^31);
{setup}
F() = {call};
t0 = getabstime(); r = alarm({int(tlimit)}, F()); t1 = (getabstime() - t0) / 1000.;
if(type(r) == "t_ERROR", print("RES\\tTIMEOUT\\t", t1, "\\t"); quit);
\\\\ the clock ticks in ms: repeat each timed run until it spans >= 100 ms
best = t1;
if(t1 < 5, best = 10^9; for(k = 1, 3, n = 0; t0 = getabstime(); until(getabstime() - t0 >= 100, g = F(); n++); dt = (getabstime() - t0) / 1000. / n; if(dt < best, best = dt)));
print("RES\\tOK\\t", best, "\\t", Str({lift}(numerator(r))), " ## ", Str({lift}(denominator(r))));
quit;
"""


def to_maxima(s):
    return s.replace("I", "%i").replace("LG", "log(x)").replace("EX", "exp(x)")


def maxima_script(case, tlimit, outpath, engine):
    o = case["options"]
    call = f"partfrac(q, {case['var']})" if case["op"] == "apart" else "ratsimp(q)"
    split = ("if mapatom(e) then [e] else if op(e) = \"+\" then args(e) else [e]"
             if case["op"] == "apart" else "[e]")
    mod_on = f"modulus:{o['modulus']}$" if "modulus" in o else ""
    alg = "algebraic:true$" if is_alg(case) else ""
    return ".mac", f"""linel:1000000$
display2d:false$
{alg}
qqin : {to_maxima(case['input'])}$
FF(q) := {call}$
TL(e) := {split}$
{mod_on}
t0 : elapsed_real_time()$
r : FF(qqin)$
t1 : elapsed_real_time() - t0$
best : t1$
if t1 < 5 then (best : 10^9, for k : 1 thru 3 do (t0 : elapsed_real_time(), r2 : FF(qqin), dt : elapsed_real_time() - t0, if dt < best then best : dt))$
ss : simplode(map(lambda([t], concat(string(num(t)), " ## ", string(denom(t)))), TL(r)), " ;; ")$
modulus:false$
printf(true, "RES~aOK~a~a~a~a~%", ascii(9), ascii(9), float(best), ascii(9), ss)$
"""


def sympy_script(case, tlimit, outpath, engine):
    o = case["options"]
    if "modulus" in o and case["op"] == "apart":
        return None                       # sympy's apart has no modulus
    kw = ""
    if "modulus" in o:
        kw = f", modulus={o['modulus']}"
    elif o.get("extension") == "auto":
        kw = ", extension=True"
    if case["op"] == "apart":
        call, split = f"apart(q, {case['var']})", "Add.make_args(r)"
    else:
        call, split = f"cancel(q{kw})", "[r]"
    syms = " ".join(v for v in case["vars"] if v not in ("LG", "EX"))
    pieces = case["terms"] if case["op"] == "together" else None
    return ".py", f"""
import time
from sympy import Add, symbols, sympify, cancel, apart, fraction, sqrt, I, log, exp
from sympy.core.cache import clear_cache
gens = symbols({syms!r}, seq=True)
loc = {{str(g): g for g in gens}}
loc['x'] = symbols('x'); loc['LG'] = log(loc['x']); loc['EX'] = exp(loc['x'])
def poly_text(t):   # term by term: one big sympify of a long sum overflows the compiler
    return Add(*[sympify(u, locals=loc) for u in t.replace(' - ', ' + -').split(' + ')])
pieces = {pieces!r}
if pieces is not None:
    q = Add(*[sympify(t.replace('^', '**'), locals=loc) for t in pieces])
elif {case['op'] == 'cancel'!r}:
    q = poly_text({case.get('num', '').replace('^', '**')!r}) / poly_text({case.get('den', '').replace('^', '**')!r})
else:
    n, d = {case['input'].replace('^', '**')!r}[1:-1].split(')/(')
    q = poly_text(n) / poly_text(d)
F = lambda q: {call}
clear_cache(); t0 = time.perf_counter(); r = F(q); t1 = time.perf_counter() - t0
best = t1
if t1 < 5:
    best = 1e9
    for _ in range(3):
        clear_cache(); t0 = time.perf_counter(); F(q); best = min(best, time.perf_counter() - t0)
s = " ;; ".join("%s ## %s" % tuple(map(str, fraction(t))) for t in {split})
print("RES\\tOK\\t%.9f\\t%s" % (best, s))
"""


ENGINES = {
    "mathilda":  (wl_script,     lambda p: [MATHILDA, "-file", p]),
    "reference": (wl_script,     lambda p: [REFERENCE, "-file", p]),
    "maxima":    (maxima_script, lambda p: ["maxima", "--very-quiet", "-b", p]),
    "sympy":     (sympy_script,  lambda p: [PY, p]),
    "pari":      (gp_script,     lambda p: ["gp", "-q", "-f", p]),
}
SHARED_SCRIPT_ENGINES = {"mathilda", "reference"}


def run_case(engine, case, tlimit, wall):
    """Run one case in a fresh process; return (status, seconds, answer text)."""
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
# Parsing: both parsers build a small AST, evaluated three ways --
#   * ev_mod:   a number modulo VP at a random point (values, Schwartz-Zippel);
#   * ev_upoly: a univariate rational function over GF(VP) in one variable,
#               the others at random points (reducedness: gcd(N, D) = 1);
#   * ev_exact: an exact sparse polynomial pair (Apart denominators, Modulus).
# VP = 2305843009213693921 is prime, = 1 mod 8, so I and Sqrt[2] have images
# in GF(VP) (a ring homomorphism from Q(i, Sqrt[2])[vars]); LG = log(x) and
# EX = exp(x) are algebraically independent of x, so they are free variables.
# Nodes: ("n", Fraction) ("v", name) ("+", [..]) ("*", [..]) ("^", b, int)
#        ("inv", a) ("i",) ("s2",) ("exp", k)   (exp(k x), k an integer)
# --------------------------------------------------------------------------
VP = 2305843009213693921
VI = 583529827753931384               # sqrt(-1) mod VP
VS2 = 432741127753990578              # sqrt(2) mod VP


def ast_const(a):
    """Exact value of a numeric subtree (exponents), else None."""
    t = a[0]
    if t == "n":
        return a[1]
    if t in ("+", "*"):
        vals = [ast_const(u) for u in a[1]]
        if any(v is None for v in vals):
            return None
        out = Fraction(0) if t == "+" else Fraction(1)
        for v in vals:
            out = out + v if t == "+" else out * v
        return out
    if t == "inv":
        v = ast_const(a[1])
        return None if v is None else 1 / v
    if t == "^":
        v = ast_const(a[1])
        return None if v is None else v ** a[2]
    return None


def mk_pow(base, ev):
    """base^ev with ev a Fraction; n^(k/2) for n = 2 becomes Sqrt[2] powers."""
    if ev.denominator == 2:
        n = ast_const(base)
        if n != 2:
            raise ValueError("unsupported radical")
        k2 = ev.numerator
        out = ("*", [("n", Fraction(2) ** ((k2 - 1) // 2)), ("s2",)]) if k2 % 2 else ("n", Fraction(2) ** (k2 // 2))
        return out
    if ev.denominator != 1:
        raise ValueError("non-integer exponent")
    return ("^", base, int(ev))


def mk_exp(arg):
    """E^arg for arg = k*x."""
    k = None
    if arg == ("v", "x"):
        k = 1
    elif arg[0] == "*" and len(arg[1]) == 2 and arg[1][1] == ("v", "x") and arg[1][0][0] == "n":
        k = arg[1][0][1]
    elif arg[0] == "*" and len(arg[1]) == 2 and arg[1][0] == ("v", "x") and arg[1][1][0] == "n":
        k = arg[1][1][1]
    if k is None or Fraction(k).denominator != 1:
        raise ValueError(f"unsupported exponential {arg}")
    return ("exp", int(k))


INFIX = re.compile(r"\s*(?:(\d+)|(%?[A-Za-z_][A-Za-z0-9_]*)|(\*\*|[-+*/^(),]))")


def parse_infix(s):
    """Infix text (corpus / Maxima / sympy / PARI) -> AST."""
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
        t = term()
        items = [("*", [("n", Fraction(-1)), t]) if neg else t]
        while peek() in (("o", "+"), ("o", "-")):
            sign = take()[1]
            t = term()
            items.append(t if sign == "+" else ("*", [("n", Fraction(-1)), t]))
        return items[0] if len(items) == 1 else ("+", items)

    def term():
        items = [power()]
        while peek() in (("o", "*"), ("o", "/")):
            op = take()[1]
            rhs = power()
            items.append(rhs if op == "*" else ("inv", rhs))
        return items[0] if len(items) == 1 else ("*", items)

    def exponent():
        neg = False
        if peek() == ("o", "-"):
            take()
            neg = True
        if peek() == ("o", "("):
            take()
            e = expr()
            assert take() == ("o", ")")
        else:
            e = atom()
        return e, neg

    def power():
        if peek() == ("v", "%e"):
            take()
            assert take() == ("o", "^")
            e, neg = exponent()
            if neg:
                e = ("*", [("n", Fraction(-1)), e])
            return mk_exp(e)
        base = atom()
        if peek() == ("o", "^"):
            take()
            e, neg = exponent()
            ev = ast_const(e)
            if ev is None:
                raise ValueError("symbolic exponent")
            return mk_pow(base, -ev if neg else ev)
        return base

    def atom():
        t = take()
        if t[0] == "n":
            return ("n", Fraction(t[1]))
        if t[0] == "v":
            name = t[1]
            if peek() == ("o", "("):               # sqrt(2), log(x), exp(k*x)
                take()
                arg = expr()
                assert take() == ("o", ")")
                if name == "sqrt":
                    return mk_pow(arg, Fraction(1, 2))
                if name == "log":
                    if arg != ("v", "x"):
                        raise ValueError("unsupported log")
                    return ("v", "LG")
                if name == "exp":
                    return mk_exp(arg)
                raise ValueError(f"function {name}")
            if name in ("I", "%i"):
                return ("i",)
            return ("v", name)
        if t == ("o", "("):
            e = expr()
            assert take() == ("o", ")")
            return e
        if t == ("o", "-"):
            return ("*", [("n", Fraction(-1)), power()])
        raise ValueError(f"unexpected {t}")

    out = expr()
    if peek() != ("o", "END"):
        raise ValueError(f"trailing input at token {k}")
    return out


TOK = re.compile(r"\s*(?:(-?\d+)|([A-Za-z$][A-Za-z0-9$`]*)|(\[)|(\])|(,))")


def parse_fullform(s):
    """FullForm text -> AST.  Sqrt[2], I, Log[x], E^(k x) handled."""
    pos = 0

    def node():
        nonlocal pos
        m = TOK.match(s, pos)
        if not m:
            raise ValueError(f"FullForm parse error at {pos}: {s[pos:pos + 40]!r}")
        pos = m.end()
        if m.group(1):
            return ("n", Fraction(int(m.group(1))))
        name = m.group(2)
        m2 = TOK.match(s, pos)
        if not (m2 and m2.group(3)):
            return ("E",) if name == "E" else ("v", name)
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

    def build(h, a):
        if h == "Rational":
            return ("n", a[0][1] / a[1][1])
        if h == "Complex":
            return ("+", [a[0], ("*", [a[1], ("i",)])])
        if h == "Log":
            if a[0] != ("v", "x"):
                raise ValueError("unsupported Log")
            return ("v", "LG")
        if h == "Plus":
            return ("+", a)
        if h == "Times":
            return ("*", a)
        if h == "Power":
            b, e = a
            if b == ("E",):
                return mk_exp(e)
            ev = ast_const(e)
            if ev is None:
                raise ValueError("symbolic exponent")
            return mk_pow(b, ev)
        raise ValueError(f"unexpected head {h}")

    out = node()
    if out == ("E",):
        raise ValueError("bare E")
    return out


# ---- evaluation modulo VP ---------------------------------------------------

def inv_mod(a, p):
    a %= p
    if a == 0:
        raise ZeroDivisionError
    return pow(a, p - 2, p)


def ev_mod(a, env, p=VP):
    t = a[0]
    if t == "n":
        return a[1].numerator % p * inv_mod(a[1].denominator, p) % p
    if t == "v":
        return env[a[1]]
    if t == "+":
        return sum(ev_mod(u, env, p) for u in a[1]) % p
    if t == "*":
        out = 1
        for u in a[1]:
            out = out * ev_mod(u, env, p) % p
        return out
    if t == "^":
        v = ev_mod(a[1], env, p)
        return pow(v, a[2], p) if a[2] >= 0 else pow(inv_mod(v, p), -a[2], p)
    if t == "inv":
        return inv_mod(ev_mod(a[1], env, p), p)
    if t == "i":
        return VI
    if t == "s2":
        return VS2
    if t == "exp":
        v = env["EX"]
        return pow(v, a[1], p) if a[1] >= 0 else pow(inv_mod(v, p), -a[1], p)
    raise ValueError(f"node {t}")


# ---- univariate rational functions over GF(VP) -------------------------------
# A polynomial is a list of coefficients (low degree first), no trailing zeros.

def p_trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def p_add(a, b, p=VP):
    n = max(len(a), len(b))
    return p_trim([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p for i in range(n)])


def p_mul(a, b, p=VP):
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] = (out[i + j] + x * y) % p
    return p_trim(out)


def p_pow(a, e, p=VP):
    out, base = [1], a
    while e:
        if e & 1:
            out = p_mul(out, base, p)
        e >>= 1
        if e:
            base = p_mul(base, base, p)
    return out


def p_mod(a, b, p=VP):
    a = list(a)
    il = inv_mod(b[-1], p)
    db = len(b) - 1
    while len(a) - 1 >= db and a:
        c = a[-1] * il % p
        s = len(a) - 1 - db
        for i in range(db + 1):
            a[s + i] = (a[s + i] - c * b[i]) % p
        p_trim(a)
    return a


def p_gcd(a, b, p=VP):
    while b:
        a, b = b, p_mod(a, b, p)
    return a


def ev_upoly(a, v, env, p=VP):
    """(num, den) univariate in variable v over GF(p); other symbols from env."""
    t = a[0]
    if t == "v" and a[1] == v:
        return [0, 1], [1]
    if t == "+":
        n, d = [], [1]
        for u in a[1]:
            un, ud = ev_upoly(u, v, env, p)
            n, d = p_add(p_mul(n, ud, p), p_mul(un, d, p), p), p_mul(d, ud, p)
        return n, d
    if t == "*":
        n, d = [1], [1]
        for u in a[1]:
            un, ud = ev_upoly(u, v, env, p)
            n, d = p_mul(n, un, p), p_mul(d, ud, p)
        return n, d
    if t == "^":
        n, d = ev_upoly(a[1], v, env, p)
        if a[2] < 0:
            n, d = d, n
        return p_pow(n, abs(a[2]), p), p_pow(d, abs(a[2]), p)
    if t == "inv":
        n, d = ev_upoly(a[1], v, env, p)
        if not n:
            raise ZeroDivisionError
        return d, n
    if t == "exp" and v == "EX":
        k = a[1]
        mon = [0] * abs(k) + [1]
        return (mon, [1]) if k >= 0 else ([1], mon)
    c = ev_mod(a, env, p)
    return ([c] if c else []), [1]


# ---- exact sparse polynomials (Apart denominators, Modulus) ------------------

def ev_exact(a, names):
    """Exact (num, den) as sympy expressions (Sqrt[2], I kept algebraic)."""
    from sympy import Rational, Integer, exp as sexp
    t = a[0]
    if t == "n":
        return Rational(a[1].numerator, a[1].denominator), Integer(1)
    if t == "v":
        return Symbol(a[1]), Integer(1)
    if t == "+":
        n, d = Integer(0), Integer(1)
        for u in a[1]:
            un, ud = ev_exact(u, names)
            n, d = n * ud + un * d, d * ud
        return n, d
    if t == "*":
        n, d = Integer(1), Integer(1)
        for u in a[1]:
            un, ud = ev_exact(u, names)
            n, d = n * un, d * ud
        return n, d
    if t == "^":
        n, d = ev_exact(a[1], names)
        return (n ** a[2], d ** a[2]) if a[2] >= 0 else (d ** -a[2], n ** -a[2])
    if t == "inv":
        n, d = ev_exact(a[1], names)
        return d, n
    if t == "i":
        return SI, Integer(1)
    if t == "s2":
        return ssqrt(2), Integer(1)
    if t == "exp":
        return (Symbol("EX") ** a[1], Integer(1)) if a[1] >= 0 else (Integer(1), Symbol("EX") ** -a[1])
    raise ValueError(f"node {t}")


# --------------------------------------------------------------------------
# Verification
# --------------------------------------------------------------------------

def parse_terms(engine, text):
    out = []
    for item in text.split(" ;; "):
        if not item.strip():
            continue
        if " ## " not in item:
            raise ValueError(f"malformed term {item[:60]!r}")
        n, d = item.split(" ## ", 1)
        if engine in SHARED_SCRIPT_ENGINES:
            out.append((parse_fullform(n.strip()), parse_fullform(d.strip())))
        else:
            out.append((parse_infix(n.strip().replace("%i", "I")),
                        parse_infix(d.strip().replace("%i", "I"))))
    return out


def field_kw(case):
    o = case["options"]
    if "modulus" in o:
        return {"modulus": o["modulus"]}
    if "sqrt(" in case["input"]:
        return {"extension": ssqrt(2)}
    if "I" in case["input"]:
        return {"extension": SI}
    return {}


def rand_env(names, rng, p):
    env = {v: rng.randrange(2, p - 1) for v in names}
    env["EX"] = env.get("EX", rng.randrange(2, p - 1))
    return env


def check_value(case, inp, terms, rng):
    """Sum of terms == input, at 3 random points mod VP (exact for Modulus)."""
    o = case["options"]
    if "modulus" in o:
        p = o["modulus"]
        x = case["vars"][0]
        n_in, d_in = ev_exact(inp, case["vars"])
        num, den = 0, 1
        for n, d in terms:
            (nn, nd), (dn, dd) = ev_exact(n, case["vars"]), ev_exact(d, case["vars"])
            num, den = num * dn * nd + nn * dd * den, den * dn * nd
        diff = Poly((num * d_in - n_in * den).expand(), Symbol(x), modulus=p)
        return diff.is_zero
    for _ in range(3):
        while True:
            env = rand_env(case["vars"], rng, VP)
            try:
                lhs = ev_mod(inp, env)
                rhs = 0
                for n, d in terms:
                    rhs = (rhs + ev_mod(n, env) * inv_mod(ev_mod(d, env), VP)) % VP
                break
            except ZeroDivisionError:
                continue
        if lhs != rhs:
            return False
    return True


def check_reduced(case, n, d, rng):
    """gcd(N, D) = 1: in each variable v, univariate gcd of the images after a
    random substitution of the others (a true common factor survives it)."""
    o = case["options"]
    p = o.get("modulus", VP)
    names = list(case["vars"])
    for v in names:
        for _ in range(4):
            env = rand_env(names, rng, p)
            try:
                (a, b), (c, e) = ev_upoly(n, v, env, p), ev_upoly(d, v, env, p)
            except ZeroDivisionError:
                continue
            N, D = p_mul(a, e, p), p_mul(b, c, p)
            if not D:
                continue
            if len(p_gcd(N, D, p)) > 1:
                if p < 1000:          # small field: a spurious root is plausible, retry
                    continue
                return False
            break
    return True


def verify(case, text, engine):
    """(ok, term count, reason)."""
    rng = random.Random(hash(case["id"]) & 0xffffffff)
    terms = parse_terms(engine, text)
    if not terms:
        return False, 0, "no terms"
    inp = parse_infix(case["input"])
    if not check_value(case, inp, terms, rng):
        return False, None, "value differs"
    if case["op"] in ("together", "cancel"):
        if len(terms) != 1:
            return False, None, f"{len(terms)} terms, not one fraction"
        if not check_reduced(case, terms[0][0], terms[0][1], rng):
            return False, None, "not reduced"
        return True, 1, ""
    var = case["var"]
    gens = [Symbol(v) for v in case["vars"]]
    kw = field_kw(case)
    seen, count = set(), 0
    for n, d in terms:
        (nn, nd), (dn, dd) = ev_exact(n, case["vars"]), ev_exact(d, case["vars"])
        N, D = (nn * dd).expand(), (dn * nd).expand()
        if not D.has(Symbol(var)):
            continue                       # polynomial part
        if N.has(Symbol(var)) and Poly(N, Symbol(var)).degree() >= Poly(D, Symbol(var)).degree():
            return False, None, "improper term"
        _, fl = factor_list(D, *gens, **kw)
        dep = [(f, e) for f, e in fl if f.has(Symbol(var))]
        if len(dep) != 1:
            return False, None, f"denominator with {len(dep)} distinct factors in {var}"
        f, e = dep[0]
        key = (str(Poly(f, *gens, **kw).monic().as_expr()), e)
        if key in seen:
            return False, None, "repeated (base, exponent)"
        seen.add(key)
        count += 1
    return True, count, ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engines", default="mathilda,reference,maxima,sympy,pari")
    ap.add_argument("--only", default="", help="comma-separated families")
    ap.add_argument("--cases", default="", help="comma-separated case ids")
    ap.add_argument("--tag", default="run")
    ap.add_argument("--timeout", type=float, default=60, help="per-case limit, seconds")
    a = ap.parse_args()
    corpus = json.load(open(os.path.join(HERE, "corpus.json")))
    fams = set(a.only.split(",")) if a.only else None
    ids = set(a.cases.split(",")) if a.cases else None
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    wall = 1.5 * a.timeout + 30
    for eng in a.engines.split(","):
        outpath = os.path.join(HERE, "results", f"{a.tag}-{eng}.json")
        results = json.load(open(outpath)) if os.path.exists(outpath) else {}
        for c in corpus["cases"]:
            if (fams and c["family"] not in fams) or (ids and c["id"] not in ids):
                continue
            t0 = time.time()
            st, secs, txt = run_case(eng, c, a.timeout, wall)
            if st == "OK" and secs is not None and secs > a.timeout:
                st = "TIMEOUT"
            ok, count, err = None, None, ""
            if st == "OK":
                try:
                    ok, count, err = verify(c, txt, eng)
                    if not ok:
                        st = "CHECK-FAIL"
                except Exception as ex:
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
