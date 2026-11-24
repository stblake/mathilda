"""crosscheck.py -- re-check, in Mathilda, every answer that failed its own system's check.

Each system first verifies its own answer (adapters.py: D[F] - f at three points
to 10^-20, plus FriCAS/Axiom's exact normalize).  An answer that passes stays
verified, and nothing more is done with it.  An answer that FAILS its own check
(a "bad" point) or could not be checked ("na" at every point) is not called
wrong on that evidence alone: the point check is the system's own arithmetic,
and Mathilda's answers are judged by the generator's test, PossibleZeroQ[D[F, x]
- f].  So the answer text is translated into Mathilda and judged by that same
test, giving every system the same criterion for "wrong":

    zero       PossibleZeroQ says D[F, x] - f == 0: the answer is correct
               (the system's own check failed on a branch, or could not evaluate)
    nonzero    PossibleZeroQ says it is not zero: WRONG
    undecided  PossibleZeroQ timed out, the answer would not parse in Mathilda, or
               D[F, x] - f does not evaluate to a number at two test points (then
               PossibleZeroQ's True is not evidence: it answers True for an
               expression it cannot evaluate, e.g. x b - b with b a Root over
               algebraic coefficients)
    untranslatable  the answer uses a construct outside to_wl (recorded, not guessed)

The last two are counted UNVERIFIED, never WRONG.

to_wl parses the infix output of FriCAS / Axiom (unparse of InputForm), Maxima
(string()) and SymPy (str()) with one precedence-climbing parser and prints
Mathilda InputForm.  rootOf(p(%%H0), %%H0) -- FriCAS's generic root, which the
answer is valid for whichever root is taken -- becomes Root[p &, 1]; SymPy's
RootSum(poly, Lambda(t, e)) becomes RootSum[poly &, e &].
"""
import hashlib
import os
import re

import adapters as A

PZQ_TIMEOUT = int(os.environ.get("SHOOTOUT_PZQ_TIMEOUT", "60"))


class Untranslatable(Exception):
    pass


TOK = re.compile(r"\s*(?:(\d+\.\d*(?:[eEb][-+]?\d+)?|\d+)|(%%?[A-Za-z][A-Za-z0-9]*|[A-Za-z_][A-Za-z0-9_]*)"
                 r"|(\*\*|[-+*/^(),\[\]=<>]))")

# function name -> Mathilda head, per source system (all share most names)
FUNCS = {
    "log": "Log", "exp": "Exp", "sqrt": "Sqrt", "abs": "Abs", "Abs": "Abs",
    "sin": "Sin", "cos": "Cos", "tan": "Tan", "cot": "Cot", "sec": "Sec", "csc": "Csc",
    "asin": "ArcSin", "acos": "ArcCos", "atan": "ArcTan", "acot": "ArcCot", "asec": "ArcSec",
    "acsc": "ArcCsc", "sinh": "Sinh", "cosh": "Cosh", "tanh": "Tanh", "coth": "Coth",
    "sech": "Sech", "csch": "Csch", "asinh": "ArcSinh", "acosh": "ArcCosh", "atanh": "ArcTanh",
    "acoth": "ArcCoth", "asech": "ArcSech", "acsch": "ArcCsch",
    "erf": "Erf", "erfi": "Erfi", "erfc": "Erfc", "Ei": "ExpIntegralEi", "li": "LogIntegral",
    "Si": "SinIntegral", "Ci": "CosIntegral", "Shi": "SinhIntegral", "Chi": "CoshIntegral",
    "gamma": "Gamma", "Gamma": "Gamma", "polylog": "PolyLog", "LambertW": "ProductLog",
    "lambert_w": "ProductLog", "dilog": None, "exp_polar": "Exp",
}
CONSTS = {"pi": "Pi", "%pi": "Pi", "E": "E", "%e": "E", "I": "I", "%i": "I", "x": "x"}


def tokens(s):
    out, i = [], 0
    s = s.strip()
    while i < len(s):
        m = TOK.match(s, i)
        if not m or m.end() == i:
            raise Untranslatable(f"token at {s[i:i + 20]!r}")
        i = m.end()
        num, name, op = m.groups()
        out.append(("num", num) if num else ("name", name) if name else ("op", op))
    return out


def to_wl(system, text):
    """Infix answer text of `system` -> Mathilda InputForm (fully parenthesised)."""
    toks = tokens(text)
    pos = [0]
    bound = {}          # bound variable name -> Mathilda slot name (rootOf / Lambda)

    def peek():
        return toks[pos[0]] if pos[0] < len(toks) else ("end", "")

    def take(v=None):
        t = peek()
        if v is not None and t[1] != v:
            raise Untranslatable(f"expected {v!r}, got {t[1]!r}")
        pos[0] += 1
        return t

    def var(name):
        # %%H0, _t, _z ... : bound or generated symbols, renamed to plain Mathilda symbols
        return "xc" + re.sub(r"[^A-Za-z0-9]", "", name)

    def atom_name(n):
        if n in bound:
            return bound[n]
        if n in CONSTS:
            return CONSTS[n]
        if n.startswith("_") and system == "sympy":   # RootSum's polynomial variable
            return var(n)
        raise Untranslatable(f"symbol {n}")

    def args():
        take("(")
        a = []
        if peek()[1] != ")":
            a.append(expr(0))
            while peek()[1] == ",":
                take(",")
                a.append(expr(0))
        take(")")
        return a

    def tuple_list(t):
        # a parsed SymPy tuple arrives as ("tuple", [...]); a lone element is a 1-tuple
        return "{" + ", ".join(t[1]) + "}" if isinstance(t, tuple) else "{" + t + "}"

    def call(name):
        if name == "rootOf":
            # rootOf(p, v): bind v while reading p
            take("(")
            save = pos[0]
            depth = 0
            while True:      # find the bound variable (last argument) first
                t = toks[pos[0]]
                if t[1] in "([":
                    depth += 1
                elif t[1] in ")]":
                    if depth == 0:
                        break
                    depth -= 1
                elif t[1] == "," and depth == 0:
                    vtok = toks[pos[0] + 1]
                pos[0] += 1
            v = vtok[1]
            pos[0] = save
            if "#1" in bound.values():
                # a rootOf inside another's polynomial (FriCAS's second root of a
                # deflated polynomial): Root over algebraic coefficients does not
                # evaluate numerically in Mathilda, so it is not translated
                raise Untranslatable("nested rootOf")
            bound[v] = "#1"
            p = expr(0)
            if "Root[" in p:
                raise Untranslatable("nested rootOf")
            take(",")
            take()
            take(")")
            del bound[v]
            return f"Root[({p}) &, 1]"
        if name == "Lambda":
            take("(")
            v = take()[1]
            take(",")
            bound[v] = "#1"
            body = expr(0)
            del bound[v]
            take(")")
            return ("lambda", body)
        a = args()
        if name == "RootSum":
            if len(a) != 2 or not isinstance(a[1], tuple) or a[1][0] != "lambda":
                raise Untranslatable("RootSum form")
            poly = a[0]
            # SymPy prints the polynomial in its own variable (_z): make it the slot
            m = re.findall(r"\bxc[a-z]\b", poly)
            if len(set(m)) != 1:
                raise Untranslatable("RootSum polynomial variable")
            poly = re.sub(r"\b" + m[0] + r"\b", "#1", poly)
            return f"RootSum[({poly}) &, ({a[1][1]}) &]"
        if name == "hyper" and system == "sympy":
            return f"HypergeometricPFQ[{tuple_list(a[0])}, {tuple_list(a[1])}, {a[2]}]"
        if name == "atan2" or (name == "atan" and len(a) == 2):
            return f"ArcTan[{a[1]}, {a[0]}]"
        if any(isinstance(v, tuple) for v in a):
            raise Untranslatable(f"tuple argument to {name}")
        head = FUNCS.get(name)
        if head is None:
            raise Untranslatable(f"function {name}")
        return f"{head}[{', '.join(a)}]"

    def prefix():
        t = take()
        if t[0] == "num":
            v = t[1]
            if "." in v or "e" in v.lower() or "b" in v:
                raise Untranslatable("inexact number")
            return v
        if t[1] == "-":
            return f"(-({expr(25)}))"
        if t[1] == "+":
            return expr(25)
        if t[1] == "(":
            first = expr(0)
            if peek()[1] == ",":          # a SymPy tuple
                items = [first]
                while peek()[1] == ",":
                    take(",")
                    if peek()[1] == ")":
                        break
                    items.append(expr(0))
                take(")")
                return ("tuple", items)
            take(")")
            return first if isinstance(first, tuple) else f"({first})"
        if t[0] == "name":
            if peek()[1] == "(":
                return call(t[1])
            return atom_name(t[1])
        raise Untranslatable(f"unexpected {t[1]!r}")

    BIN = {"+": (10, "L"), "-": (10, "L"), "*": (20, "L"), "/": (20, "L"), "^": (30, "R"), "**": (30, "R")}

    def expr(minp):
        lhs = prefix()
        while True:
            op = peek()[1]
            if peek()[0] != "op" or op not in BIN:
                return lhs
            p, assoc = BIN[op]
            if p < minp:
                return lhs
            take()
            rhs = expr(p + 1 if assoc == "L" else p)
            if isinstance(lhs, tuple) or isinstance(rhs, tuple):
                raise Untranslatable("tuple in arithmetic")
            o = "^" if op == "**" else op
            lhs = f"({lhs}){o}({rhs})"

    out = expr(0)
    if peek()[0] != "end":
        raise Untranslatable(f"trailing {peek()[1]!r}")
    if isinstance(out, tuple):
        raise Untranslatable("tuple result")
    return out


def needs_check(rec):
    """An ANSWER whose own check did not pass at every point."""
    nums = rec.get("nums") or []
    return rec.get("status") == "ANSWER" and not (nums and all(v == "ok" for v in nums))


def ckey(base, text):
    return "xcheck|" + base + "|" + hashlib.sha1(text.encode()).hexdigest()[:12]


XC_TMPL = r"""
xcRun[id_String, f_, s_String] := Module[{F, r},
  F = Quiet[Check[ToExpression[s], $Failed]];
  r = If[F === $Failed, "undecided-parse",
    Replace[TimeConstrained[Quiet[PossibleZeroQ[D[F, x] - f]], TMO, $Aborted],
      {True -> "zero", False -> "nonzero", _ -> "undecided"}]];
  (* True is evidence only where the expression evaluates numerically *)
  If[r === "zero" && ! TrueQ[TimeConstrained[Quiet[And @@ Table[
        NumberQ[N[(D[F, x] - f) /. x -> p, 30]], {p, {37/19, 1/3 + 2 I/7}}]], TMO, False]],
    r = "undecided-nonnumeric"];
  Print["@@X\t", id, "\t", r]];
"""


def run_batch(items, tmpdir):
    """items: [(id, f FullForm, wl text)] -> {id: verdict}"""
    import tempfile
    def runner(chunk):
        body = "".join(f'xcRun["{i}", {f}, "{w}"];\n' for i, f, w in chunk)
        fd, path = tempfile.mkstemp(suffix=".m", dir=tmpdir)
        with os.fdopen(fd, "w") as fh:
            fh.write(XC_TMPL.replace("TMO", str(PZQ_TIMEOUT)) + body)
        out, err, to, wall = A.run_proc([A.MATHILDA, "-file", path], 20 + len(chunk) * (PZQ_TIMEOUT + 15))
        os.unlink(path)
        got = {}
        for ln in out.splitlines():
            p = ln.split("\t")
            if len(p) == 3 and p[0] == "@@X":
                got[p[1]] = p[2]
        return got
    res = runner(items)
    for it in items:                      # a crash or a runaway takes its batch: rerun alone
        if it[0] not in res:
            res.update(runner([it]))
            res.setdefault(it[0], "undecided")
    return res
