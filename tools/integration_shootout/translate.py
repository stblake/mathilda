"""translate.py -- FullForm parser and per-system emitters for the integration
shootout.

The corpus carries every integrand as Mathilda FullForm (`Plus[Times[x, Power[
Sec[x], 2]], Tan[x]]`): a grammar of heads, symbols, integers and strings that
parses without ambiguity.  `parse` turns it into nested tuples (head, args...)
-- a str for a symbol, an int for an integer -- and `emit(tree, system)` prints
the integrand in the syntax of sympy, Maxima, FriCAS or the Wolfram Language.

Every head the corpus can contain is in HEADS; an unknown head raises
Untranslatable, so a case is never silently mistranslated (the driver marks
that run ERROR "untranslatable" instead).
"""
import re

TOKEN = re.compile(r"\s*(?:(-?\d+)|([A-Za-z$][A-Za-z0-9$`]*)|(\[)|(\])|(,))")


class Untranslatable(Exception):
    pass


def parse(s):
    pos = 0
    toks = []
    while pos < len(s):
        m = TOKEN.match(s, pos)
        if not m or m.end() == pos:
            if s[pos:].strip() == "":
                break
            raise Untranslatable(f"bad FullForm near {s[pos:pos + 20]!r}")
        pos = m.end()
        if m.group(1) is not None:
            toks.append(int(m.group(1)))
        elif m.group(2) is not None:
            toks.append(("sym", m.group(2)))
        elif m.group(3):
            toks.append("[")
        elif m.group(4):
            toks.append("]")
        else:
            toks.append(",")
    i = 0

    def expr():
        nonlocal i
        t = toks[i]
        i += 1
        if isinstance(t, int):
            return t
        if not (isinstance(t, tuple) and t[0] == "sym"):
            raise Untranslatable(f"unexpected token {t!r}")
        name = t[1]
        if i < len(toks) and toks[i] == "[":
            i += 1
            args = []
            if toks[i] == "]":
                i += 1
                return (name,)
            while True:
                args.append(expr())
                if toks[i] == ",":
                    i += 1
                    continue
                if toks[i] == "]":
                    i += 1
                    break
                raise Untranslatable("expected , or ]")
            return (name, *args)
        return name

    tree = expr()
    if i != len(toks):
        raise Untranslatable("trailing tokens")
    return tree


# head -> {system: name}; None = same spelling as the Wolfram one lowercased
FUNCS = {
    "Log": ("log", "log", "log", "Log"),
    "Sin": ("sin", "sin", "sin", "Sin"), "Cos": ("cos", "cos", "cos", "Cos"),
    "Tan": ("tan", "tan", "tan", "Tan"), "Cot": ("cot", "cot", "cot", "Cot"),
    "Sec": ("sec", "sec", "sec", "Sec"), "Csc": ("csc", "csc", "csc", "Csc"),
    "ArcTan": ("atan", "atan", "atan", "ArcTan"), "ArcSin": ("asin", "asin", "asin", "ArcSin"),
    "ArcCos": ("acos", "acos", "acos", "ArcCos"),
    "Sinh": ("sinh", "sinh", "sinh", "Sinh"), "Cosh": ("cosh", "cosh", "cosh", "Cosh"),
    "Tanh": ("tanh", "tanh", "tanh", "Tanh"), "ArcSinh": ("asinh", "asinh", "asinh", "ArcSinh"),
    "ArcTanh": ("atanh", "atanh", "atanh", "ArcTanh"),
}
SYSTEMS = ("sympy", "maxima", "fricas", "wolfram")
SYMS = {  # symbol -> per-system spelling
    "x": ("x", "x", "x", "x"),
    "E": ("E", "%e", "%e", "E"),
    "Pi": ("pi", "%pi", "%pi", "Pi"),
    "I": ("I", "%i", "%i", "I"),
}
HEADS = set(FUNCS) | {"Plus", "Times", "Power", "Rational", "Complex"}


def emit(tree, system):
    k = SYSTEMS.index(system)

    def go(t):
        if isinstance(t, int):
            return f"({t})" if t < 0 else str(t)
        if isinstance(t, str):
            if t not in SYMS:
                raise Untranslatable(f"symbol {t}")
            return SYMS[t][k]
        h, args = t[0], t[1:]
        if h not in HEADS:
            raise Untranslatable(f"head {h}")
        if h == "Plus":
            return "(" + " + ".join(go(a) for a in args) + ")"
        if h == "Times":
            return "(" + "*".join(go(a) for a in args) + ")"
        if h == "Rational":
            p, q = args
            if system == "sympy":
                return f"Rational({p}, {q})"
            return f"({p}/{q})"
        if h == "Complex":
            a, b = (go(v) for v in args)
            return f"({a} + {b}*{SYMS['I'][k]})"
        if h == "Power":
            b, e = go(args[0]), go(args[1])
            op = "**" if system == "sympy" else "^"
            return f"({b}{op}{e})"
        name = FUNCS[h][k]
        if system == "wolfram":
            return name + "[" + ", ".join(go(a) for a in args) + "]"
        return name + "(" + ", ".join(go(a) for a in args) + ")"

    return go(tree)


def heads(tree, acc=None):
    """Every head (and symbol) in tree -- for scope checks."""
    acc = set() if acc is None else acc
    if isinstance(tree, tuple):
        acc.add(tree[0])
        for a in tree[1:]:
            heads(a, acc)
    elif isinstance(tree, str):
        acc.add(tree)
    return acc


# Non-elementary function names in each system's ANSWER text: an answer that
# uses one is SPECIAL (a closed form, but not an elementary one).
SPECIAL_RE = {
    "mathilda": re.compile(r"\b(Erf|Erfi|Erfc|ExpIntegralEi|ExpIntegralE|LogIntegral|SinIntegral|"
                           r"CosIntegral|SinhIntegral|CoshIntegral|PolyLog|Elliptic\w+|Hypergeometric\w+|"
                           r"AppellF1|Gamma|Fresnel\w+|ProductLog|MeijerG|Beta|InverseJacobi\w+|Jacobi\w+|"
                           r"Weierstrass\w+|Bessel\w+|Zeta|Piecewise|DifferentialRoot)\["),
    "sympy": re.compile(r"\b(erf|erfi|erfc|Ei|li|Li|Si|Ci|Shi|Chi|expint|polylog|elliptic_\w+|hyper|"
                        r"meijerg|gamma|uppergamma|lowergamma|fresnel\w|LambertW|besse\w+|zeta|Piecewise|"
                        r"appellf1)\("),
    "maxima": re.compile(r"\b(erf|erfi|erfc|expintegral_\w+|gamma_incomplete\w*|gamma|elliptic_\w+|"
                         r"hypergeometric|fresnel_\w|lambert_w|bessel_\w+|li\[|%f|beta)\b"),
    "fricas": re.compile(r"\b(erf|erfi|Ei|li|Si|Ci|Shi|Chi|dilog|polylog|ellipticF|ellipticE|ellipticPi|"
                         r"ellipticK|fresnelS|fresnelC|Gamma|lambertW|besselJ|besselY|besselI|besselK|"
                         r"hypergeometricF|meijerG|weierstrass\w*)\b"),
}
SPECIAL_RE["wolfram"] = SPECIAL_RE["mathilda"]

UNEVAL_RE = {
    "mathilda": re.compile(r"\bIntegrate\["),
    "wolfram": re.compile(r"\bIntegrate\["),
    "sympy": re.compile(r"\bIntegral\("),
    "maxima": re.compile(r"\bintegrate\("),
    "fricas": re.compile(r"\bintegral\(|%A|%B"),
}


def special_q(system, text):
    return bool(SPECIAL_RE[system].search(text or ""))


def uneval_q(system, text):
    return bool(UNEVAL_RE[system].search(text or ""))
