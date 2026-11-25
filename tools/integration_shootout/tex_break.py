"""tex_break.py -- break a long TeX formula into multline* lines.

    display(lhs, tex)   ->  TeX for "lhs = tex" (integral(tex): for "\\int tex dx"): one inline line when
                            it fits the page, else a multline* environment

A formula is broken only at a binary + or - that sits at the top level: outside
every {...} group and every \\left ... \\right pair, since TeX cannot end a line
inside either.  The terms between such breaks are packed greedily into lines of
at most WIDTH units (one unit is roughly the width of a math letter at 10pt; the
page holds about 95).  A single term wider than a line is opened up when its
widest \\left X ... \\right Y pair is top-level in the term: the pair becomes
\\Bigl X ... \\Bigr Y (\\bigl/\\bigr when no fraction is inside), which has the
same look without the pairing, so the + and - inside it become break points.
A line that is still too wide (one huge fraction, a radical over a long
radicand) is scaled to the line width with \\resizebox: smaller, but on the page.

The widths are estimates computed from the token stream (a fraction is as wide
as the wider of its numerator and denominator, a script is 0.7 of its base, ...),
which is all the packing needs.
"""
import re

WIDTH = 72           # units per multline line (a unit is ~6pt of displaystyle math; the text block is ~82)
INLINE_WIDTH = 78    # an unbroken formula up to this wide stays on one inline line

TOKEN = re.compile(r"\\[A-Za-z]+\*?|\\.|[{}]|\s+|.", re.S)
FUNC_WORDS = {"log", "exp", "sin", "cos", "tan", "cot", "sec", "csc", "sinh", "cosh", "tanh", "coth",
              "sech", "csch", "arcsin", "arccos", "arctan", "ln", "erf", "operatorname", "text", "mathrm"}


def tokens(tex):
    return [t for t in TOKEN.findall(tex) if not t.isspace()]


def _group(toks, i):
    """toks[i] starts an argument: a {...} group or a single token -> (inner tokens, next index)"""
    if i >= len(toks):
        return [], i
    if toks[i] != "{":
        return [toks[i]], i + 1
    d, j = 0, i
    while j < len(toks):
        if toks[j] == "{":
            d += 1
        elif toks[j] == "}":
            d -= 1
            if d == 0:
                return toks[i + 1:j], j + 1
        j += 1
    return toks[i + 1:], len(toks)


def width(toks):
    """estimated width of a token list, in math-letter units"""
    w, i, n = 0.0, 0, len(toks)
    while i < n:
        t = toks[i]
        if t in ("\\frac", "\\dfrac", "\\tfrac"):
            a, i = _group(toks, i + 1)
            b, i = _group(toks, i)
            w += max(width(a), width(b)) + 0.6
            continue
        if t in ("^", "_"):
            a, i = _group(toks, i + 1)
            w += 0.7 * width(a)
            continue
        if t == "\\sqrt":
            i += 1
            if i < n and toks[i] == "[":
                j = toks.index("]", i) if "]" in toks[i:] else i
                w += 0.5 * width(toks[i + 1:j])
                i = j + 1
            a, i = _group(toks, i)
            w += width(a) + 1.2
            continue
        if t in ("\\left", "\\right") or re.fullmatch(r"\\[Bb]ig+[lrm]?", t):
            w += 0.6
            i += 2
            continue
        if t in ("{", "}"):
            i += 1
            continue
        if t in ("+", "-", "=", "<", ">"):
            w += 1.8
        elif t in ("\\ ", "\\,", "\\;", "\\:", "\\!", "~"):
            w += 0.3
        elif t.startswith("\\") and t[1:].rstrip("*") in FUNC_WORDS:
            w += 0.8 * len(t) - 0.8
        elif t.startswith("\\"):
            w += 1.0
        else:
            w += 1.0
        i += 1
    return w


SIZED = re.compile(r"\\[Bb]igg?[lrm]?$")


def _scan(toks):
    """yield (index, token, top) where top says the token is outside every group, every
    \\left pair and every plain ( ) pair (TeXForm writes a function's argument in plain
    parentheses; a break inside one is a last resort, made by open_widest_pair)"""
    brace = left = paren = 0
    i = 0
    while i < len(toks):
        t = toks[i]
        yield i, t, brace == 0 and left == 0 and paren == 0
        if t == "{":
            brace += 1
        elif t == "}":
            brace -= 1
        elif t == "(":
            paren += 1
        elif t == ")":
            paren -= 1
        elif t == "\\left":
            left += 1
            i += 1          # the delimiter
        elif t == "\\right":
            left -= 1
            i += 1
        elif SIZED.match(t):
            i += 1          # \bigl( ... \bigr): sized and unpaired, so its delimiter is not nesting
        i += 1


def split_terms(toks):
    """split at top-level binary + / -; the sign stays with the term it begins"""
    cuts, prev = [], None
    for i, t, top in _scan(toks):
        if top and t in ("+", "-") and i > 0 and prev not in ("+", "-", "=", "(", "[", "{", "\\left", ","):
            cuts.append(i)
        prev = t
    out, a = [], 0
    for c in cuts:
        out.append(toks[a:c])
        a = c
    out.append(toks[a:])
    return [p for p in out if p]


def open_widest_pair(toks):
    """replace the widest top-level pair of toks -- \\left X ... \\right Y, or a plain ( ... ) --
    by \\Bigl X ... \\Bigr Y, so the + and - inside it become break points; None when there
    is no such pair"""
    best, stack, brace = None, [], 0
    i = 0
    while i < len(toks):
        t = toks[i]
        if t == "{":
            brace += 1
        elif t == "}":
            brace -= 1
        elif t in ("\\left", "("):
            stack.append((i, brace, t))
            if t == "\\left":
                i += 1
        elif t in ("\\right", ")"):
            if stack:
                j, b, kind = stack.pop()
                if not stack and b == 0 and brace == 0:
                    wd = width(toks[j:i + (2 if t == "\\right" else 1)])
                    if best is None or wd > best[0]:
                        best = (wd, j, i, kind)
            if t == "\\right":
                i += 1
        elif SIZED.match(t):
            i += 1
        i += 1
    if best is None:
        return None
    _, j, k, kind = best
    if kind == "(":
        inner, lo, hi, rest = toks[j + 1:k], "(", ")", toks[k + 1:]
    else:
        inner, lo, hi, rest = toks[j + 2:k], toks[j + 1], toks[k + 1], toks[k + 2:]
    size = "Big" if any(t in ("\\frac", "\\dfrac") for t in inner) else "big"
    # \Bigl. is a null delimiter, as \left. is
    return toks[:j] + [f"\\{size}l", lo] + inner + [f"\\{size}r", hi] + rest


def slash_widest_frac(toks):
    """replace the widest top-level \\frac{A}{B} of toks by \\Bigl(A\\Bigr)/\\Bigl(B\\Bigr): the
    same expression, but no longer one unbreakable box (TeX cannot end a line inside
    \\frac, and a fraction with a 5000-character numerator is wider than TeX's largest
    dimension, so not even \\resizebox can set it); None when there is none"""
    best = None
    for i, t, top in _scan(toks):
        if top and t in ("\\frac", "\\dfrac", "\\tfrac"):
            a, j = _group(toks, i + 1)
            b, k = _group(toks, j)
            wd = max(width(a), width(b))
            if best is None or wd > best[0]:
                best = (wd, i, a, b, k)
    if best is None:
        return None
    _, i, a, b, k = best
    return (toks[:i] + ["\\Bigl", "("] + a + ["\\Bigr", ")", "\\Big", "/", "\\Bigl", "("] + b +
            ["\\Bigr", ")"] + toks[k:])


def _pieces(toks, budget):
    """break toks into pieces no wider than budget where the structure allows"""
    out = []
    for term in split_terms(toks):
        while width(term) > budget:
            opened = open_widest_pair(term)
            if opened is None and width(term) > 4 * budget:
                opened = slash_widest_frac(term)    # last resort: a fraction far wider than a line
            if opened is None:
                break
            sub = split_terms(opened)
            if len(sub) > 1:
                out.extend(_pieces(opened, budget))
                term = None
                break
            term = opened
        if term:
            out.append(term)
    return out


def render(toks):
    s = ""
    for t in toks:
        if s and re.match(r"[A-Za-z]", t) and re.search(r"\\[A-Za-z]+$", s):
            s += " "
        s += t
    return s


def integral(tex):
    """TeX for the integral of tex dx: one inline line, or multline* lines; a sum is
    bracketed (\\Bigl( ... \\Bigr), unpaired, so a line may end inside it)"""
    toks = tokens(tex)
    if len(split_terms(toks)) > 1:
        if width(toks) + 6 <= INLINE_WIDTH:
            return _lines(["\\int", "\\left", "("], toks, ["\\right", ")", "\\,", "dx"])
        return _lines(["\\int", "\\Bigl", "("], toks, ["\\Bigr", ")", "\\,", "dx"])
    return _lines(["\\int"], toks, ["\\,", "dx"])


def display(lhs, tex):
    """TeX for lhs = tex, broken into multline* lines when it does not fit"""
    return _lines(tokens(lhs) + ["="], tokens(tex), [])


def _lines(head, toks, tail):
    """head toks tail on one inline line when it fits, else multline* lines broken in toks"""
    if width(head) + width(toks) + width(tail) + 2 <= INLINE_WIDTH:
        return r"\par\hspace*{1em}$\displaystyle " + render(head + toks + tail) + "$"
    lines, cur = [], list(head)
    for p in _pieces(toks, WIDTH):
        if len(cur) > len(head) and width(cur) + width(p) > WIDTH:
            lines.append(cur)
            cur = []
        cur = cur + p
    lines.append(cur + tail)
    body = []
    for ln in lines:
        s = render(ln)
        if width(ln) > WIDTH * 1.08:
            s = r"\resizebox{0.92\linewidth}{!}{$\displaystyle " + s + "$}"
        body.append(s)
    return "\\begin{multline*}\n" + " \\\\\n".join(body) + "\n\\end{multline*}"
