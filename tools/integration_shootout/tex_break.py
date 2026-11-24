"""tex_break.py -- break a long TeX formula into multline* lines.

    display(lhs, tex)   ->  TeX for "lhs = tex": one inline \\displaystyle line when
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


def _scan(toks):
    """yield (index, token, top) where top says the token is outside every group and \\left pair"""
    brace = left = 0
    i = 0
    while i < len(toks):
        t = toks[i]
        yield i, t, brace == 0 and left == 0
        if t == "{":
            brace += 1
        elif t == "}":
            brace -= 1
        elif t == "\\left":
            left += 1
            i += 1          # the delimiter
        elif t == "\\right":
            left -= 1
            i += 1
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
    """replace the widest top-level \\left X ... \\right Y of toks by \\Bigl X ... \\Bigr Y;
    None when there is no such pair"""
    best, stack, brace = None, [], 0
    i = 0
    while i < len(toks):
        t = toks[i]
        if t == "{":
            brace += 1
        elif t == "}":
            brace -= 1
        elif t == "\\left":
            stack.append((i, brace))
            i += 1
        elif t == "\\right":
            if stack:
                j, b = stack.pop()
                if not stack and b == 0 and brace == 0:
                    wd = width(toks[j:i + 2])
                    if best is None or wd > best[0]:
                        best = (wd, j, i)
            i += 1
        i += 1
    if best is None:
        return None
    _, j, k = best
    inner = toks[j + 2:k]
    size = "Big" if any(t in ("\\frac", "\\dfrac") for t in inner) else "big"
    # \Bigl. is a null delimiter, as \left. is
    return toks[:j] + [f"\\{size}l", toks[j + 1]] + inner + [f"\\{size}r", toks[k + 1]] + toks[k + 2:]


def _pieces(toks, budget):
    """break toks into pieces no wider than budget where the structure allows"""
    out = []
    for term in split_terms(toks):
        while width(term) > budget:
            opened = open_widest_pair(term)
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


def display(lhs, tex):
    """TeX for lhs = tex, broken into multline* lines when it does not fit"""
    toks = tokens(tex)
    if width(toks) + width(tokens(lhs)) + 2 <= INLINE_WIDTH:
        return r"\par\hspace*{1em}$\displaystyle " + lhs + " = " + tex + "$"
    head = tokens(lhs) + ["="]
    lines, cur = [], list(head)
    for p in _pieces(toks, WIDTH):
        if len(cur) > len(head) and width(cur) + width(p) > WIDTH:
            lines.append(cur)
            cur = []
        cur = cur + p
    lines.append(cur)
    body = []
    for ln in lines:
        s = render(ln)
        if width(ln) > WIDTH * 1.08:
            s = r"\resizebox{0.92\linewidth}{!}{$\displaystyle " + s + "$}"
        body.append(s)
    return "\\begin{multline*}\n" + " \\\\\n".join(body) + "\n\\end{multline*}"
