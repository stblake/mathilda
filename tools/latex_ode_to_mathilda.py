#!/usr/bin/env python3
"""
latex_ode_to_mathilda.py — convert a 12000.org "Solving ODEs" section table
(tex4ht HTML) into a Mathilda/Mathematica-syntax ODE corpus (.m).

The 12000.org notes (Nasser M. Abbasi) render each ODE as a MathJax LaTeX block
inside an 8-column table (# | ID | ODE | CAS classification | Maple | Mma | Sympy |
time). This tool extracts every row, converts the LaTeX ODE to Mathilda input
syntax (`y'[x]`, `Sqrt[..]`, `Tan[..]`, `E^..`), auto-detects the dependent
function(s) and the independent variable, and emits a `List` of records

    {"<label>-N", <equation(s)>, <function(s)>, <indvar>, "<CAS classif>", <sympy?>}

suitable for `DSolve[...]` and for the fork-per-case corpus harness
(`tests/test_dsolve_corpus.c`). Systems (multiple dependent functions) are emitted
with a `List` equation and a `List` of functions and are tagged
`system_of_ODEs` in the classification field so the scalar harness can skip them.

Usage:
    python3 tools/latex_ode_to_mathilda.py <section.html> <out.m> [--label 2.1.2]

The output is validated to PARSE (round-tripped through the Mathilda binary's
NDJSON pipe mode) by the companion check in the campaign; this tool itself only
performs the textual conversion so it has no build/runtime dependency.
"""
import re, html, sys, argparse

GREEK = {r'\nu':'nu', r'\lambda':'lam', r'\mu':'mu', r'\alpha':'alpha',
         r'\beta':'beta', r'\gamma':'gam', r'\delta':'delta', r'\omega':'om',
         r'\sigma':'sig', r'\theta':'theta', r'\vartheta':'theta', r'\phi':'phi',
         r'\varphi':'phi', r'\rho':'rho', r'\kappa':'kap', r'\epsilon':'eps',
         r'\varepsilon':'eps', r'\pi':'Pi', r'\tau':'tau', r'\xi':'xi',
         r'\zeta':'zeta', r'\eta':'eta', r'\chi':'chi', r'\psi':'psi'}
FUNCS = {r'\arcsin':'ArcSin', r'\arccos':'ArcCos', r'\arctan':'ArcTan',
         r'\sinh':'Sinh', r'\cosh':'Cosh', r'\tanh':'Tanh', r'\coth':'Coth',
         r'\sech':'Sech', r'\csch':'Csch', r'\sin':'Sin', r'\cos':'Cos',
         r'\tan':'Tan', r'\cot':'Cot', r'\sec':'Sec', r'\csc':'Csc',
         r'\ln':'Log', r'\log':'Log', r'\exp':'Exp', r'\sign':'Sign',
         r'\wp':'WeierstrassP'}
ARBFUN = set('fgh')                 # conventional arbitrary-function letters
INDVAR_PREF = ['x', 't', 'z', 's', 'r', 'u', 'v']


def strip_array(tex):
    tex = re.sub(r'\\begin\s*\{array\}', '', tex)
    tex = re.sub(r'\\end\s*\{array\}', '', tex)
    tex = re.sub(r'\{>\{\\displaystyle\s*\}[a-z]\s*@\{\\;\}\s*>\{\\displaystyle\s*\}[a-z]\s*\}', '', tex)
    return tex


def find_matching(s, i):
    op = s[i]; cl = {'(': ')', '{': '}', '[': ']'}[op]; depth = 0; j = i
    while j < len(s):
        if s[j] == op: depth += 1
        elif s[j] == cl:
            depth -= 1
            if depth == 0: return j
        j += 1
    return -1


def replace_frac(s):
    while True:
        m = re.search(r'\\frac\s*', s)
        if not m: break
        i = m.end()
        while i < len(s) and s[i] == ' ': i += 1
        if i >= len(s) or s[i] != '{': s = s[:m.start()] + ' ' + s[m.end():]; continue
        a0 = i; a1 = find_matching(s, a0); j = a1 + 1
        while j < len(s) and s[j] == ' ': j += 1
        if j >= len(s) or s[j] != '{': s = s[:m.start()] + ' ' + s[m.end():]; continue
        b0 = j; b1 = find_matching(s, b0); A = s[a0 + 1:a1]; B = s[b0 + 1:b1]
        s = s[:m.start()] + '((' + A + ')/(' + B + '))' + s[b1 + 1:]
    return s


def replace_sqrt(s):
    while True:
        m = re.search(r'\\sqrt\s*', s)
        if not m: break
        i = m.end()
        while i < len(s) and s[i] == ' ': i += 1
        if i < len(s) and s[i] == '[':
            n1 = find_matching(s, i); nexpr = s[i + 1:n1]; i = n1 + 1
            while i < len(s) and s[i] == ' ': i += 1
            if i < len(s) and s[i] == '{':
                b1 = find_matching(s, i); B = s[i + 1:b1]
                s = s[:m.start()] + '((' + B + ')^(1/(' + nexpr + ')))' + s[b1 + 1:]; continue
        if i < len(s) and s[i] == '{':
            b1 = find_matching(s, i); B = s[i + 1:b1]
            s = s[:m.start()] + 'Sqrt[' + B + ']' + s[b1 + 1:]; continue
        s = s[:m.start()] + ' ' + s[m.end():]
    return s


def fn_paren_to_bracket(s, name):
    """Turn `Name( ... )` into `Name[ ... ]` (immediately-following paren only)."""
    out = s; idx = 0
    while True:
        k = out.find(name, idx)
        if k < 0: break
        if k > 0 and out[k - 1].isalpha(): idx = k + len(name); continue
        j = k + len(name)
        while j < len(out) and out[j] == ' ': j += 1
        if j < len(out) and out[j] == '(':
            cl = find_matching(out, j)
            if cl > 0:
                out = out[:j] + '[' + out[j + 1:cl] + ']' + out[cl + 1:]
                idx = k + len(name); continue
        idx = k + len(name)
    return out


def normalize_subscripts(tex):
    """x_{1} / x_1 -> x1 (valid Mathematica symbol; used by subscripted systems)."""
    tex = re.sub(r'([A-Za-z])_\s*\{\s*([0-9]+)\s*\}', r'\1\2', tex)
    tex = re.sub(r'([A-Za-z])_\s*([0-9])', r'\1\2', tex)
    return tex


def detect_symbols(rows_tex):
    joined = ' '.join(rows_tex)
    primed = set(re.findall(r'([A-Za-z][0-9]*)\s*\^\s*\{\s*(?:\\prime\s*)+\}', joined))
    mains = sorted(c for c in primed if c not in ARBFUN)
    if not mains: mains = ['y']
    arbs = set(c for c in primed if c in ARBFUN)
    for c in ARBFUN:
        if re.search(r'(?<![A-Za-z])' + c + r'\s*\\left', joined) or \
           re.search(r'(?<![A-Za-z])' + c + r'\s*\(', joined):
            arbs.add(c)
    main_set = set(mains); indvar = None
    # The independent variable is chosen from the MAIN (non-condition) rows only, so
    # that parameters appearing solely in an initial condition (y(a)=b) never win, and
    # a swapped-variable ODE (x = x(y), y independent) is detected from its own body.
    main_rows = [r for r in rows_tex if not is_condition_row(r, mains)]
    mjoined = ' '.join(main_rows) if main_rows else joined
    for cand in INDVAR_PREF:                             # preferred standard letter, present
        if cand in main_set: continue
        if re.search(r'(?<![A-Za-z0-9])' + cand + r'(?![A-Za-z0-9])', mjoined):
            indvar = cand; break
    if indvar is None:                                   # any present non-main letter (x=x(y))
        present = re.findall(r'(?<![A-Za-z0-9\\])([A-Za-z])(?![A-Za-z0-9])', mjoined)
        for cand in sorted(set(present)):
            if cand in main_set or cand in ARBFUN or cand == 'e': continue
            indvar = cand; break
    if indvar is None:                                   # autonomous: fresh standard letter
        for cand in INDVAR_PREF:
            if cand not in main_set: indvar = cand; break
    if indvar is None: indvar = 't'
    return mains, sorted(arbs), indvar


def _apply_arbfun(s, f, mains, indvar, protected):
    out = ''; i = 0
    while i < len(s):
        if s[i] == f and (i == 0 or not s[i - 1].isalpha()) and (i + 1 >= len(s) or not s[i + 1].isalpha()):
            j = i + 1; primes = 0
            m = re.match(r'\s*\^\s*\{\s*((?:\\prime\s*)+)\}', s[j:])
            if m: primes = m.group(1).count(r'\prime'); j += m.end()
            k = j
            while k < len(s) and s[k] == ' ': k += 1
            if k < len(s) and s[k] == '(':
                cl = find_matching(s, k); arg = s[k + 1:cl]; j = cl + 1
                argc = convert_side(arg, mains, [], indvar)
                final = f + ("'" * primes) + '[' + argc + ']'
            else:
                final = f + ("'" * primes) + '[' + indvar + ']'
            protected.append(final); out += '\x07%d\x07' % (len(protected) - 1); i = j
        else:
            out += s[i]; i += 1
    return out


def convert_side(expr, mains, arbs, indvar):
    protected = []; s = expr; extra = set()
    for nm in re.findall(r'\\operatorname\s*\{([^{}]*)\}', s): extra.add(nm.strip())
    s = re.sub(r'\\operatorname\s*\{([^{}]*)\}', lambda m: m.group(1).strip(), s)
    for mm in re.findall(r'\\textit\s*\{\s*\\?_?(F\d+)\s*\}', s): extra.add('Maple' + mm)
    s = re.sub(r'\\textit\s*\{\s*\\?_?(F\d+)\s*\}', lambda m: 'Maple' + m.group(1), s)
    s = re.sub(r'\\textit\s*\{([^{}]*)\}', r'\1', s)
    s = s.replace(r'\left', '').replace(r'\right', '')
    s = re.sub(r'\\[,;!> ]', ' ', s)
    s = s.replace(r'\cdot', '*').replace(r'\times', '*')
    s = s.replace(r'{\mathrm e}', 'E').replace(r'\mathrm{e}', 'E').replace(r'{\rm e}', 'E')
    s = re.sub(r'\\mathrm\s*\{([^{}]*)\}', r'\1', s)
    s = replace_frac(s); s = replace_sqrt(s)
    for f in arbs: s = _apply_arbfun(s, f, mains, indvar, protected)
    for f in sorted(mains, key=len, reverse=True):
        s = re.sub(r'(?<![A-Za-z0-9\x00])' + re.escape(f) + r'\s*\^\s*\{\s*((?:\\prime\s*)+)\}',
                   lambda m, f=f: '\x00' + f + '\x00' + str(m.group(1).count(r'\prime')) + '\x00', s)
        s = re.sub(r'(?<![A-Za-z0-9\x00])' + re.escape(f) + r'(?![A-Za-z0-9])',
                   '\x00' + f + '\x000\x00', s)
    for k in sorted(FUNCS, key=len, reverse=True): s = s.replace(k, FUNCS[k] + ' ')
    for k in sorted(GREEK, key=len, reverse=True): s = re.sub(re.escape(k) + r'(?![A-Za-z])', GREEK[k], s)
    s = re.sub(r'\\([A-Za-z]+)', r'\1', s)
    s = s.replace('{', '(').replace('}', ')')
    s = re.sub('\x00([A-Za-z][0-9]*)\x00([0-9]+)\x00',
               lambda m: m.group(1) + ("'" * int(m.group(2))) + '[' + indvar + ']', s)
    for h in (set(FUNCS.values()) | {'Sqrt', 'Exp'} | extra): s = fn_paren_to_bracket(s, h)
    for idx, val in enumerate(protected): s = s.replace('\x07%d\x07' % idx, val)
    s = re.sub(r'\s*\[\s*', '[', s); s = re.sub(r'\s*\]', ']', s)
    s = re.sub(r'\(\s+', '(', s); s = re.sub(r'\s+\)', ')', s)
    s = re.sub(r'\s+', ' ', s).strip()
    return s


def is_condition_row(row, mains):
    """True iff the row is an initial/boundary condition `y(P)=V` / `y'(P)=V`.

    A condition row's LHS is *solely* the function application.  The earlier
    loose `<main>\\left(`-anywhere search misread multiplication (`y^2 (y' x+y)`,
    `x (5-x)`) as an application and silently dropped the whole ODE row
    (§2.2.2-109/119/129/166/175-178).  So anchor to the LHS (before the first
    `=`) and require that, after removing `\\left`/`\\right`, the entire LHS is
    exactly `y(...)` / `y'(...)` — the matching close-paren must end the LHS."""
    r = row.replace('&', '')
    if '=' not in r:
        return False
    lhs = r.split('=', 1)[0].replace(r'\left', '').replace(r'\right', '').strip()
    for f in mains:
        m = re.match(r'\s*' + re.escape(f) + r'\s*(?:\^\s*\{\s*(?:\\prime\s*)+\}\s*)?\(', lhs)
        if m:
            cl = find_matching(lhs, m.end() - 1)
            if cl >= 0 and lhs[cl + 1:].strip() == '':
                return True
    return False


def convert_condition(row, mains, indvar):
    """Convert an initial-condition row `y(P)=V` / `y'(P)=V` into a Mathilda
    point equation `y[P] == V` / `y'[P] == V`. The point P and value V are
    constants/parameters (numbers, `a`, `Pi/2`, `2 E`, ...) converted WITHOUT
    the y->y[iv] substitution (mains=[]), so `0`/`a`/`Pi/2` stay literal."""
    r = row.replace('&', '')
    if '=' not in r: return None
    lhs, rhs = r.split('=', 1)
    lhs = lhs.replace(r'\left', '').replace(r'\right', '')
    m = re.match(r'\s*([A-Za-z][0-9]*)\s*((?:\^\s*\{\s*(?:\\prime\s*)+\}\s*)?)\(', lhs)
    if not m or m.group(1) not in mains: return None
    sym = m.group(1); primes = m.group(2).count(r'\prime')
    popen = m.end() - 1; pclose = find_matching(lhs, popen)
    if pclose < 0: return None
    point = convert_side(lhs[popen + 1:pclose], [], [], indvar)
    value = convert_side(rhs, [], [], indvar)
    return '%s%s[%s] == %s' % (sym, "'" * primes, point, value)


def convert_row(tex):
    tex = normalize_subscripts(strip_array(tex))
    rows = [r for r in re.split(r'\\\\', tex) if r.strip()]
    mains, arbs, indvar = detect_symbols(rows)
    eqs = []; conds = []
    for r in rows:
        if '=' not in r.replace('&', ''): continue
        if is_condition_row(r, mains):
            c = convert_condition(r, mains, indvar)
            if c: conds.append(c)
        else:
            lhs, rhs = r.replace('&', '').split('=', 1)
            eq = convert_side(lhs, mains, arbs, indvar) + ' == ' + convert_side(rhs, mains, arbs, indvar)
            eqs.append(eq)
    return mains, arbs, indvar, eqs, conds


def parse_table(html_text):
    """Yield dicts {n, tex, classif, sympy} from the tex4ht problems table.

    tex4ht numbers each section's table differently (`TBL-4-...` for \u00a72.1.2,
    `TBL-12-...` for \u00a72.2.1) and the column order also varies between sections
    (\u00a72.1.2 carries an extra ID column, so ODE sits at col 3; \u00a72.2.1 has ODE at
    col 2). So we (1) auto-select the table with the most `\\[..\\]` ODE blocks,
    and (2) read the column indices for #, ODE, classification and Sympy from
    that table's header row rather than hard-coding them."""
    cell_re = re.compile(r"id='TBL-(\d+)-(\d+)-(\d+)'[^>]*>(.*?)</td>", re.S)
    tables = {}                                # tnum -> {row -> {col -> cell}}
    for m in cell_re.finditer(html_text):
        t, r, c = int(m.group(1)), int(m.group(2)), int(m.group(3))
        tables.setdefault(t, {}).setdefault(r, {})[c] = m.group(4)
    if not tables: return

    def clean(x):
        x = re.sub(r'<!--.*?-->', ' ', x, flags=re.S)
        x = re.sub(r'<[^>]+>', ' ', x)
        return re.sub(r'\s+', ' ', html.unescape(x)).strip()

    def rawtex(cell):
        mm = re.search(r'\\\[(.*?)\\\]', cell, re.S)
        return html.unescape(re.sub(r'<[^>]+>', ' ', mm.group(1))).strip() if mm else ''

    # the ODE table is the one with the most rows carrying a LaTeX ODE block
    def n_odes(tbl):
        return sum(1 for row in tbl.values() for cell in row.values() if '\\[' in cell)
    rows = tables[max(tables, key=lambda t: n_odes(tables[t]))]

    # detect columns from the header row (the one whose cells name #, ODE, ...)
    col = {'n': 1, 'ode': 2, 'classif': 3, 'sympy': 7}
    for r in sorted(rows):
        labels = {clean(v).lower(): c for c, v in rows[r].items()}
        if '#' in labels and 'ode' in labels:
            col['n'] = labels['#']; col['ode'] = labels['ode']
            for txt, key in (('classification', 'classif'), ('sympy', 'sympy')):
                hit = next((c for lab, c in labels.items() if txt in lab), None)
                if hit is not None: col[key] = hit
            break

    for r in sorted(rows):
        cells = rows[r]; idx = clean(cells.get(col['n'], ''))
        if idx == '#' or not idx or not re.search(r'\d', idx): continue
        yield dict(n=idx, tex=rawtex(cells.get(col['ode'], '')),
                   classif=clean(cells.get(col['classif'], '')),
                   sympy=('\u2713' in cells.get(col['sympy'], '')))


def mm_escape(s):
    return s.replace('\\', '\\\\').replace('"', '\\"')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('html'); ap.add_argument('out')
    ap.add_argument('--label', default='2.1.2')
    ap.add_argument('--url', default='https://12000.org/my_notes/solving_ODE/current_version/indexsubsection2.htm')
    args = ap.parse_args()

    text = open(args.html, encoding='utf-8', errors='replace').read()
    recs = list(parse_table(text))
    out_lines = []; n_scalar = 0; n_system = 0; n_ivp = 0
    for rec in recs:
        mains, arbs, indvar, eqs, conds = convert_row(rec['tex'])
        label = '%s-%s' % (args.label, re.sub(r'[^0-9]', '', rec['n']))
        classif = rec['classif']
        sy = 'True' if rec['sympy'] else 'False'
        if len(mains) > 1 or len(eqs) > 1:       # system (or multi-equation row)
            n_system += 1
            eqn = '{' + ', '.join(eqs + conds) + '}'
            fns = '{' + ', '.join(mains) + '}'
            cl = classif if classif else 'system_of_ODEs'
            out_lines.append('  {"%s", %s, %s, %s, "%s", %s}' %
                             (label, eqn, fns, indvar, mm_escape(cl), sy))
        else:                                    # scalar
            n_scalar += 1
            ode = eqs[0] if eqs else '(* NO ODE ROW: %s *) True == True' % mm_escape(rec['tex'][:40])
            # IVP: emit the DSolve-native list {ode, ic1, ...}; else the bare equation.
            eqn = '{' + ', '.join([ode] + conds) + '}' if conds else ode
            if conds: n_ivp += 1
            out_lines.append('  {"%s", %s, %s, %s, "%s", %s}' %
                             (label, eqn, mains[0], indvar, mm_escape(classif), sy))

    header = (
        '(* DE_examples_%s.m --- ODE corpus from 12000.org "Solving ODEs" section %s\n'
        '   Source: %s (Table, %d rows).\n'
        '   Generated by tools/latex_ode_to_mathilda.py --- DO NOT hand-edit; regenerate.\n\n'
        '   Record: {"label", equation(s), function(s), indVar, "MapleClassification", sympySolved}.\n'
        '   Scalar ODE: equation is one ODE, function is a symbol.\n'
        '   Scalar IVP (%d): equation is the DSolve-native list {ode, ic1, ...} (each ic a\n'
        '     point equation y[x0]==v / y\'[x0]==v); function stays a symbol, so the scalar\n'
        '     harness still handles it. The harness verifies the ODE residual AND every ic.\n'
        '   Systems (%d): equation is a List, function is a List; classification\n'
        '   contains "system_of_ODEs" so the scalar harness skips them.\n'
        '   Consumed (parsed, NOT evaluated) by tests/test_dsolve_corpus.c. *)\n\n'
        '{\n' % (args.label, args.label, args.url, len(recs), n_ivp, n_system))
    with open(args.out, 'w') as f:
        f.write(header + ',\n'.join(out_lines) + '\n}\n')
    sys.stderr.write('wrote %s: %d records (%d scalar [%d IVP], %d systems)\n' %
                     (args.out, len(recs), n_scalar, n_ivp, n_system))


if __name__ == '__main__':
    main()
