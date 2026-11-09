#!/usr/bin/env python3
"""
triage.py -- classify every G&R corpus record into Can-do / Algorithmic family /
Table, from a Mathilda Integrate run.

Usage:
    python3 -I tools/gr/triage.py --runs run1.jsonl[,run2.jsonl...] \
        --corpus GR2IntegralsCorpus.m[,GR3...] --out triage.json [--md triage.md]

<run*.jsonl> are tests/test_mathar_corpus.c outputs (MATHAR_CORPUS_JSON), one
line per (record, mode).  The corpus files supply integrand, limits, kind and
the book's value is read from the record's symbolic Integrate (Hold) is not
needed: the reference value text is taken from the verified JSON when given
(--verified py.json[,...]), which carries the value.

Buckets (per record):
  CAN-DO        S mode CORRECT: Integrate already returns a correct closed form
                with symbolic parameters.
  INSTANCE-ONLY N mode CORRECT, S not: the method exists, only the parametric
                generalisation (or Assumptions handling) is missing.
  BUG           WRONG / NONNUMERIC / CRASH in any mode (fix before features).
  <family>      algorithmic family (see FAMILIES): a general procedure closes
                the class; ranked by count in the report.
  TABLE         no general procedure applies with reasonable effort: isolated
                closed forms (special constants, values in functions outside
                the Meijer-G class, symbolic-order families whose answer is a
                finite sum with no hypergeometric closed form).

The family test is structural (heads of the integrand and value, interval,
kind).  It is a first cut: docs/design/gr_integrate_triage.md reviews each
family by sampling its members.
"""
import argparse
import json
import re
from collections import Counter, defaultdict

ELEM_HEADS = {'Sin', 'Cos', 'Tan', 'Cot', 'Sec', 'Csc', 'Sinh', 'Cosh', 'Tanh', 'Coth', 'Sech',
              'Csch', 'ArcSin', 'ArcCos', 'ArcTan', 'ArcCot', 'ArcSec', 'ArcCsc', 'ArcSinh',
              'ArcCosh', 'ArcTanh', 'ArcCoth', 'ArcSech', 'ArcCsch', 'Log', 'Exp', 'Sqrt',
              'Abs', 'Sign'}
STRUCT_HEADS = {'Integrate', 'Hold', 'List', 'Plus', 'Times', 'Power', 'Rational', 'Complex'}
ELLIPTIC = {'EllipticK', 'EllipticE', 'EllipticF', 'EllipticPi', 'JacobiSN', 'JacobiCN',
            'JacobiDN', 'JacobiAmplitude', 'MatharK', 'MatharE', 'MatharK1', 'MatharE1',
            'MatharF2', 'MatharE2', 'MatharPi'}
LIOUVILLIAN_SPECIAL = {'Erf', 'Erfc', 'Erfi', 'ExpIntegralEi', 'ExpIntegralE', 'SinIntegral',
                       'CosIntegral', 'SinhIntegral', 'CoshIntegral', 'LogIntegral', 'PolyLog',
                       'FresnelS', 'FresnelC', 'Gamma'}
MEIJER_CLASS = {'BesselJ', 'BesselY', 'BesselI', 'BesselK', 'HankelH1', 'HankelH2', 'StruveH',
                'StruveL', 'Hypergeometric0F1', 'Hypergeometric1F1', 'Hypergeometric2F1',
                'HypergeometricU', 'HypergeometricPFQ', 'MeijerG', 'WhittakerM', 'WhittakerW',
                'ParabolicCylinderD', 'LegendreP', 'LegendreQ', 'LaguerreL', 'HermiteH',
                'ChebyshevT', 'ChebyshevU', 'GegenbauerC', 'JacobiP', 'AiryAi', 'AiryBi',
                'Erf', 'Erfc', 'Gamma', 'ExpIntegralEi', 'ExpIntegralE', 'SinIntegral',
                'CosIntegral', 'KelvinBer', 'KelvinBei', 'KelvinKer', 'KelvinKei', 'AngerJ',
                'WeberE', 'SphericalBesselJ', 'SphericalBesselY', 'PolyLog', 'LerchPhi'}
NON_MEIJER = {'LommelS1', 'LommelS2', 'LommelU', 'LommelV', 'ThetaN', 'MathieuCe', 'MathieuSe',
              'AppellF1', 'AppellF2', 'AppellF3', 'AppellF4', 'LobachevskyL', 'VolterraNu',
              'VolterraMu', 'VolterraLambda', 'WeierstrassP', 'WeierstrassZeta',
              'WeierstrassSigma', 'BarnesG', 'MacRobertE', 'HumbertPhi1', 'NeumannO',
              'SchlaefliS', 'BatemanK'}
NUM_CONST = {'Pi', 'E', 'EulerGamma', 'Catalan', 'Zeta', 'Glaisher', 'Khinchin'}
GREEKISH = {'n', 'm', 'k', 'l', 'j', 'p', 'q', 'mu', 'nu', 'lam', 'alpha', 'beta', 'rho', 'sigma'}

FAMILIES = [
    # (name, description) -- order of the classifier below
    ('risch-elementary', 'indefinite, elementary integrand and elementary antiderivative (Risch, incl. algebraic)'),
    ('liouvillian-special', 'indefinite, elementary integrand, antiderivative needs Erf/Ei/Si/Li/PolyLog (Risch + special-function extension)'),
    ('elliptic-reduction', 'integrand algebraic of genus 1 / trig sqrt; value in F, E, Pi (Legendre/Carlson reduction)'),
    ('symbolic-order', 'parameter in an exponent/order (x^n, sin^n, Sum in value): closed form for symbolic n'),
    ('holonomic-indefinite', 'indefinite integral of special functions (Bessel x power, orthogonal polynomials)'),
    ('meijerg-mellin', 'definite over (0,inf) or (0,1) of Meijer-G-representable factors: Mellin convolution / Slater'),
    ('residue-rational-trig', 'definite rational x trig/exp(I) over (0,inf)/(-inf,inf), or trig-rational over a period: residues'),
    ('log-polylog', 'definite with Log/ArcTan/PolyLog of rational functions: differentiation under the integral, polylog ladders'),
    ('beta-gamma', 'definite of powers x exp/algebraic on (0,1)/(0,inf): Beta/Gamma (Euler, Mellin of elementary)'),
    ('trig-power-period', 'definite of trig/hyperbolic powers over (0, pi/2), (0, pi), (0, inf): Wallis/Beta, Fourier'),
    ('other-definite-elementary', 'definite elementary integrand not in the families above'),
]


def heads(s):
    return set(re.findall(r'([A-Z][A-Za-z0-9]*)\[', s or ''))


def interval(lim):
    l = (lim or '').replace(' ', '')
    if not l:
        return 'indef'
    if l in ('{0,Infinity}',):
        return '(0,inf)'
    if l in ('{-Infinity,Infinity}',):
        return '(-inf,inf)'
    if l in ('{0,1}',):
        return '(0,1)'
    if re.match(r'\{0,Pi(/\d+)?\}$', l) or l in ('{0,2*Pi}', '{-Pi,Pi}'):
        return 'period'
    if l.endswith(',Infinity}'):
        return '(a,inf)'
    return 'finite'


def symbolic_order(f, v):
    """x^n, Sin[x]^n, Sum[..] -- a parameter in an exponent or order."""
    if 'Sum[' in v or 'Product[' in v:
        return True
    for m in re.finditer(r'\^\(?([A-Za-z][A-Za-z0-9]*)', f):
        if m.group(1) in GREEKISH:
            return True
    return False


def classify(kind, f, lim, v):
    hf, hv = heads(f) - STRUCT_HEADS, heads(v) - STRUCT_HEADS
    if hv & NON_MEIJER or hf & NON_MEIJER:
        return 'TABLE', 'value/integrand outside the Meijer-G class (%s)' % ', '.join(sorted((hv | hf) & NON_MEIJER))
    if hv & ELLIPTIC or hf & ELLIPTIC:
        return 'elliptic-reduction', ''
    elementary_f = hf <= ELEM_HEADS
    if kind == 'indefinite':
        if symbolic_order(f, v):
            return 'symbolic-order', ''
        if not elementary_f:
            return 'holonomic-indefinite', ''
        if hv <= ELEM_HEADS:
            return 'risch-elementary', ''
        if hv & LIOUVILLIAN_SPECIAL:
            return 'liouvillian-special', ''
        return 'TABLE', 'indefinite, value in %s' % ', '.join(sorted(hv))
    iv = interval(lim)
    if hf & MEIJER_CLASS or (hf - ELEM_HEADS):
        if iv in ('(0,inf)', '(0,1)', '(a,inf)', 'finite', '(-inf,inf)'):
            return 'meijerg-mellin', ''
        return 'TABLE', 'special-function integrand over %s' % iv
    if symbolic_order(f, v) and iv in ('period', '(0,1)', 'finite'):
        return 'trig-power-period' if hf & {'Sin', 'Cos', 'Tan', 'Cot', 'Sec', 'Csc'} else 'beta-gamma', ''
    if hf & {'Log', 'ArcTan', 'ArcSin', 'ArcCos', 'ArcCot', 'ArcTanh', 'ArcSinh', 'ArcCosh'}:
        return 'log-polylog', ''
    trig = hf & {'Sin', 'Cos', 'Tan', 'Cot', 'Sec', 'Csc', 'Sinh', 'Cosh', 'Tanh', 'Coth', 'Sech', 'Csch'}
    if iv == 'period' and trig:
        return 'trig-power-period', ''
    if iv in ('(0,inf)', '(-inf,inf)') and (trig or 'Exp' in hf or 'E^' in f):
        return 'residue-rational-trig', ''
    if iv in ('(0,1)', '(0,inf)', 'finite', '(a,inf)'):
        return 'beta-gamma', ''
    return 'other-definite-elementary', ''


REC_RE = re.compile(r'\{"(GR[^"]+)", "(definite|indefinite)",\n  Hold\[Integrate\[(.*?)\]\],\n', re.S)


def load_corpus(path):
    txt = open(path).read()
    out = {}
    for m in REC_RE.finditer(txt):
        rid, kind, body = m.groups()
        body = re.sub(r', Assumptions -> .*$', '', body, flags=re.S)
        if kind == 'definite':
            mm = re.match(r'(.*), \{(\w+), (.*)\}$', body, re.S)
            f, lim = (mm.group(1), '{%s}' % mm.group(3)) if mm else (body, '')
        else:
            f, lim = body.rsplit(', ', 1)[0], ''
        out[rid] = {'kind': kind, 'f': f, 'lim': lim}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--runs', required=True)
    ap.add_argument('--corpus', required=True)
    ap.add_argument('--verified', default='')
    ap.add_argument('--out', required=True)
    ap.add_argument('--md')
    a = ap.parse_args()
    recs = {}
    for p in a.corpus.split(','):
        recs.update(load_corpus(p))
    vals = {}
    for p in [x for x in a.verified.split(',') if x]:
        for r in json.load(open(p)):
            vals[r['id']] = r.get('value', '')
    res = defaultdict(dict)
    for p in a.runs.split(','):
        for line in open(p):
            r = json.loads(line)
            res[r['id']][r['mode']] = (r['status'], r.get('result', ''))
    out = []
    for rid, rec in recs.items():
        st = res.get(rid, {})
        s, n = st.get('S', ('MISSING', ''))[0], st.get('N', ('MISSING', ''))[0]
        if 'WRONG' in (s, n) or 'NONNUMERIC' in (s, n) or 'CRASH' in (s, n):
            bucket, why = 'BUG', '%s/%s' % (s, n)
        elif s == 'CORRECT':
            bucket, why = 'CAN-DO', ''
        else:
            fam, why = classify(rec['kind'], rec['f'], rec['lim'], vals.get(rid, ''))
            bucket = ('INSTANCE-ONLY:' + fam) if n == 'CORRECT' else fam
        out.append({'id': rid, 'kind': rec['kind'], 'S': s, 'N': n, 'bucket': bucket,
                    'why': why, 'f': rec['f'], 'lim': rec['lim'], 'value': vals.get(rid, '')})
    json.dump(out, open(a.out, 'w'), indent=1)
    c = Counter(r['bucket'] for r in out)
    print('records:', len(out))
    for k, v in c.most_common():
        print('  %-40s %d' % (k, v))
    if a.md:
        L = ['| bucket | records |', '|---|---|'] + ['| %s | %d |' % kv for kv in c.most_common()]
        open(a.md, 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    main()
