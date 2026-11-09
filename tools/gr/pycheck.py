#!/usr/bin/env python3
"""
pycheck.py -- independent Python cross-check of the G&R oracle and of
Mathilda's numerics (mpmath + scipy; no Mathilda in the loop).

Usage (needs a Python with mpmath and scipy; autodetected, or GR_PYTHON):
    python3.11 -I tools/gr/pycheck.py <verified.json> <pychecked.json> [--jobs 8]

tools/gr/verify.m prints, for every record, up to three admissible parameter
instantiations ("probes") in FullForm: the instantiated integrand, limits and
the book's value, together with Mathilda's N[value] and its NIntegrate verdict.
This script re-decides every probe without Mathilda:

  py_value   the book's value, evaluated by mpmath (30 digits; nsum for series);
  py_int     definite: the integral, by mpmath quad (tanh-sinh; quadosc on an
             infinite range with an oscillatory kernel: tanh-sinh head + quadosc
             tail, at two split points that must agree) and, as a second
             opinion, scipy.integrate.quad; reliable when the two agree to 1e-6
             (or only mpmath converged with a small error estimate / agreeing
             splits);
  py_deriv   indefinite: |d/dx value - f| / max(1, |f|) at x = 0.31, 0.77, 1.43
             (mpmath numerical differentiation).

Per probe it records whether the BOOK agrees with Python, whether Mathilda's
N[value] agrees with mpmath's value of the same expression (a Mathilda special-
function evaluation check), and whether Mathilda's NIntegrate number agrees
with Python's integral (a Mathilda NIntegrate check).  Per record:

  py_status   VERIFIED (book = Python integral on >= 2 probes, or the only one),
              MISMATCH (Python integral reliable and != book on every probe),
              UNDECIDED (Python could not evaluate or integrate reliably),
              UNSUPPORTED (a head this evaluator does not map).
  bugs        list of Mathilda numeric defects found on the probes:
              NINT_WRONG   NIntegrate returned a number != the Python integral
                           (book and Python agree);
              N_WRONG      Mathilda's N[value] != mpmath's value of it;
              NOVALUE      Mathilda could not numericize a value mpmath can.

The output keeps every input field and adds py_status, py_probes and bugs;
tools/gr/build_corpus.py combines it with the Mathilda verdict.
"""
import argparse
import json
import math
import os
import re
import signal
import sys
from concurrent.futures import ProcessPoolExecutor

import mpmath as mp

try:
    import scipy.integrate as si
except ImportError:  # pragma: no cover - scipy is required, reported clearly
    sys.exit('pycheck.py needs scipy: run it with a Python that has mpmath and scipy (GR_PYTHON)')

mp.mp.dps = 20
TOL = mp.mpf('1e-6')
PROBE_SECONDS = 60


# ---------------------------------------------------------------- FullForm parser
TOK = re.compile(r'\s*(?:(?P<num>-?\d+\.\d*(?:`[\d.]*)?(?:\*\^-?\d+)?|-?\d+(?:`[\d.]*)?(?:\*\^-?\d+)?)'
                 r'|(?P<str>"(?:[^"\\]|\\.)*")|(?P<sym>[A-Za-z$][A-Za-z0-9$`]*)|(?P<p>[\[\],]))')


class Unsupported(Exception):
    pass


def parse(text):
    """FullForm text -> nested tuples ('call', head, [args]) / ('num', str) /
    ('sym', name) / ('str', s)."""
    toks = []
    pos = 0
    text = text.strip()
    while pos < len(text):
        m = TOK.match(text, pos)
        if not m or m.end() == pos:
            raise ValueError('cannot tokenize at %r' % text[pos:pos + 30])
        pos = m.end()
        for k in ('num', 'str', 'sym', 'p'):
            if m.group(k) is not None:
                toks.append((k, m.group(k)))
                break
    i = 0

    def expr():
        nonlocal i
        k, v = toks[i]
        i += 1
        if k == 'num':
            node = ('num', v)
        elif k == 'str':
            node = ('str', v[1:-1])
        elif k == 'sym':
            node = ('sym', v)
        else:
            raise ValueError('unexpected %r' % v)
        while i < len(toks) and toks[i] == ('p', '['):   # head[...][...] (Derivative[..][f][..])
            i += 1
            args = []
            if toks[i] != ('p', ']'):
                while True:
                    args.append(expr())
                    if toks[i] == ('p', ','):
                        i += 1
                        continue
                    break
            if toks[i] != ('p', ']'):
                raise ValueError('expected ]')
            i += 1
            node = ('call', node, args)
        return node
    node = expr()
    if i != len(toks):
        raise ValueError('trailing input')
    return node


def num(s):
    s = re.sub(r'`[\d.]*', '', s)
    if '*^' in s:
        a, b = s.split('*^')
        return mp.mpf(a) * mp.mpf(10) ** int(b)
    return mp.mpf(s) if '.' in s else mp.mpf(int(s))


# ---------------------------------------------------------------- evaluator
CONST = {'Pi': lambda: mp.pi, 'E': lambda: mp.e, 'EulerGamma': lambda: mp.euler,
         'Catalan': lambda: mp.catalan, 'GoldenRatio': lambda: mp.phi, 'I': lambda: mp.j,
         'Infinity': lambda: mp.inf, 'True': lambda: True, 'False': lambda: False}


def _pfq(a, b, z):
    return mp.hyper(list(a), list(b), z)


def _legendre_p(*a):
    if len(a) == 2:
        return mp.legendre(a[0], a[1])
    if len(a) == 3:
        return mp.legenp(a[0], a[1], a[2], type=2)
    return mp.legenp(a[0], a[1], a[3], type=int(a[2]))


def _legendre_q(*a):
    if len(a) == 2:
        return mp.legenq(a[0], 0, a[1], type=2)
    if len(a) == 3:
        return mp.legenq(a[0], a[1], a[2], type=2)
    return mp.legenq(a[0], a[1], a[3], type=int(a[2]))


def _gamma(*a):
    if len(a) == 1:
        return mp.gamma(a[0])
    if len(a) == 2:
        return mp.gammainc(a[0], a[1], mp.inf)
    return mp.gammainc(a[0], a[1], a[2])


def _beta(*a):
    if len(a) == 2:
        return mp.beta(a[0], a[1])
    return mp.betainc(a[1], a[2], 0, a[0])


def _polygamma(*a):
    return mp.digamma(a[0]) if len(a) == 1 else mp.psi(int(a[0]), a[1])


def _log(*a):
    return mp.log(a[0]) if len(a) == 1 else mp.log(a[1]) / mp.log(a[0])


def _arctan(*a):
    return mp.atan(a[0]) if len(a) == 1 else mp.arg(a[0] + mp.j * a[1])


def _zeta(*a):
    return mp.zeta(a[0]) if len(a) == 1 else mp.zeta(a[0], a[1])


def _ellipf(phi, m):
    return mp.ellipf(phi, m)


def _ellipe(*a):
    return mp.ellipe(a[0]) if len(a) == 1 else mp.ellipe(a[0], a[1])


def _ellippi(*a):
    return mp.ellippi(a[0], a[1]) if len(a) == 2 else mp.ellippi(a[0], a[1], a[2])


def _sign(z):
    return mp.sign(z)


def _real(z):
    return mp.re(z)


FUN = {
    'Exp': mp.exp, 'Log': _log, 'Sqrt': mp.sqrt, 'Abs': abs, 'Sign': _sign,
    'Re': _real, 'Im': mp.im, 'Arg': mp.arg, 'Conjugate': mp.conj,
    'Sin': mp.sin, 'Cos': mp.cos, 'Tan': mp.tan, 'Cot': mp.cot, 'Sec': mp.sec, 'Csc': mp.csc,
    'Sinh': mp.sinh, 'Cosh': mp.cosh, 'Tanh': mp.tanh, 'Coth': mp.coth, 'Sech': mp.sech, 'Csch': mp.csch,
    'ArcSin': mp.asin, 'ArcCos': mp.acos, 'ArcTan': _arctan, 'ArcCot': mp.acot, 'ArcSec': mp.asec,
    'ArcCsc': mp.acsc, 'ArcSinh': mp.asinh, 'ArcCosh': mp.acosh, 'ArcTanh': mp.atanh,
    'ArcCoth': mp.acoth, 'ArcSech': mp.asech, 'ArcCsch': mp.acsch,
    'Gamma': _gamma, 'LogGamma': mp.loggamma, 'Beta': _beta, 'PolyGamma': _polygamma,
    'Zeta': _zeta, 'HurwitzZeta': lambda s, a: mp.zeta(s, a), 'PolyLog': mp.polylog,
    'LerchPhi': mp.lerchphi, 'Erf': mp.erf, 'Erfc': mp.erfc, 'Erfi': mp.erfi,
    'ExpIntegralEi': mp.ei, 'ExpIntegralE': mp.expint, 'SinIntegral': mp.si, 'CosIntegral': mp.ci,
    'SinhIntegral': mp.shi, 'CoshIntegral': mp.chi, 'LogIntegral': mp.li,
    'FresnelS': mp.fresnels, 'FresnelC': mp.fresnelc,
    'BesselJ': mp.besselj, 'BesselY': mp.bessely, 'BesselI': mp.besseli, 'BesselK': mp.besselk,
    'HankelH1': mp.hankel1, 'HankelH2': mp.hankel2, 'StruveH': mp.struveh, 'StruveL': mp.struvel,
    'AngerJ': mp.angerj, 'WeberE': mp.webere, 'KelvinBer': mp.ber, 'KelvinBei': mp.bei,
    'KelvinKer': mp.ker, 'KelvinKei': mp.kei, 'AiryAi': mp.airyai, 'AiryBi': mp.airybi,
    'SphericalBesselJ': lambda n, z: mp.sqrt(mp.pi / (2 * z)) * mp.besselj(n + mp.mpf(1) / 2, z),
    'SphericalBesselY': lambda n, z: mp.sqrt(mp.pi / (2 * z)) * mp.bessely(n + mp.mpf(1) / 2, z),
    'EllipticK': mp.ellipk, 'EllipticE': _ellipe, 'EllipticF': _ellipf, 'EllipticPi': _ellippi,
    'JacobiSN': lambda u, m: mp.ellipfun('sn', u, m=m), 'JacobiCN': lambda u, m: mp.ellipfun('cn', u, m=m),
    'JacobiDN': lambda u, m: mp.ellipfun('dn', u, m=m),
    'Hypergeometric0F1': mp.hyp0f1, 'Hypergeometric1F1': mp.hyp1f1, 'Hypergeometric2F1': mp.hyp2f1,
    'HypergeometricU': mp.hyperu, 'HypergeometricPFQ': _pfq,
    'WhittakerM': mp.whitm, 'WhittakerW': mp.whitw, 'ParabolicCylinderD': mp.pcfd,
    'LegendreP': _legendre_p, 'LegendreQ': _legendre_q, 'LaguerreL': lambda *a: mp.laguerre(a[0], 0 if len(a) == 2 else a[1], a[-1]),
    'HermiteH': mp.hermite, 'ChebyshevT': mp.chebyt, 'ChebyshevU': mp.chebyu, 'JacobiP': mp.jacobi,
    'GegenbauerC': lambda *a: mp.gegenbauer(a[0], a[1], a[2]) if len(a) == 3 else (2 / a[0]) * mp.chebyt(a[0], a[1]),
    'Factorial': mp.factorial, 'Factorial2': mp.fac2, 'Binomial': mp.binomial, 'Pochhammer': mp.rf,
    'BernoulliB': lambda *a: mp.bernoulli(int(a[0])) if len(a) == 1 else mp.bernpoly(int(a[0]), a[1]),
    'EulerE': lambda *a: mp.eulernum(int(a[0])) if len(a) == 1 else mp.eulerpoly(int(a[0]), a[1]),
    'Floor': mp.floor, 'Ceiling': mp.ceil, 'Max': max, 'Min': min,
}


def ev(node, env):
    """Evaluate a parsed FullForm node with mpmath.  env maps symbol -> value."""
    t = node[0]
    if t == 'num':
        return num(node[1])
    if t == 'str':
        return node[1]
    if t == 'sym':
        name = node[1]
        if name in env:
            return env[name]
        if name in CONST:
            return CONST[name]()
        raise Unsupported('symbol ' + name)
    head, args = node[1], node[2]
    if head[0] == 'call':
        return ev_derivative(head, args, env)
    h = head[1]
    if h == 'List':
        return [ev(a, env) for a in args]
    if h == 'Rational':
        return mp.mpf(ev(args[0], env)) / ev(args[1], env)
    if h == 'Complex':
        return mp.mpc(ev(args[0], env), ev(args[1], env))
    if h == 'Plus':
        s = 0
        for a in args:
            s = s + ev(a, env)
        return s
    if h == 'Times':
        p = 1
        for a in args:
            p = p * ev(a, env)
        return p
    if h == 'Power':
        b, e = ev(args[0], env), ev(args[1], env)
        if b == 0 and mp.re(e) < 0:
            return mp.inf
        return mp.power(b, e)
    if h in ('Sum', 'Product'):
        return ev_sum(h, args, env)
    if h == 'Boole':
        return 1 if ev(args[0], env) else 0
    if h == 'If':
        return ev(args[1], env) if ev(args[0], env) else ev(args[2], env)
    if h in ('Equal', 'Unequal', 'Less', 'LessEqual', 'Greater', 'GreaterEqual', 'And', 'Or', 'Not',
             'OddQ', 'EvenQ'):
        return ev_logic(h, args, env)
    if h == 'MeijerG':
        a = ev(args[0], env)
        b = ev(args[1], env)
        return mp.meijerg(a, b, ev(args[2], env))
    if h == 'DirectedInfinity':
        return mp.inf if not args else mp.inf * ev(args[0], env)
    if h in FUN:
        return FUN[h](*[ev(a, env) for a in args])
    raise Unsupported('head ' + h)


def ev_logic(h, args, env):
    v = [ev(a, env) for a in args]
    if h == 'Equal':
        return all(v[0] == x for x in v[1:])
    if h == 'Unequal':
        return v[0] != v[1]
    if h == 'Less':
        return all(v[i] < v[i + 1] for i in range(len(v) - 1))
    if h == 'LessEqual':
        return all(v[i] <= v[i + 1] for i in range(len(v) - 1))
    if h == 'Greater':
        return all(v[i] > v[i + 1] for i in range(len(v) - 1))
    if h == 'GreaterEqual':
        return all(v[i] >= v[i + 1] for i in range(len(v) - 1))
    if h == 'And':
        return all(v)
    if h == 'Or':
        return any(v)
    if h == 'Not':
        return not v[0]
    if h == 'OddQ':
        return int(v[0]) % 2 == 1
    return int(v[0]) % 2 == 0


def ev_sum(h, args, env):
    body = args[0]
    it = args[1]
    if it[0] != 'call' or it[1] != ('sym', 'List'):
        raise Unsupported('iterator')
    var = it[2][0][1]
    lo = ev(it[2][1], env) if len(it[2]) > 2 else 1
    hi = ev(it[2][-1], env)

    def term(k):
        e = dict(env)
        e[var] = k
        return ev(body, e)
    if h == 'Sum':
        if hi == mp.inf:
            return mp.nsum(term, [lo, mp.inf])
        s = 0
        for k in range(int(lo), int(hi) + 1):
            s += term(mp.mpf(k))
        return s
    p = 1
    if hi == mp.inf:
        return mp.nprod(term, [lo, mp.inf])
    for k in range(int(lo), int(hi) + 1):
        p *= term(mp.mpf(k))
    return p


def ev_derivative(head, args, env):
    """Derivative[n1, n2, ...][f][x1, x2, ...] by mpmath numerical differentiation."""
    orders = [int(ev(a, env)) for a in head[2]]
    fnode = head[1]
    if fnode[0] != 'call' or fnode[1] != ('sym', 'Derivative'):
        raise Unsupported('compound head')
    fname = fnode[2][0][1] if fnode[2] and fnode[2][0][0] == 'sym' else None
    if fname not in FUN:
        raise Unsupported('derivative of ' + str(fname))
    pt = [ev(a, env) for a in args]
    return mp.diff(lambda *z: FUN[fname](*z), pt, tuple(orders))


# ---------------------------------------------------------------- integration
OSC = {'Sin', 'Cos', 'BesselJ', 'BesselY', 'StruveH', 'AngerJ', 'WeberE', 'HankelH1', 'HankelH2'}


def omega(node, var):
    """Largest |linear frequency| of an oscillatory factor in var, or None."""
    best = None

    def walk(n):
        nonlocal best
        if n[0] != 'call':
            return
        h = n[1][1] if n[1][0] == 'sym' else None
        if h in OSC or h == 'Power' or h == 'Exp':
            arg = n[2][-1]
            if h in ('Power',):
                if n[2][0] != ('sym', 'E'):
                    arg = None
            if arg is not None:
                try:
                    a0 = ev(arg, {var: mp.mpf(0)})
                    a1 = ev(arg, {var: mp.mpf(1)})
                    a2 = ev(arg, {var: mp.mpf(2)})
                    w = a1 - a0
                    if abs((a2 - a1) - w) < 1e-12:   # linear in var
                        w = abs(mp.im(w)) if h in ('Power', 'Exp') else abs(w)
                        if w > 0 and (best is None or w > best):
                            best = w
                except Exception:
                    pass
        for c in n[2]:
            walk(c)
    walk(node)
    return best


def py_integral(fnode, var, lo, hi):
    """(value, reliable, how) -- mpmath, with scipy as a second opinion."""
    f = lambda t: ev(fnode, {var: t})
    infinite = lo == -mp.inf or hi == mp.inf
    w = omega(fnode, var) if infinite else None
    vm = em = None
    osc_ok = False
    try:
        if w:
            # quadosc alone drifted by 1e-3..1e-5 on integrands singular at the
            # finite end (Y(a/x), log x) or oscillating in sqrt(x): integrate a
            # finite head with tanh-sinh, hand only the tail to quadosc, and do
            # it at two different split points -- reliable only if they agree.
            def osc_split(n_periods):
                A = n_periods * 2 * mp.pi / w

                def half(a, b_inf):     # [a, inf) with a finite
                    head = mp.quad(f, [a, a + A], maxdegree=8)
                    return head + mp.quadosc(f, [a + A, mp.inf], omega=w)
                if lo == -mp.inf and hi == mp.inf:
                    g = lambda t: f(-t)
                    hneg = mp.quad(g, [0, A], maxdegree=8) + mp.quadosc(g, [A, mp.inf], omega=w)
                    return half(0, True) + hneg
                if hi == mp.inf:
                    return half(lo, True)
                g = lambda t: f(-t)
                return mp.quad(g, [-hi, -hi + A], maxdegree=8) + mp.quadosc(g, [-hi + A, mp.inf], omega=w)
            v1 = osc_split(6)
            v2 = osc_split(11)
            vm = v1
            osc_ok = abs(v1 - v2) <= mp.mpf('1e-9') * max(1, abs(v1))
            em = None
        else:
            vm, em = mp.quad(f, [lo, hi], error=True, maxdegree=8)
    except Exception:
        vm = None
    vs = None
    try:
        def fr(t):
            return float(mp.re(f(mp.mpf(t))))

        def fi(t):
            return float(mp.im(f(mp.mpf(t))))
        a = float(lo) if lo != -mp.inf else -math.inf
        b = float(hi) if hi != mp.inf else math.inf
        kw = {'limit': 400}
        if w and (math.isinf(a) or math.isinf(b)) and not (math.isinf(a) and math.isinf(b)):
            pass   # scipy's QAWF needs a pure sin/cos weight; skip, mpmath quadosc decides
        else:
            r = si.quad(fr, a, b, **kw)[0]
            im = si.quad(fi, a, b, **kw)[0]
            vs = mp.mpc(r, im)
    except Exception:
        vs = None
    if vm is not None and vs is not None:
        if abs(vm - vs) <= mp.mpf('1e-6') * max(1, abs(vm)):
            return vm, True, 'mpmath+scipy'
        return vm, False, 'mpmath!=scipy'
    if vm is not None:
        ok = (em is not None and em < mp.mpf('1e-10') * max(1, abs(vm))) or (w and osc_ok)
        return vm, bool(ok), 'mpmath' + ('-osc' if w else '')
    if vs is not None:
        return vs, False, 'scipy-only'
    return None, False, 'none'


def close(a, b, tol=TOL):
    return abs(a - b) <= tol * max(1, abs(b))


def mathilda_number(s):
    """Mathilda InputForm number (possibly complex) -> mpc, else None."""
    if not s or s in ('novalue', '""'):
        return None
    t = re.sub(r'`[\d.]*', '', s).replace('*^', 'e').replace('*I', 'j').replace(' ', '')
    t = t.replace('I', '1j') if t in ('I', '-I') else t
    try:
        return mp.mpc(complex(t))
    except ValueError:
        return None


def judge_number(s):
    """Mathilda grJudge verdict string, e.g. {good, 0.4, strict} -> (kind, number)."""
    m = re.match(r'\{\s*"?(good|bad|skip)"?\s*(?:,\s*([^,}]+))?', s or '')
    if not m:
        return None, None
    return m.group(1), (mathilda_number(m.group(2)) if m.group(2) else None)


def check_probe(p):
    out = {'inst': p['inst']}
    try:
        fnode = parse(p['f'])
        vnode = parse(p['value'])
        var = p['var']
    except ValueError as ex:
        out['status'] = 'PARSE'
        out['why'] = str(ex)
        return out
    try:
        if p['lo'] == '' and p['hi'] == '':     # indefinite: the derivative test
            errs = []
            for x0 in ('0.31', '0.77', '1.43'):
                x = mp.mpf(x0)
                dv = mp.diff(lambda t: ev(vnode, {var: t}), x)
                fv = ev(fnode, {var: x})
                errs.append(abs(dv - fv) / max(1, abs(fv)))
            e = max(errs)
            out['py_deriv'] = mp.nstr(e, 5)
            out['status'] = 'AGREE' if e < mp.mpf('1e-8') else 'DISAGREE'
            return out
        val = ev(vnode, {})
        out['py_value'] = mp.nstr(val, 20)
        mN = mathilda_number(p.get('mathilda_N'))
        if p.get('mathilda_N') == 'novalue':
            out['bug'] = 'NOVALUE'
        elif mN is not None and not close(mN, val, mp.mpf('1e-8')):
            out['bug'] = 'N_WRONG'
        lo, hi = ev(parse(p['lo']), {}), ev(parse(p['hi']), {})
        iv, reliable, how = py_integral(fnode, var, lo, hi)
        out['py_int'] = mp.nstr(iv, 20) if iv is not None else ''
        out['py_how'] = how
        if iv is None or not reliable:
            out['status'] = 'UNRELIABLE'
            return out
        if abs(iv) < 1e-6:
            out['status'] = 'TINY'
            return out
        out['status'] = 'AGREE' if close(iv, val) else 'DISAGREE'
        kind, mnum = judge_number(p.get('mathilda_judge'))
        if mnum is not None and not close(mnum, iv, mp.mpf('1e-5')) and out['status'] == 'AGREE':
            out['bug'] = 'NINT_WRONG'
        return out
    except Unsupported as ex:
        out['status'] = 'UNSUPPORTED'
        out['why'] = str(ex)
        return out
    except Exception as ex:
        out['status'] = 'ERROR'
        out['why'] = '%s: %s' % (type(ex).__name__, str(ex)[:120])
        return out


class ProbeTimeout(BaseException):
    """BaseException, not Exception: mpmath's inner retry loops (and the generic
    handler in check_probe) catch Exception, which swallowed a one-shot alarm and
    left workers spinning for minutes.  The timer also re-fires every 2 s until
    cleared, so a handler that catches it once cannot stall the probe."""
    pass


def _alarm(signum, frame):
    raise ProbeTimeout()


def check_record(r):
    signal.signal(signal.SIGALRM, _alarm)
    probes = []
    for p in r.get('probes', []):
        signal.setitimer(signal.ITIMER_REAL, PROBE_SECONDS, 2)
        try:
            probes.append(check_probe(p))
        except ProbeTimeout:
            probes.append({'inst': p['inst'], 'status': 'TIMEOUT'})
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
    st = [p['status'] for p in probes]
    agree = st.count('AGREE')
    if any(s == 'UNSUPPORTED' for s in st) and agree == 0:
        verdict = 'UNSUPPORTED'
    elif agree >= 2 or (agree == 1 and len(probes) == 1):
        verdict = 'VERIFIED'
    elif agree == 0 and st.count('DISAGREE') >= 1 and all(s in ('DISAGREE', 'TINY') for s in st):
        verdict = 'MISMATCH'
    elif agree == 1:
        verdict = 'VERIFIED1'
    else:
        verdict = 'UNDECIDED'
    bugs = sorted(set(p['bug'] for p in probes if p.get('bug')))
    return dict(r, py_status=verdict, py_probes=probes, bugs=bugs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('verified')
    ap.add_argument('out')
    ap.add_argument('--jobs', type=int, default=8)
    a = ap.parse_args()
    recs = json.load(open(a.verified))
    # Resumable and crash-tolerant: each finished record is appended to
    # <out>.cache.jsonl; a worker that dies (an mpmath/scipy segfault or the
    # OOM killer) breaks the whole pool, so the pool is rebuilt and the record
    # it was holding is marked CRASHED instead of losing the run.
    from concurrent.futures.process import BrokenProcessPool
    cache_path = a.out + '.cache.jsonl'
    done = {}
    if os.path.exists(cache_path):
        for line in open(cache_path):
            try:
                r = json.loads(line)
                done[r['id']] = r
            except ValueError:
                pass
    todo = [r for r in recs if r['id'] not in done]
    print('pycheck: %d records, %d cached, %d to run' % (len(recs), len(done), len(todo)), flush=True)
    cache = open(cache_path, 'a')
    while todo:
        ex = ProcessPoolExecutor(a.jobs)
        futs = {ex.submit(check_record, r): r for r in todo}
        broken = False
        progress = 0
        for fu in futs:
            r = futs[fu]
            try:
                res = fu.result()
            except BrokenProcessPool:
                broken = True
                continue
            done[r['id']] = res
            progress += 1
            cache.write(json.dumps(res) + '\n')
            cache.flush()
        ex.shutdown(cancel_futures=True)
        left = [r for r in todo if r['id'] not in done]
        if broken and left and progress == 0:
            # the pool died before finishing anything: one of the first
            # records is the culprit -- run the head of the queue one per
            # single-worker pool, recording a record that kills it as CRASHED
            for r in left[:a.jobs]:
                ex1 = ProcessPoolExecutor(1)
                try:
                    res = ex1.submit(check_record, r).result()
                except BrokenProcessPool:
                    res = dict(r, py_status='UNDECIDED', py_probes=[], bugs=[], py_crash=True)
                ex1.shutdown(cancel_futures=True)
                done[r['id']] = res
                cache.write(json.dumps(res) + '\n')
                cache.flush()
            left = [r for r in left if r['id'] not in done]
        todo = left
    cache.close()
    out = [done[r['id']] for r in recs]
    json.dump(out, open(a.out, 'w'), indent=1)
    cnt, bugs = {}, {}
    for r in out:
        k = '%s/%s' % (r.get('status'), r['py_status'])
        cnt[k] = cnt.get(k, 0) + 1
        for b in r['bugs']:
            bugs[b] = bugs.get(b, 0) + 1
    print('mathilda/python verdicts:', dict(sorted(cnt.items())))
    print('Mathilda numeric defects found:', bugs)


if __name__ == '__main__':
    main()
