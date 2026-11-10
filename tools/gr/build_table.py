#!/usr/bin/env python3
"""
build_table.py -- emit the lazily loaded G&R integral table
(src/internal/GRTables/*.m, read by src/internal/GRTable.m).

Usage:
    python3 -I tools/gr/build_table.py --batch <entries.json>:<pychecked.json>[,...] \
        [--runs run1.jsonl[,run2.jsonl...]] --names <names.txt> --out src/internal/GRTables

Inputs:
    entries.json     tools/gr/merge.py output (the book's side conditions)
    pychecked.json   tools/gr/pycheck.py output: the Mathilda oracle verdict plus
                     the Python mpmath+scipy verdict, with the verified value
    run*.jsonl       Mathilda Integrate runs over the corpora (the
                     tests/test_mathar_corpus.c MATHAR_CORPUS_JSON lines)
    names.txt        Mathilda Names["*"]: every head of a record must exist

A formula becomes a table record when all of these hold:
  - the Mathilda oracle verified it, either strictly (no NIntegrate message)
    or loosely with the Python cross-check also verifying it, and no Python
    probe contradicts the book;
  - Integrate does not already return a correct closed form for the symbolic
    integral (S mode).  A formula Integrate already handles stays with the
    method that handles it, so the table never replaces a general algorithm.
    Without --runs this filter is skipped (only for testing the table engine);
  - its integrand, value and conditions use only heads Mathilda defines;
  - the value has no Sum, Product, Integrate or Limit (re-evaluated on every
    use, which costs more than looking it up), and no test of its own
    parameters (If, Boole, OddQ, ...), which a symbolic query would decide on
    the bare symbol;
  - an antiderivative follows Mathilda's convention Integrate[1/x, x] = Log[x]:
    the book's Log[Abs[u]] is rewritten Log[u] (tools/gr/logconv.py);
  - its conditions survive the condition audit (tools/gr/condaudit.py), with
    the corrections --cond-fix supplies (derived from first principles:
    tools/gr/cond_fixes.tsv by hand, condaudit's .fix output automatically);
  - the book's value passes the corpus runner's checks AS STATED (--book-wrong,
    from build_corpus.py --book): an incorrect G&R formula is corrected from
    first principles or left out, never served.

Output files are split by interval (see GRTable.m): Indefinite.m,
DefHalfLine.m, DefWholeLine.m, DefUnit.m, DefTrig.m, DefOther.m.
"""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from logconv import convert_record  # noqa: E402  (Log[Abs[u]] -> Log[u])

NON_HEADS = {'Hold', 'List', 'Rule', 'Integrate'}
BANNED_VALUE_HEADS = {'Sum', 'Product', 'Integrate', 'Limit', 'NIntegrate', 'D', 'Derivative',
                      # a value that tests its own parameters decides the test on
                      # the bare symbol when the query is symbolic: OddQ[m] is
                      # False for a symbol m, so G&R 2.421.1's Boole[OddQ[m] && ..]
                      # Log term vanished -- a wrong antiderivative
                      'If', 'Boole', 'OddQ', 'EvenQ', 'IntegerQ', 'Which', 'Switch',
                      'PrimeQ', 'Positive', 'Negative', 'NonNegative', 'NonPositive'}
TRIG_END = re.compile(r'^-?(0|(\d+\*?)?Pi(/\d+)?)$')

FILES = ['Indefinite.m', 'DefHalfLine.m', 'DefWholeLine.m', 'DefUnit.m', 'DefTrig.m', 'DefOther.m']


def heads(s):
    return set(re.findall(r'([A-Z][A-Za-z0-9]*)\[', s or ''))


def file_for(kind, limits):
    if kind == 'indefinite':
        return 'Indefinite.m'
    lim = (limits or '').replace(' ', '')
    if lim == '{0,Infinity}':
        return 'DefHalfLine.m'
    if lim == '{-Infinity,Infinity}':
        return 'DefWholeLine.m'
    if lim == '{0,1}':
        return 'DefUnit.m'
    m = re.match(r'^\{(.*),(.*)\}$', lim)
    if m and TRIG_END.match(m.group(1)) and TRIG_END.match(m.group(2)):
        return 'DefTrig.m'
    return 'DefOther.m'


def read_cond_fix(paths):
    """id -> [extra condition, ...] from id<TAB>condition files (# comments)."""
    fix = defaultdict(list)
    for p in filter(None, paths.split(',')):
        for line in open(p):
            line = line.split('#', 1)[0].rstrip()
            if line.strip():
                rid, _, c = line.partition('\t')
                fix[rid.strip()].append(c.strip())
    return fix


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--batch', required=True)
    ap.add_argument('--runs', default='')
    ap.add_argument('--names', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--mathilda', default='./Mathilda')
    ap.add_argument('--book-wrong', default='', help='id<TAB>reason lines to leave out (comma-separated files)')
    ap.add_argument('--cond-fix', default='',
                    help='id<TAB>condition lines (comma-separated files): a book condition '
                         'corrected from first principles, conjoined to the stated ones')
    a = ap.parse_args()

    names = set(l.strip() for l in open(a.names) if l.strip())
    s_status = {}
    for p in filter(None, a.runs.split(',')):
        for line in open(p):
            try:
                r = json.loads(line)
            except ValueError:
                continue
            if r.get('mode') == 'S':
                s_status[r['id']] = r.get('status')

    book_wrong = set()
    for p in filter(None, a.book_wrong.split(',')):
        book_wrong |= set(l.split('\t')[0] for l in open(p) if l.strip())
    cond_fix = read_cond_fix(a.cond_fix)
    out = defaultdict(list)
    why = Counter()
    seen = set()      # batches overlap at their page boundaries
    for pair in a.batch.split(','):
        ep, vp = pair.split(':')
        entries = {e['id']: e for e in json.load(open(ep))}
        for v in json.load(open(vp)):
            e = entries.get(v['id'])
            if e is None or e.get('status') != 'ok' or v['id'] in seen:
                continue
            seen.add(v['id'])
            if v['id'] in book_wrong:
                why['book-wrong'] += 1
                continue
            m_ok = v['status'] in ('VERIFIED', 'VERIFIED1')
            p_ok = v.get('py_status') in ('VERIFIED', 'VERIFIED1')
            py_disagree = (v.get('py_status') == 'MISMATCH' or
                           any(q.get('status') == 'DISAGREE' for q in v.get('py_probes', [])))
            if not m_ok:
                why['not-verified'] += 1
                continue
            if py_disagree:
                why['python-disagrees'] += 1
                continue
            if v.get('tier') == 'loose' and not p_ok:
                why['loose-unconfirmed'] += 1
                continue
            if a.runs:
                st = s_status.get(v['id'])
                if st is None:
                    why['not-run'] += 1
                    continue
                if st == 'CORRECT':
                    why['integrate-does-it'] += 1
                    continue
            v = convert_record(v)     # Integrate[1/x, x] is Log[x]
            f, val = v['integrand'], v['value']
            conds = e.get('conds', []) + cond_fix.get(v['id'], [])
            hs = heads(f) | heads(val) | heads(' '.join(conds))
            if hs & BANNED_VALUE_HEADS:
                why['value-sum-or-integral'] += 1
                continue
            missing = sorted(h for h in hs if h not in names and h not in NON_HEADS)
            if missing:
                why['missing-head'] += 1
                continue
            kind = v['kind']
            lim = 'None' if kind == 'indefinite' else v['limits'].strip()
            body = '"%s", %s, %s, %s, Hold[%s], {%s}, Hold[{%s}]' % (
                v['id'], f, v['var'], lim, val, ', '.join(v.get('params', [])), ', '.join(conds))
            cmt = '(* %s p%s %s *)' % (e['gr'], e.get('page', '?'),
                                      (e.get('cite') or '').replace('*)', '* )'))
            out[file_for(kind, v.get('limits'))].append((e.get('key', v['id']), v['id'], cmt, body))
            why['admitted'] += 1

    keys = dispatch_keys(a.mathilda, [r for lst in out.values() for r in lst])

    os.makedirs(a.out, exist_ok=True)
    for fn in FILES:
        recs = []
        for k, rid, cmt, body in sorted(out.get(fn, [])):
            if rid not in keys:
                why['no-key'] += 1
                continue
            recs.append((k, '%s\nGRTable`Private`add["%s", %s];' % (cmt, keys[rid], body)))
        txt = ('(* %s -- generated by tools/gr/build_table.py; DO NOT EDIT BY HAND.\n'
               '   Gradshteyn & Ryzhik, Table of Integrals, Series, and Products,\n'
               '   7th ed. (2007).  %d records, each verified numerically by the\n'
               '   G&R oracle (tools/gr/verify.m).  Read by src/internal/GRTable.m;\n'
               '   the parameter symbols live in their own context GRTable`R`. *)\n\n'
               '%s\n' % (fn, len(recs), HEADER))
        txt += '\n'.join(r for _, r in recs)
        txt += '\n' + FOOTER
        open(os.path.join(a.out, fn), 'w').write(txt)
        print('%-16s %5d records' % (fn, len(recs)))
    print('selection:', dict(why))


# The records' symbols (the book's parameter names a, b, s, nu, and x) live
# in a context of their own, GRTable`R`, not in the engine's
# GRTable`Private`: a Do iterator of the engine is not renamed the way a
# Module local is, and tryRecord's {s, sols} was the same symbol as G&R
# 17.13.30's parameter s -- s's value then contained s and the evaluation ran
# into $RecursionLimit.
HEADER = ('GRTable`Private`savedPath = $ContextPath;\n$ContextPath = {"System`"};\n'
          'Begin["GRTable`R`"];\n')
FOOTER = '\nEnd[];\n$ContextPath = GRTable`Private`savedPath;\n'


def dispatch_keys(mathilda, recs):
    """id -> dispatch key, computed by GRTable.m's own keyOf so that the key a
    record is filed under is exactly the key a query computes."""
    fd, raw = tempfile.mkstemp(suffix='.m')
    with os.fdopen(fd, 'w') as fh:
        fh.write(HEADER)
        fh.write('GRTable`Private`addraw[id_String, f_, x_Symbol, lim_, ___] := '
                 'Print["@@K\\t", id, "\\t", GRTable`Private`keyOf['
                 'GRTable`Private`splitConst[f, x][[2]], x, lim]];\n')
        for _, rid, _, body in recs:
            fh.write('GRTable`Private`addraw[%s];\n' % body)
        fh.write(FOOTER)
    fd, drv = tempfile.mkstemp(suffix='.m')
    with os.fdopen(fd, 'w') as fh:
        fh.write('Mathilda`LoadModule["GRTable.m"];\nGet["%s"];\n' % raw)
    p = subprocess.run([mathilda, '-file', drv], capture_output=True, text=True, timeout=3600)
    os.unlink(raw)
    os.unlink(drv)
    keys = {}
    for line in p.stdout.splitlines():
        if line.startswith('@@K\t'):
            _, rid, key = line.split('\t', 2)
            keys[rid] = key
    return keys


if __name__ == '__main__':
    main()
