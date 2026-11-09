#!/usr/bin/env python3
"""
build_corpus.py -- emit the Gradshteyn & Ryzhik Integrate corpora.

Usage:
    python3 -I tools/gr/build_corpus.py <entries.json> <verified.json>[,more.json...] \
        <out-dir> [<mathilda names file>]

Inputs:
    entries.json        tools/gr/merge.py output (every transcribed formula)
    verified*.json      tools/gr/pycheck.py output: the Mathilda oracle verdict
                        (tools/gr/verify.m) plus the Python mpmath+scipy verdict;
                        later files override earlier ones by id
    names file          optional: Mathilda `Names["*"]`, one per line; a record
                        whose integrand uses a head Mathilda does not define is
                        skipped as SKIP-UNSUPPORTED-HEAD

Outputs (in <out-dir>):
    GR<ch>IntegralsCorpus.m   one per chapter (2, 3, 4, 5, 6, 7, 17); record
                              shape identical to the Mathar corpora, so
                              tests/test_mathar_corpus.c runs them unchanged
    GRSkipped.md              every formula not in a corpus, with its reason
    GRNIntegrateGaps.md       Mathilda NIntegrate gaps the oracle met: records
                              admitted only on a "loose" instantiation (NIntegrate
                              converged with a message) and records the oracle
                              could not decide (NOCONV) -- an NIntegrate work list

A record is admitted when the Mathilda oracle or the Python cross-check
verifies the book's value and neither contradicts it (a contradiction is
SKIP-CONFLICT).  Which judge admitted it -- mathilda-strict / mathilda-loose
(+python), or python-only -- is recorded next to each record as a comment.
    GRNumericDefects.md       Mathilda numeric defects pycheck found
"""
import importlib.util
import json
import os
import re
import sys
from collections import Counter, defaultdict

_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    'mathar_build', os.path.join(_here, '..', 'mathar', 'build_corpus.py'))
mb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mb)


def chapter(gr):
    return gr.split('.')[0]


def section(gr):
    m = re.match(r'(\d+\.\d\d?)', gr)
    return m.group(1) if m else gr


def main():
    entries = json.load(open(sys.argv[1]))
    ver = {}
    for path in sys.argv[2].split(','):
        for r in json.load(open(path)):
            ver[r['id']] = r
    outdir = sys.argv[3]
    names = None
    if len(sys.argv) > 4:
        names = set(l.strip() for l in open(sys.argv[4]) if l.strip())

    recs = defaultdict(list)
    skipped = []
    gaps = []
    numbugs = []      # Mathilda numeric defects pycheck found (NINT_WRONG, N_WRONG, NOVALUE)
    conflicts = []    # Mathilda and Python disagree on whether the book is right
    for e in entries:
        st = e['status']
        if st != 'ok':
            skipped.append((e, st, e.get('reason', '')))
            continue
        v = ver.get(e['id'])
        if v is None:
            skipped.append((e, 'SKIP-UNVERIFIED', 'not run through the oracle'))
            continue
        m_ok = v['status'] in ('VERIFIED', 'VERIFIED1')
        ps = v.get('py_status', 'UNCHECKED')
        p_ok = ps in ('VERIFIED', 'VERIFIED1')
        for b in v.get('bugs', []):
            numbugs.append((e, v, b))
        if v['status'] == 'NOCONV' or (m_ok and v.get('tier') == 'loose'):
            gaps.append((e, v, 'loose-silentwrong' if v.get('silentwrong') else
                         ('loose' if m_ok else 'NOCONV')))
        # Admission: at least one numeric judge (Mathilda oracle, or Python
        # mpmath+scipy) verifies the book, and none contradicts it.
        # Any instantiation where Python's integral is reliable and disagrees
        # with the book blocks admission: the formula fails at an admissible
        # point (misreading, or a condition the book omits) -- even if another
        # instantiation verified (6.513.2 verified once, loosely, and failed
        # reliably at two others).
        py_disagree = any(q.get('status') == 'DISAGREE' for q in v.get('py_probes', []))
        if (m_ok and (ps == 'MISMATCH' or py_disagree)) or (p_ok and v['status'] == 'MISMATCH'):
            skipped.append((e, 'SKIP-CONFLICT', 'Mathilda %s / Python %s' % (v['status'], ps)))
            conflicts.append((e, v))
            continue
        if not (m_ok or p_ok):
            why = v.get('why') or v.get('bad') or ''
            tag = v['status'] if v['status'] != 'NOCONV' or ps in ('UNCHECKED', 'UNDECIDED') else 'NOCONV-PY-' + ps
            skipped.append((e, 'SKIP-ORACLE-' + tag, str(why)[:200]))
            continue
        f = v['integrand']
        if names is not None:
            missing = sorted(h for h in mb.heads_of(f) if h not in names and h not in mb.NON_HEADS)
            if missing:
                skipped.append((e, 'SKIP-UNSUPPORTED-HEAD', ', '.join(missing)))
                continue
        kind, var = v['kind'], v['var']
        if m_ok:
            good = v['good']
            insts = [g['inst'] for g in good]
            refs = [mb.plain_number(g['ref']) for g in good] if kind == 'definite' else []
            inst_f, inst_lim = v.get('inst_integrand'), v.get('inst_limits')
            judge = 'mathilda-' + v.get('tier', 'strict') + ('+python' if p_ok else '')
        else:
            # Python-only: Mathilda could not decide (NIntegrate failure, or the
            # value does not numericize in Mathilda); the probes Python agreed
            # on supply instantiations, instance and mpmath reference values.
            pr = [(p, q) for p, q in zip(v.get('probes', []), v.get('py_probes', []))
                  if q.get('status') == 'AGREE']
            insts = [p['inst'] for p, q in pr]
            refs = [q['py_value'] for p, q in pr] if kind == 'definite' else []
            inst_f = pr[0][0]['f']
            inst_lim = '{%s, %s}' % (pr[0][0]['lo'], pr[0][0]['hi']) if kind == 'definite' else ''
            judge = 'python-only'
        assum = mb.assumptions_for(e, dict(v, good=[{'inst': x} for x in insts]))
        spec = var if kind == 'indefinite' else '{%s, %s}' % (var, v['limits'].strip()[1:-1])
        opt = ', Assumptions -> %s' % assum if assum else ''
        sym = 'Integrate[%s, %s%s]' % (f, spec, opt)
        if inst_f:
            ispec = var if kind == 'indefinite' else '{%s, %s}' % (var, inst_lim.strip()[1:-1])
            num = 'Integrate[%s, %s]' % (inst_f, ispec)
        else:
            num = 'Integrate[%s, %s]' % (f, spec)
        nv = judge
        rec = '(* %s p%s %s; oracle %s *)\n{"%s", "%s",\n  Hold[%s],\n  Hold[%s],\n  {%s},\n  {%s},\n  Hold[%s], %s}' % (
            e['gr'], e.get('page', '?'), e.get('cite', ''), nv,
            e['id'], kind, sym, num, ', '.join(insts), ', '.join(refs), f, var)
        recs[chapter(e['gr'])].append((e, rec))

    for ch, lst in recs.items():
        lst.sort(key=lambda t: t[0]['key'])
        body, cur = [], None
        for e, rec in lst:
            sec = section(e['gr'])
            if sec != cur:
                title = re.sub(r'^\s*\d+\.\d+\s*', '', e.get('section', ''))
                body.append('(* ---- %s %s ---- *)' % (sec, title.replace('*)', '* )')))
                cur = sec
            body.append(rec)
        head = (
            '(* GR%sIntegralsCorpus.m -- generated by tools/gr/build_corpus.py.\n'
            '   DO NOT EDIT BY HAND: fix tools/gr/overrides.json and regenerate.\n\n'
            '   Source: I. S. Gradshteyn & I. M. Ryzhik, Table of Integrals, Series,\n'
            '   and Products, 7th ed. (2007), chapter %s.  id = GR<formula>.<item>.\n'
            '   Transcribed from the page images (tools/gr/TRANSCRIBE.md); every\n'
            '   reference value comes from the Mathilda oracle tools/gr/verify.m: the\n'
            '   book\'s value had to agree with NIntegrate at two parameter\n'
            '   instantiations (definite) or differentiate back to the integrand\n'
            '   (indefinite).  The comment above each record gives the formula, PDF\n'
            '   page, the book\'s source citation and the oracle tier (strict: no\n'
            '   NIntegrate message; loose: agreement to 1e-5 after a convergence\n'
            '   message).\n\n'
            '   Record: {id, kind, Hold[symbolic], Hold[instance], {inst..}, {ref..},\n'
            '            Hold[integrand], var}   -- runner: tests/test_mathar_corpus.c\n'
            '*)\n' % (ch, ch))
        txt = head + '{\n' + _join(body) + '\n}\n'
        open(os.path.join(outdir, 'GR%sIntegralsCorpus.m' % ch), 'w').write(txt)

    by = defaultdict(list)
    for e, st, why in skipped:
        by[st].append((e, why))
    lines = ['# Gradshteyn & Ryzhik integrals not in the corpus', '',
             'Generated by `tools/gr/build_corpus.py`. Every transcribed formula is either',
             'a corpus record or listed here with the reason it is not one.', '',
             '| Status | Count |', '|---|---|']
    for st in sorted(by):
        lines.append('| %s | %d |' % (st, len(by[st])))
    lines.append('')
    for st in sorted(by):
        lines += ['## %s' % st, '']
        for e, why in sorted(by[st], key=lambda t: t[0]['key']):
            tex = ' '.join((e.get('tex') or e.get('lhs') or '').split())[:160].replace('|', '\\|')
            lines.append('- `%s` (p%s, %s): %s — `%s`' % (
                e['gr'], e.get('page', '?'), e.get('cite', ''), str(why).replace('|', '\\|')[:200], tex))
        lines.append('')
    open(os.path.join(outdir, 'GRSkipped.md'), 'w').write('\n'.join(lines))

    dl = ['# Mathilda NIntegrate gaps met by the G&R oracle', '',
          'Generated by `tools/gr/build_corpus.py`.  `loose`: admitted, but NIntegrate only',
          'reached the book\'s value (to 1e-5) with a convergence message.  `NOCONV`: the',
          'oracle could not decide the reading because NIntegrate never converged.',
          '`loose-silentwrong`: the default NIntegrate returned a WRONG number without a',
          'message (the second strategy matched the book) -- an NIntegrate correctness bug.', '',
          '| id | kind | integrand | limits |', '|---|---|---|---|']
    for e, v, why in sorted(gaps, key=lambda t: t[0]['key']):
        dl.append('| %s | %s | `%s` | `%s` |' % (
            e['id'], why, (v.get('integrand') or e.get('lhs', '')).replace('|', '\\|'),
            v.get('limits', '')))
    open(os.path.join(outdir, 'GRNIntegrateGaps.md'), 'w').write('\n'.join(dl) + '\n')

    bl = ['# Mathilda numeric defects found by the G&R cross-check', '',
          'Generated by `tools/gr/build_corpus.py` from `tools/gr/pycheck.py` (mpmath + scipy).',
          '`NINT_WRONG`: NIntegrate returned a number that disagrees with the Python integral (the',
          'book agrees with Python).  `N_WRONG`: Mathilda\'s N of the book\'s value disagrees with',
          'mpmath.  `NOVALUE`: Mathilda cannot numericize a value mpmath can.  CONFLICT: Mathilda\'s',
          'oracle and Python disagree on the book itself.', '',
          '| id | defect | integrand | limits | value |', '|---|---|---|---|---|']
    for e, v, b in sorted(numbugs, key=lambda t: (t[2], t[0]['key'])):
        bl.append('| %s | %s | `%s` | `%s` | `%s` |' % (
            e['id'], b, (v.get('integrand') or '').replace('|', '\\|'), v.get('limits', ''),
            (v.get('value') or '')[:120].replace('|', '\\|')))
    for e, v in conflicts:
        bl.append('| %s | CONFLICT | `%s` | `%s` | `%s` |' % (
            e['id'], (v.get('integrand') or '').replace('|', '\\|'), v.get('limits', ''),
            (v.get('value') or '')[:120].replace('|', '\\|')))
    open(os.path.join(outdir, 'GRNumericDefects.md'), 'w').write('\n'.join(bl) + '\n')

    print('records per chapter:', {k: len(v) for k, v in sorted(recs.items())},
          'skipped:', len(skipped), dict(Counter(st for _, st, _ in skipped)),
          'NIntegrate gaps:', len(gaps), 'numeric defects:', len(numbugs), 'conflicts:', len(conflicts))


def _join(body):
    recs_idx = [i for i, b in enumerate(body) if not b.startswith('(* ----')]
    last = recs_idx[-1] if recs_idx else -1
    return '\n'.join(b + (',' if (not b.startswith('(* ----') and i != last) else '')
                     for i, b in enumerate(body))


if __name__ == '__main__':
    main()
