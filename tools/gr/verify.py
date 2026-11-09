#!/usr/bin/env python3
"""
verify.py -- drive the Mathilda-side G&R tools (tools/gr/verify.m, check.m).

Usage:
    python3 -I tools/gr/verify.py check  <batch.json>
    python3 -I tools/gr/verify.py verify <entries.json> <verified.json> \
        [--ids id,id,...] [--jobs 8] [--chunk 10] [--mathilda ./Mathilda]

check   the transcription self-check (tools/gr/TRANSCRIBE.md): every record of
        one batch through check.m; prints errors/warnings, exit 1 on an error.
verify  every status-"ok" record of a merged entries file (tools/gr/merge.py)
        through the numeric oracle verify.m; writes the result list that
        tools/gr/build_corpus.py reads.

Mathilda runs .m files: this script only writes `grVerify[...]` / `grCheck[...]`
calls (strings escaped for Mathilda's lexer) into a temporary .m after a
Get of the tool file, runs `./Mathilda -file` on it, and parses the tagged
lines the tool prints.  Chunks run in parallel, one process each; a chunk
whose process dies (crash, hard timeout) is re-run one record per process so
a single crash costs a single record (reported as CRASH).
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))


def mstr(s):
    """Python str -> Mathilda string literal."""
    return '"' + s.replace('\\', '\\\\').replace('"', '\\"').replace('\n', ' ') + '"'


def mlist(xs):
    return '{' + ', '.join(mstr(x) for x in xs) + '}'


def run_mathilda(mathilda, tool, calls, hard_timeout):
    fd, path = tempfile.mkstemp(suffix='.m')
    with os.fdopen(fd, 'w') as fh:
        fh.write('Get[%s];\n' % mstr(os.path.join(HERE, tool)))
        fh.write('\n'.join(calls) + '\n')
    try:
        p = subprocess.run([mathilda, '-file', path], capture_output=True, text=True,
                           timeout=hard_timeout)
        return p.stdout, p.returncode == 0
    except subprocess.TimeoutExpired as ex:
        out = ex.stdout or ''
        return (out.decode(errors='replace') if isinstance(out, bytes) else out), False
    finally:
        os.unlink(path)


def preflight(mathilda, tool):
    """Load the tool alone; a syntax error in it would otherwise read as every
    record producing no output."""
    fd, path = tempfile.mkstemp(suffix='.m')
    with os.fdopen(fd, 'w') as fh:
        fh.write('Get[%s];\nPrint["@@LOADED"];\n' % mstr(os.path.join(HERE, tool)))
    try:
        p = subprocess.run([mathilda, '-file', path], capture_output=True, text=True, timeout=120)
    finally:
        os.unlink(path)
    text = p.stdout + p.stderr
    if '@@LOADED' not in p.stdout or 'sntx' in text or 'syntax error' in text:
        sys.exit('%s does not load cleanly in Mathilda:\n%s' % (tool, text.strip()[:2000]))


def cmd_check(a):
    preflight(a.mathilda, 'check.m')
    recs = json.load(open(a.batch))
    calls = ['grCheck[%s, %s, %s, %s, %s];' % (
        mstr(r.get('id', '?')), mstr(r.get('status', '')), mstr(r.get('lhs', '')),
        mlist(r.get('rhs', [])), mlist(r.get('conds', []))) for r in recs]
    out, ok = run_mathilda(a.mathilda, 'check.m', calls, 600)
    errs = warns = 0
    seen = {}
    for r in recs:
        seen[r.get('id')] = seen.get(r.get('id'), 0) + 1
        if not isinstance(r.get('gr'), str) or not isinstance(r.get('page'), int):
            print('ERROR %s: missing gr/page' % r.get('id')); errs += 1
    for k, v in seen.items():
        if v > 1:
            print('warn  %s: duplicate id' % k); warns += 1
    for line in out.splitlines():
        if line.startswith('@@E\t'):
            _, rid, msg = line.split('\t', 2); print('ERROR %s: %s' % (rid, msg)); errs += 1
        elif line.startswith('@@W\t'):
            _, rid, msg = line.split('\t', 2); print('warn  %s: %s' % (rid, msg)); warns += 1
    if not ok:
        print('ERROR: Mathilda did not finish the check (crash or timeout)'); errs += 1
    print('%d records, %d ok; %d errors, %d warnings' % (
        len(recs), sum(r.get('status') == 'ok' for r in recs), errs, warns))
    sys.exit(1 if errs else 0)


def parse_verify(out):
    res = {}
    for line in out.splitlines():
        if line.startswith('@@R\t'):
            f = line.split('\t')
            f += [''] * (11 - len(f))
            _, rid, st, kind, var, integrand, limits, value, params, variant, why = f[:11]
            res[rid] = {'id': rid, 'status': st, 'kind': kind, 'var': var,
                        'integrand': integrand, 'limits': limits, 'value': value,
                        'params': [p for p in params.split(',') if p],
                        'variant': variant, 'why': why, 'good': []}
        elif line.startswith('@@F\t'):
            f = line.split('\t')
            f += [''] * (10 - len(f))
            _, rid, inst, ff, fvar, flo, fhi, fval, mval, mjudge = f[:10]
            r = res.setdefault(rid, {'id': rid, 'good': []})
            r.setdefault('probes', []).append({'inst': inst, 'f': ff, 'var': fvar, 'lo': flo, 'hi': fhi,
                                               'value': fval, 'mathilda_N': mval, 'mathilda_judge': mjudge})
        elif line.startswith('@@G\t'):
            f = line.split('\t')
            f += [''] * (8 - len(f))
            _, rid, inst, ref, num, tier, ii, il = f[:8]
            g = {'inst': inst, 'ref': ref if ref != '""' else '', 'num': num if num != '""' else '',
                 'tier': tier}
            r = res.setdefault(rid, {'id': rid, 'good': []})
            r['good'].append(g)
            if 'inst_integrand' not in r:   # first verified instantiation = numeric instance
                r['inst_integrand'] = ii
                r['inst_limits'] = il
    return res


def cmd_verify(a):
    preflight(a.mathilda, 'verify.m')
    ents = [e for e in json.load(open(a.entries)) if e.get('status') == 'ok']
    if a.ids:
        want = set(a.ids.split(','))
        ents = [e for e in ents if e['id'] in want]
    call = {e['id']: 'grVerify[%s, %s, %s, %s];' % (
        mstr(e['id']), mstr(e['lhs']), mlist(e.get('rhs', [])), mlist(e.get('conds', [])))
        for e in ents}
    ids = [e['id'] for e in ents]
    # Resumable: every finished chunk is appended to <verified>.cache.jsonl
    # (keyed by id and by the record's lhs/rhs/conds text, so an edited record
    # is re-run); a restarted run skips what the cache already holds.
    cache_path = a.verified + '.cache.jsonl'
    sig = {e['id']: json.dumps([e['lhs'], e.get('rhs', []), e.get('conds', [])]) for e in ents}
    cached = {}
    if os.path.exists(cache_path) and not a.fresh:
        for line in open(cache_path):
            try:
                r = json.loads(line)
            except ValueError:
                continue
            if r.get('_sig') == sig.get(r['id']):
                cached[r['id']] = r
    todo = [i for i in ids if i not in cached]
    print('verify: %d records, %d cached, %d to run' % (len(ids), len(cached), len(todo)), flush=True)
    chunks = [todo[i:i + a.chunk] for i in range(0, len(todo), a.chunk)]
    import threading
    lock = threading.Lock()
    cache_fh = open(cache_path, 'a')

    def work(ch):
        # verify.m bounds each record at 300 s with TimeConstrained; the hard
        # process timeout is the backstop for work TimeConstrained cannot
        # interrupt, after which the chunk is re-run one record per process.
        out, ok = run_mathilda(a.mathilda, 'verify.m', [call[i] for i in ch], 320 * len(ch) + 60)
        got = parse_verify(out)
        if not ok:
            for i in ch:
                if i not in got or 'status' not in got[i]:
                    o1, ok1 = run_mathilda(a.mathilda, 'verify.m', [call[i]], 360)
                    g1 = parse_verify(o1)
                    got[i] = g1.get(i) if i in g1 and 'status' in g1[i] else {
                        'id': i, 'status': 'CRASH' if not ok1 else 'NOOUTPUT', 'good': []}
        with lock:
            for i in ch:
                if i in got:
                    cache_fh.write(json.dumps(dict(got[i], _sig=sig[i])) + '\n')
            cache_fh.flush()
        return got

    res = dict(cached)
    with ThreadPoolExecutor(a.jobs) as ex:
        for g in ex.map(work, chunks):
            res.update(g)
    cache_fh.close()
    out = [{k: v for k, v in res.get(i, {'id': i, 'status': 'NOOUTPUT', 'good': []}).items() if k != '_sig'}
           for i in ids]
    for r in out:   # an instantiation verified only "loose" marks the record;
        # "loose-silentwrong": the default NIntegrate returned a wrong number
        # without a message -- an NIntegrate bug, surfaced in the gaps list
        r['silentwrong'] = any(g.get('tier') == 'loose-silentwrong' for g in r.get('good', []))
        r['tier'] = 'loose' if any(g.get('tier', '').startswith('loose') for g in r.get('good', [])) else 'strict'
    json.dump(out, open(a.verified, 'w'), indent=1)
    cnt = {}
    for r in out:
        k = r['status'] + ('-loose' if r['status'].startswith('VERIFIED') and r['tier'] == 'loose' else '')
        cnt[k] = cnt.get(k, 0) + 1
    print('verified %d records:' % len(out), dict(sorted(cnt.items())))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    c = sub.add_parser('check'); c.add_argument('batch')
    v = sub.add_parser('verify'); v.add_argument('entries'); v.add_argument('verified')
    v.add_argument('--ids'); v.add_argument('--jobs', type=int, default=8)
    v.add_argument('--chunk', type=int, default=6)
    v.add_argument('--fresh', action='store_true', help='ignore the result cache')
    for p in (c, v):
        p.add_argument('--mathilda', default=os.path.join(HERE, '..', '..', 'Mathilda'))
    a = ap.parse_args()
    cmd_check(a) if a.cmd == 'check' else cmd_verify(a)


if __name__ == '__main__':
    main()
