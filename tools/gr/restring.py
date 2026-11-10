#!/usr/bin/env python3
"""
restring.py -- repair the printed strings of oracle records (no integration).

Usage:
    python3 -I tools/gr/restring.py <entries.json> <verified.json> <out.json> \
        [--mathilda ./Mathilda] [--jobs 8]

Before v0.6488, InputForm printed a held -1/q^(n + 1) as -q^(-n + 1).  The
oracle (tools/gr/verify.m) evaluates the held expressions, so its verdicts are
sound, but the integrand / limits / value STRINGS it wrote into
<verified.json> could be a different expression.  Those strings feed the
corpus (build_corpus.py) and the GR table (build_table.py).

This re-prints every record with the current printer (tools/gr/restring.m)
and writes <out.json>: the same records, with integrand, limits and value
replaced by the re-printed text.  The value is the rhs alternative whose
re-print equals the stored string.  When no alternative matches (the stored
string was corrupted), it is the single alternative, or with several
alternatives the record is marked "restring": "ambiguous" and left out of the
corpus by status RESTRING-AMBIGUOUS.  The record's "restring" field says
"same" / "fixed" / "ambiguous".
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify import mstr, mlist  # noqa: E402  (Mathilda string escaping)


def run_chunk(mathilda, calls):
    fd, path = tempfile.mkstemp(suffix='.m')
    with os.fdopen(fd, 'w') as fh:
        fh.write('Get["tools/gr/restring.m"];\n')
        fh.write('\n'.join(calls) + '\n')
    try:
        p = subprocess.run([mathilda, '-file', path], capture_output=True, text=True, timeout=1800)
        out = p.stdout
    finally:
        os.unlink(path)
    res = {}
    for line in out.splitlines():
        if line.startswith('@@S\t'):
            f = line.split('\t')
            if len(f) >= 7:
                _, rid, var, integ, lim, k, val = f[:7]
                res.setdefault(rid, []).append({'var': var, 'integrand': integ, 'limits': lim,
                                                'k': int(k), 'value': val})
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('entries')
    ap.add_argument('verified')
    ap.add_argument('out')
    ap.add_argument('--mathilda', default='./Mathilda')
    ap.add_argument('--jobs', type=int, default=8)
    a = ap.parse_args()
    ents = {e['id']: e for e in json.load(open(a.entries))}
    recs = json.load(open(a.verified))
    calls = []
    for r in recs:
        e = ents.get(r['id'])
        if e is None or not r.get('integrand'):
            continue
        try:
            variant = int(r.get('variant') or 1)
        except ValueError:
            variant = 1
        calls.append('grRestring[%s, %s, %s, %d];' % (mstr(r['id']), mstr(e['lhs']),
                                                      mlist(e.get('rhs', [])), variant))
    chunks = [calls[i:i + 200] for i in range(0, len(calls), 200)]
    got = {}
    with ThreadPoolExecutor(a.jobs) as ex:
        for res in ex.map(lambda c: run_chunk(a.mathilda, c), chunks):
            got.update(res)
    stats = {'same': 0, 'fixed': 0, 'ambiguous': 0, 'unprinted': 0}
    for r in recs:
        alts = got.get(r['id'])
        if not r.get('integrand'):
            continue
        if not alts:
            stats['unprinted'] += 1
            continue
        match = [x for x in alts if x['value'] == r['value']]
        pick = match[0] if match else (alts[0] if len(alts) == 1 else None)
        changed = pick is not None and (pick['integrand'] != r['integrand'] or
                                        pick['limits'] != (r.get('limits') or '') or
                                        pick['value'] != r['value'])
        if pick is None:
            r['restring'] = 'ambiguous'
            stats['ambiguous'] += 1
            continue
        if changed:
            r['restring_old'] = {'integrand': r['integrand'], 'limits': r.get('limits'),
                                 'value': r['value']}
            r['integrand'], r['limits'], r['value'] = pick['integrand'], pick['limits'], pick['value']
            r['restring'] = 'fixed'
            stats['fixed'] += 1
        else:
            r['restring'] = 'same'
            stats['same'] += 1
    json.dump(recs, open(a.out, 'w'), indent=1)
    print('restring:', stats)


if __name__ == '__main__':
    main()
