#!/usr/bin/env python3
"""
merge.py -- gather the G&R transcription batches into one entries file.

Usage:
    python3 -I tools/gr/merge.py <raw-dir> <overrides.json> <entries.json>

<raw-dir> holds the per-agent arrays written under tools/gr/TRANSCRIBE.md
(tools/gr/raw/p<first>-<last>.json).  <overrides.json> maps an id to a dict of
field replacements (hand/agent repairs after the oracle rejected a reading);
a missing file is treated as empty.

Output: one JSON array in the schema tools/gr/verify.py reads
(id, status, lhs, rhs, conds) plus the bookkeeping fields the G&R builder and
triage use (gr, page, section, cite, tex, reason, key).  `key` is a sortable
rendering of the G&R number, so records sort in book order.  Duplicate ids
(two batches transcribing the same boundary item) keep the first occurrence
in page order and are reported.
"""
import glob
import json
import os
import re
import sys


def gr_key(gr):
    """'3.389.2a' -> '03.38900.002a' (book order; the 3-4 digit formula number
    is left-aligned so 3.38 < 3.381 < 3.389 < 3.39)."""
    m = re.match(r'(\d+)\.(\d+)(?:\.(\d+))?([a-z]*)', gr or '')
    if not m:
        return 'zz' + (gr or '')
    ch, num, item, suf = m.groups()
    return '%02d.%-5s.%03d%s' % (int(ch), num, int(item or 0), suf)


def main():
    raw_dir, ovr_path, out_path = sys.argv[1:4]
    overrides = json.load(open(ovr_path)) if os.path.exists(ovr_path) else {}
    files = sorted(glob.glob(os.path.join(raw_dir, '*.json')),
                   key=lambda p: int(re.search(r'p(\d+)', os.path.basename(p)).group(1)))
    seen = {}
    out = []
    dups = []
    for path in files:
        for r in json.load(open(path)):
            rid = r.get('id')
            if not rid:
                print('record without id in', path, file=sys.stderr)
                continue
            if rid in seen:
                dups.append((rid, seen[rid], os.path.basename(path)))
                continue
            seen[rid] = os.path.basename(path)
            if rid in overrides:
                r = dict(r, **overrides[rid], override=True)
            r.setdefault('rhs', [])
            r.setdefault('conds', [])
            r.setdefault('reason', r.get('note', ''))
            r['key'] = gr_key(r.get('gr'))
            r['batch'] = os.path.basename(path)
            out.append(r)
    out.sort(key=lambda r: (r['key'], r['id']))
    json.dump(out, open(out_path, 'w'), indent=1)
    st = {}
    for r in out:
        st[r['status']] = st.get(r['status'], 0) + 1
    print('records:', len(out), st)
    for d in dups:
        print('duplicate id %s in %s and %s (kept first)' % d)


if __name__ == '__main__':
    main()
