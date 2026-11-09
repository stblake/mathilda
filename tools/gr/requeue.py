#!/usr/bin/env python3
"""
requeue.py -- drop records from a verify.py result cache so the next
`verify.py verify` run re-checks them (e.g. after a Mathilda fix).

Usage:
    python3 -I tools/gr/requeue.py <verified.json.cache.jsonl> STATUS[,STATUS...] [--novalue]

Removes every cached record whose oracle status is in the list.  --novalue
also requeues MISMATCH records whose evidence is a value Mathilda could not
numericize (a numeric gap, not a disagreement).
"""
import json
import sys

path, statuses = sys.argv[1], set(sys.argv[2].split(','))
novalue = '--novalue' in sys.argv
keep, dropped = [], 0
for line in open(path):
    try:
        r = json.loads(line)
    except ValueError:
        continue
    st = r.get('status')
    if st in statuses or (novalue and st == 'MISMATCH' and 'novalue' in (r.get('why') or '')):
        dropped += 1
        continue
    keep.append(line if line.endswith('\n') else line + '\n')
open(path, 'w').writelines(keep)
print('requeued %d, kept %d' % (dropped, len(keep)))
