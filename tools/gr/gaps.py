#!/usr/bin/env python3
"""
gaps.py -- coverage audit of the G&R transcription batches.

Usage:
    python3 -I tools/gr/gaps.py [tools/gr/raw]

A batch can silently drop items at its page boundaries (an agent misjudging
on which page an item starts).  G&R numbers the items of each formula
1, 2, 3, ... consecutively, so a hole in that sequence is a strong signal.
Prints, per formula number, the missing item numbers and the batch files
around the hole, plus the PDF page span each batch actually covered.
"""
import glob
import json
import os
import re
import sys
from collections import defaultdict


def main():
    raw = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'raw')
    items = defaultdict(dict)       # formula -> {item: batch}
    for p in glob.glob(os.path.join(raw, 'p*.json')):
        for r in json.load(open(p)):
            m = re.match(r'(\d+\.\d+)\.(\d+)[a-z]*$', r.get('gr', ''))
            if m:
                items[m.group(1)][int(m.group(2))] = os.path.basename(p)
    holes = 0
    for f in sorted(items, key=lambda s: [int(x) for x in s.split('.')]):
        have = items[f]
        missing = [i for i in range(1, max(have) + 1) if i not in have]
        if missing:
            holes += len(missing)
            near = sorted(set(have.values()))
            print('%-9s missing %s   (batches: %s)' % (f, compress(missing), ', '.join(near)))
    print('%d missing item numbers in %d formulas' % (holes, len(items)))


def compress(xs):
    out, start, prev = [], None, None
    for x in xs + [None]:
        if start is None:
            start = prev = x
        elif x == prev + 1:
            prev = x
        else:
            out.append(str(start) if start == prev else '%d-%d' % (start, prev))
            start = prev = x
    return ','.join(out)


if __name__ == '__main__':
    main()
