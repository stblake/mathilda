#!/usr/bin/env python3
"""
condaudit.py -- audit the book's conditions of the definite GR table records.

Usage:
    python3 -I tools/gr/condaudit.py <table-dir> <out.tsv> [--mathilda ./Mathilda] [--jobs 8] [--ids id,...]

Runs tools/gr/condaudit.m over every definite record of the generated table
files (tools/gr/build_table.py output): at up to four further parameter points
that satisfy the book's conditions, the book's value is compared with
NIntegrate.  Writes <out.tsv>, one line per record:

    id  agree  disagree  nint-failed  value-nonfinite  points  first-disagreement
        fix  fix-agree  fix-disagree

and, for tools/gr/build_table.py and tools/gr/build_corpus.py,
<out.tsv>.fix (id<TAB>the derived extra condition, --cond-fix) and
<out.tsv>.wrong (the uncorrected records, id<TAB>reason, --book-wrong): a
formula whose stated conditions admit a point where it is wrong is corrected
from first principles or left out.
"""
import argparse
import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor

FILES = ['DefHalfLine.m', 'DefWholeLine.m', 'DefUnit.m', 'DefTrig.m', 'DefOther.m']


def run_part(mathilda, tdir, k, n, ids=None):
    fd, path = tempfile.mkstemp(suffix='.m')
    with os.fdopen(fd, 'w') as fh:
        fh.write('gAuditK = %d; gAuditN = %d; gAuditIds = %s;\n' % (
            k, n, 'All' if ids is None else '{%s}' % ', '.join('"%s"' % i for i in ids[k::n])))
        fh.write('Get["tools/gr/condaudit.m"];\n')
        for f in FILES:
            fh.write('Get["%s"];\n' % os.path.join(os.path.abspath(tdir), f))
    try:
        p = subprocess.run([mathilda, '-file', path], capture_output=True, text=True, timeout=36000)
        return [l for l in p.stdout.splitlines() if l.startswith('@@A\t')]
    finally:
        os.unlink(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('tdir')
    ap.add_argument('out')
    ap.add_argument('--mathilda', default='./Mathilda')
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--ids', default='', help='comma-separated record ids (default: all)')
    a = ap.parse_args()
    with ThreadPoolExecutor(a.jobs) as ex:
        ids = [i for i in a.ids.split(',') if i] or None
        parts = list(ex.map(lambda k: run_part(a.mathilda, a.tdir, k, a.jobs, ids), range(a.jobs)))
    rows = sorted(l.split('\t')[1:] for part in parts for l in part)
    with open(a.out, 'w') as fh:
        for r in rows:
            fh.write('\t'.join(r) + '\n')
    bad = [r for r in rows if int(r[2]) > 0]
    fixed = [r for r in bad if len(r) > 9 and r[7] and int(r[8]) > 0 and int(r[9]) == 0]
    print('condaudit: %d records, %d with a disagreement (%d corrected), %d with no checkable point' % (
        len(rows), len(bad), len(fixed), sum(1 for r in rows if int(r[1]) == 0 and int(r[2]) == 0)))
    # the two outputs build_table.py / build_corpus.py read: corrections
    # (id<TAB>extra condition) and removals (id<TAB>reason)
    with open(a.out + '.fix', 'w') as fh:
        for r in fixed:
            fh.write('%s\t%s\n' % (r[0], r[7]))
    with open(a.out + '.wrong', 'w') as fh:
        for r in bad:
            if r not in fixed:
                fh.write('%s\tcondition audit: wrong at a point its conditions admit %s\n' % (r[0], r[6]))


if __name__ == '__main__':
    main()
