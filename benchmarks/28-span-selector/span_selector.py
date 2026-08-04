#!/usr/bin/env python3
"""Experiment 28 -- The start/step/n selector (numpy column).

Same six kernels as ``span_selector.m``, same order and sizes.

ROADMAP ITEM 6.  A strided span is describable in three integers; materialising a
position array first costs n boxed integers before a single element is read.  The
sieve is the extreme case -- nothing but strided writes -- and experiment 9
measured it 25x behind numpy while 1.18x AHEAD of Mathematica.  The two baselines
disagreeing is what identified this as a machine-level gap rather than a
competitive one, and is why both columns are kept.

``.copy()`` on every slice: numpy strided slices are views, and the Wolfram spans
materialise.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from harness import bench, check, require, seed

require(["numpy:arange"])

seed()
n = 4000000
v = np.random.random(n)

bench("span v[[1;;n;;2]]", lambda: v[0::2].copy())
check("span v[[1;;n;;2]]", int(np.arange(1, 11)[0::2].sum()))

bench("span v[[1;;n;;7]]", lambda: v[0::7].copy())
check("span v[[1;;n;;7]]", int(np.arange(1, 21)[0::7].sum()))

bench("span v[[2;;n;;2]] (offset)", lambda: v[1::2].copy())
check("span v[[2;;n;;2]] (offset)", int(np.arange(1, 11)[1::2].sum()))

bench("contiguous span v[[1;;n/2]]", lambda: v[:2000000].copy())
check("contiguous span v[[1;;n/2]]", int(np.arange(1, 11)[:5].sum()))

bench("Range[1, 4x10^6, 3]", lambda: np.arange(1, 4000001, 3))
check("Range[1, 4x10^6, 3]", int(np.arange(1, 21, 3).sum()))


def sieve(lim):
    """Count primes <= lim. Strided writes only -- the point of the row."""
    s = np.ones(lim + 1, dtype=np.int8)
    s[0:2] = 0
    i = 2
    while i * i <= lim:
        if s[i]:
            s[i * i::i] = 0
        i += 1
    return int(s.sum())


# The .m sieve counts over 1..lim with s[[1]] included then subtracts 1, which is
# the count of primes <= lim. Matched here by zeroing index 0 and 1.
bench("sieve to 200000 (strided writes)", lambda: sieve(200000), reps=1)
check("sieve to 200000 (strided writes)", sieve(100))
