#!/usr/bin/env python3
"""Experiment 25 -- Vector math for the transcendentals (numpy column).

Same eight kernels as ``vector_transcendentals.m``, same order and sizes.

ROADMAP ITEM 8.  Accelerate is already linked into Mathilda and ships
vvsin/vvexp/vvlog at roughly twice a scalar libm loop's throughput.  numpy calls
its own vectorised loops (and Accelerate where available), so a Mathilda gap of
~2x on these rows confirms the item and measures its size.

The ``Sin[Exp[Log[v]]]`` row is the fusion probe: one pass if the system fuses,
three passes plus two temporaries if not.  numpy does NOT fuse either -- it
materialises each intermediate -- so a large Mathilda gap there is about
per-element cost, not about fusion, and a large numpy gap would be the reverse.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

import numpy as np

from harness import bench, check, require, seed

require(["numpy:sin", "numpy:exp", "numpy:log", "numpy:sqrt", "numpy:tan",
         "numpy:arctan2", "numpy:power"])

seed()
n = 4000000
v = np.random.random(n) + 0.5


def r6(x):
    return int(math.floor(float(x) * 1e6 + 0.5))


bench("Sin over 4x10^6", lambda: np.sin(v))
check("Sin over 4x10^6", r6(math.sin(0.5)))

bench("Exp over 4x10^6", lambda: np.exp(v))
check("Exp over 4x10^6", r6(math.exp(0.5)))

bench("Log over 4x10^6", lambda: np.log(v))
check("Log over 4x10^6", r6(math.log(1.5)))

bench("Sqrt over 4x10^6", lambda: np.sqrt(v))
check("Sqrt over 4x10^6", r6(math.sqrt(1.5)))

bench("Tan over 4x10^6", lambda: np.tan(v))
check("Tan over 4x10^6", r6(math.tan(0.5)))

# Wolfram's ArcTan[x, y] == atan2(y, x): the argument order is reversed.
bench("ArcTan 2-arg over 4x10^6", lambda: np.arctan2(v + 1, v))
check("ArcTan 2-arg over 4x10^6", r6(math.atan2(1.5, 0.5)))

bench("Sin[Exp[Log[v]]] fused?", lambda: np.sin(np.exp(np.log(v))))
check("Sin[Exp[Log[v]]] fused?", r6(math.sin(math.exp(math.log(1.5)))))

bench("v^2.5 over 4x10^6", lambda: v ** 2.5)
check("v^2.5 over 4x10^6", r6(1.5 ** 2.5))
