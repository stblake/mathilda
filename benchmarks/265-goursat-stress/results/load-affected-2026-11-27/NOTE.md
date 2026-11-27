# Load-affected partial run (2026-11-27) -- NOT a baseline

The first full timed run (30 s per case, 8 jobs), stopped after ~25 min. The machine's
load average reached 46 on 16 CPUs (other agents building and testing), so its times
and timeouts cannot be interpreted. It is kept only as a record. Do not compare it with
anything, and do not move it back to results/cache.jsonl.

Contents: plain sympy integrate on ~600 cases (170 answers, 131 unevaluated, 299
timeouts at 30 s), plus a 6-case smoke test of the Mathilda methods and of goursat.py.
