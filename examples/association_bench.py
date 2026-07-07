#!/usr/bin/env python3
"""Reference benchmark: Python dict / collections.Counter vs a hash-free scan.

Companion to examples/association_bench.m (Mathilda) and the write-up in
examples/association-benchmarks.md. Builds a frequency table of N random
integers three ways and reports wall-clock time:

  * collections.Counter   — CPython's C-level hash dict  (O(n))
  * plain dict loop        — Python-level hash dict        (O(n))
  * naive list scan        — no hashing, linear membership  (O(n*k))

The naive scan mirrors what a CAS *without* a real association type is forced
to do (rescan the key list on every element), and shows the quadratic blow-up
the hash-backed structure avoids.
"""

import random
import time
from collections import Counter


def timed(fn, repeats=5):
    best = float("inf")
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best * 1000.0  # milliseconds


def counter_counts(data):
    return Counter(data)


def dict_counts(data):
    d = {}
    for x in data:
        d[x] = d.get(x, 0) + 1
    return d


def naive_list_counts(data):
    keys, vals = [], []
    for x in data:
        for i, k in enumerate(keys):      # linear scan — no hashing
            if k == x:
                vals[i] += 1
                break
        else:
            keys.append(x)
            vals.append(1)
    return dict(zip(keys, vals))


def main():
    print("== Hash-backed (O(n)) ==")
    for n in (20000, 40000, 80000, 160000):
        random.seed(1)
        data = [random.randint(1, 1000) for _ in range(n)]
        print(f"  Counter   N={n:>7}: {timed(lambda: counter_counts(data)):8.4f} ms")
        print(f"  dict-loop N={n:>7}: {timed(lambda: dict_counts(data)):8.4f} ms")

    print("== Naive list scan (O(n*k), no hashing) ==")
    for n in (1000, 2000, 4000, 8000):
        random.seed(1)
        data = [random.randint(1, 200) for _ in range(n)]
        print(f"  listscan  N={n:>7}: {timed(lambda: naive_list_counts(data), repeats=1):8.4f} ms")


if __name__ == "__main__":
    main()
