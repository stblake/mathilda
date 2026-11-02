# Association in Mathilda: audit

Mathilda v0.187 was compared against Mathematica 15.0.0 and Python 3.14 on an Apple M5 Max, 2026-09-25. Five independent agents each covered one angle; the detailed sections are linked below. I reproduced the headline claims by hand before writing this summary.

## Verdict

Reading from associations works well. Construction, lookup, grouping and merging match Mathematica 15, and a single-key lookup is a true O(1) hash probe. The problems are on the write side and in structural semantics:

- Adding keys one at a time is O(n²).
- Several ways of changing an association silently do nothing.
- Associations are not treated as indivisible, so `/.` can rewrite keys and merge two into one.

| | Result |
|---|---|
| Matches Mathematica 15 | 211 of 271 test expressions (77.9%) |
| The 60 mismatches | 32 wrong value, 16 unevaluated, 5 evaluated where Mathematica isn't, 7 formatting only |
| Wolfram functions | 81 checked: 31 full, 30 partial, 18 missing, 2 not applicable |
| Speed at 10⁶ entries | faster than Mathematica on 11 of 18 operations (`Keys` 16×, `KeyDrop` 5×, `Merge` 2.4×) |
| Worst slowdowns | update loops up to 1371× slower at 10⁴, O(n²); `PositionIndex` 18× slower at 10⁶ |
| Tests | 236 association and 19 compile tests pass, 5 of 5 benchmark gates pass, 0 leaks |
| Memory | about 350 bytes per entry, 2.8× Mathematica and 4× a Python `dict` |

## What to fix first

### P0: silent wrong answers

1. **Treat Association as indivisible in structural operations.** This one cause produces 11 of the 32 wrong values.
   - `<|a->1, b->2|> /. b->a` gives `<|a->2|>`, so a key is lost.
   - `AtomQ`, `Depth`, `Level`, `FreeQ` and `OrderedQ` all look inside the internal rules.
   - `Map` and `Apply` with a level spec build malformed `Association[f[a->1]]` nodes, and `AssociationQ` still returns `True` for them, because it only checks the head (`assoc.h:32`).
2. **Finish in-place changes and deletion.**
   - `c[k] += 1` and `c[k]++` fail with an `rvalue` message and leave the value unchanged.
   - `a[k] =.`, `Delete[a, Key[k]]` and `a[k1][k2] = v` are silent no-ops.
   - `KeyDropFrom` is interned but never registered.
   - As a result, the counting idiom `If[KeyExistsQ[c,#], c[#]+=1, c[#]=1]` returns wrong counts.
3. **Parse `#name` as `Slot["name"]`.** Today it parses as `Slot[1]*name`, so `Select[rows, #a > 2 &]` returns `{}` with no error. `Slot["a"]` is also not applied to associations.
4. **`Minus` doesn't evaluate on numbers.** This is a core bug, not an Association one. `Minus[3]` stays unevaluated, so `SortBy[{3,1,2}, Minus]` gives `{1,2,3}`.
5. **Part with spans, lists or `All` should return a smaller association.**
   - `a[[2;;3]]` and `a[[All]]` give `Missing["KeyAbsent", …]`, and `a[[{1,3}]]` gives a bare list of values.
   - Related wrong answers:
     - `Sort[a, p]` ignores `p`.
     - `KeyTake` ignores the requested order.
     - `Prepend` keeps the old value of an existing key.
     - `Catenate` of associations returns an association.

### P1: performance

1. **Make updates O(1) by changing the association in place when only one reference to it exists.**
   - An `AssociateTo` or `a[k] = v` loop takes 5.5 s for 10⁴ inserts, against 4 ms in Mathematica.
   - Every update copies all the entries, finds the key by linear scan (`part.c:283`) and rebuilds the canonical form (`assoc.c:1057–1089`).
   - Fix: add a fast path for the single-reference case, and let `AssocIndex` grow and update incrementally.
   - Add an update-loop gate to `tests/bench_assoc.c`, which currently can't catch this.
2. **Remove the hidden second canonicalisation pass.** `builtin_association` rebuilds a node that is already canonical, then throws the copy away (`assoc.c:219–226`). Every function that returns an association pays for it.
3. **Cut memory per entry.** Each rule allocates a fresh `Rule` symbol node (`assoc.c:90–93`).
4. **Add packed paths for `PositionIndex` (18× slower than Mathematica) and `GroupBy` (1.7× slower).** Model them on the machine-word `Counts`, which is 23× faster than Python's `Counter`.

### P2: coverage and cleanup

1. **Curried operator forms.** None of them evaluate: `Lookup[k]`, `Select[p]`, `KeyTake`/`KeyDrop`/`KeySelect`/`KeyMap`, `GroupBy[f]`, `CountsBy[f]`, `KeyValueMap[f]`, `AssociationMap[f]`, `Merge[f]` and `Apply[f]`.
2. **Arithmetic and equality.** `assoc + 1`, `Sqrt[assoc]` and `Mean[{a1, a2}]` stay symbolic, and so does `==` between unequal associations.
3. **Missing functions.** `KeyDropFrom`, `KeyIntersection`, `KeyComplement`, `MissingQ`, `JoinAcross`, `ApplyTo`, `Transpose` of associations, RawJSON import and export, `Query` and `Dataset`.
4. **Visible `NDArray` arguments.** `AssociationThread`, `AssociationMap` and `GroupBy` leave them unevaluated, which the project rules count as a wrong answer.
5. **Cleanup.**
   - `RuleDelayed` values are lost when a rebuild happens, for example on a duplicate key or a splice.
   - The index and the fallback scan disagree on malformed nodes: `Lookup[Association[f[1,2]], 1]` gives `2`.
   - Two comments contradict how the index is actually built and kept (`expr.h:205`, `assoc.c:291`).
   - `integrate_goursat.c:16` still says Mathilda lacks Association.
   - The `Counts` entry in the `EXEMPT` list of `compile_coverage.py` is stale.
   - `Key` and `Missing` have no docstring and no refpage.

## What already works

- **Construction.** Duplicate keys resolve as in Mathematica: the last value wins and the key keeps its first position.
- **Keys and lookup.** `1`, `1.` and `"1"` are distinct keys; `Missing["KeyAbsent", k]` is returned for absent keys; `Lookup` handles defaults and key lists.
- **Copying and assignment.** Copying works as expected: after `b = a`, changing `a` leaves `b` alone. Plain `a[k] = v` assignment also works.
- **Functions that match Mathematica.** `Merge`, every tested form of `GroupBy`, `Counts`/`CountsBy`, `AssociationThread`, `KeyValuePattern`, `SameQ`, `Hash` and `FullForm`.
- **Other subsystems.** `FindGraphIsomorphism`, `MapThread`, `Union` and `GroupBy` produce well-formed associations with O(1) lookup.
- **Docs.** All 27 Association builtins have refpages and docstrings, and all 49 doc examples reproduce.

## Sections

1. [Architecture & implementation](1-architecture.md): how an association is represented, the hash index, mutation, and complexity, with file:line citations.
2. [Feature coverage vs Mathematica 15](2-coverage.md): the 81-row inventory of Wolfram functions, checked with 333 test expressions.
3. [Correctness vs Mathematica 15](3-correctness.md): the 271-expression test battery, every mismatch, and minimal reproducers.
4. [Performance vs Mathematica 15 and Python](4-performance.md): the full timing table, scaling, memory, and root causes.
5. [Integration, tests & documentation](5-integration.md): packed arrays, `Compile`, the pattern matcher, tests, leaks, and docs.

Scripts and raw outputs are in `/tmp/assocreport/`:

- `diff/`: the correctness battery.
- `perf-bench/`: the benchmarks.
- `p*.m` and `leak.m`: the integration checks.
