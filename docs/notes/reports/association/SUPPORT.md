# Association support in Mathilda (as of v0.200)

This is the current state after the `association-overhaul` branch. The audit it responds to is in [README.md](README.md), with detailed sections 1–5. Figures come from the 271-case differential battery against Mathematica 15 (`/tmp/assocreport/diff/battery.m`) and from the test suites.

## At a glance

| | Before (v0.196) | Now (v0.200) |
|---|---|---|
| Battery cases matching Mathematica 15 | 211 / 271 (77.9%) | **250 / 271 (92.3%)**, no regressions |
| Association test assertions | 236 + 19 compile | **+691** in 5 new suites |
| Operator (curried) forms | none evaluated | **43 heads** |
| Full ctest | 52 known failures on main | same 52; `moebiusmu`/`primenu` flaky under load, pass alone |

## Supported, matching Mathematica 15

- **Construction and identity.**
  - `<|…|>`, `Association[{rules}]`, nested rule lists, `AssociationThread`, `AssociationMap` (lists and associations).
  - Duplicate keys: the last value wins and the key keeps its first position. `RuleDelayed` values are kept.
  - Keys are compared with `SameQ` (`1`, `1.` and `"1"` are distinct).
  - `AssociationQ` accepts only well-formed associations.
- **Reading.**
  - `a[k]`, `Lookup` (with a lazy default, key lists, lists of associations and operator forms), `KeyExistsQ`, and `KeyMemberQ`/`KeyFreeQ` with patterns.
  - `Missing["KeyAbsent", k]` for absent keys, and `MissingQ`.
  - Part: `a[[k]]`, `a[[Key[k]]]`, `a[[i]]`, and spans, lists and `All` returning sub-associations.
  - `Keys`/`Values`, including the function form and threading over lists.
- **Atomicity.**
  - `AtomQ`, `Depth`, `LeafCount` and `Level` treat an association as atomic; so do the level specs of `Map`/`Apply`/`Scan`/`MapIndexed`/`Replace`.
  - `FreeQ`, `Cases`, `Count`, `MemberQ`, `DeleteCases`, `Position` (returns `Key[k]`) and `OrderedQ` work on values.
  - `/.` never rewrites keys, and an association can be used as a rule set.
- **Transforming.**
  - `Map`, `KeyMap`, `KeyValueMap`, `Select`/`Discard`, `KeySelect`, `KeyTake` (in the requested order), `KeyDrop`.
  - `Sort`/`SortBy`/`ReverseSort`/`Ordering`/`KeySort`/`KeySortBy`, with ordering functions and Mathematica's tie placement.
  - `Take`/`Drop`/`First`/`Last`/`Reverse`/`Append`/`Prepend`/`Insert`/`Pick`/`Catenate`.
- **Combining.**
  - `Join`; `Merge` of associations and of rule lists; `KeyUnion` with fill; `KeyIntersection`, `KeyComplement`.
  - `JoinAcross`, with all four join types and `KeyCollisionFunction`.
  - `Transpose` of an association of associations.
- **Aggregating.**
  - `Counts`, `CountsBy`, `CountDistinct(By)`, `GroupBy` (multi-level, with a reducer), `PositionIndex` (lists and associations).
  - `Total`/`Mean` on single associations and on lists of associations.
- **Arithmetic.** Listable heads thread over the values: `<|a->1|> + 1`, `Sqrt[a]`, `a1 + a2` when the keys match.
- **Equality.** `==` and `!=` between associations; `SameQ`; `Hash`.
- **Patterns.** `_Association`, `KeyValuePattern`, literal association patterns, DownValues on `f[a_Association]`.
- **Pure functions.** Named slots: `#name` and `#"key"` read keys, so `Select[rows, #a > 2 &]` works.
- **Other heads.** `ApplyTo`, `AssociationComap`, `Splice`, `SubsetQ`, deep `Normal`, `DeleteMissing` with levels, `Key`.
- **Assignment.** Plain `a[k] = v`, `a[[k]] = v`, `a[k1, k2] = v`, `AssociateTo`, `AppendTo[a[k], …]`. Copy semantics are safe: after `b = a`, changing `a` leaves `b` unchanged.
- **Compile.** `Lookup`, `KeyExistsQ`, `Length`, `Values`, `KeyDrop`/`KeyTake` and `Counts`.

## Deliberate differences (documented in `docs/spec/builtins/data-structures.md`)

- Values rebuilt by `/.` are re-evaluated: `<|"b"->x^2|> /. x->3` gives `9`, where Mathematica gives `3^2`.
- Associations are atomic under `Hold` as well.
- Unevaluated `Part` prints as `Part[a, 5]`, not `a[[5]]`, and there are no `Part::partw` messages.
- A few extensions that Mathematica rejects: `Accumulate`, `Differences` and `Fold` on associations.

## Not yet supported (next batch)

These were in progress on separate worktrees and were stopped before merging.

| Area | Gap | Impact |
|---|---|---|
| Mutation | `a[k] += v`, `a[k]++`, `a[k] =.`, `KeyDropFrom`, `Delete[a, Key[k]]`, nested `a[k1][k2] = v` | Silent no-ops or `rvalue` errors, so counter idioms give wrong results |
| Write performance | `a[k] = v` / `AssociateTo` loops are O(n²): 10⁴ inserts take about 5.5 s, against 4 ms in Mathematica | The most common way to build an association incrementally |
| Memory | About 350 B per entry, 2.8× Mathematica | Large associations |
| Packed / NDArray | `Values` never packed; `PositionIndex` 18× slower than Mathematica; `GroupBy` 1.7× slower; `NDArray` arguments to `AssociationThread`/`AssociationMap`/`GroupBy` stay unevaluated | Numeric workloads |
| Compile | `Keys`, `p[k]` and `Part` do not lower; `Map`/`Select`/`Append` do not compose | Compiled code |
| Tabular and IO | `Query`, `Dataset`, `ImportString`/`ExportString` `"RawJSON"` | Data pipelines |
| Minor | `Catenate[<|x->1|>]`; `Equal` on lists inside values; `ReplacePart` into nested associations | Edge cases |

Worktree branches with the partial work: `worktree-agent-a9b92fd7413dec2b0` (mutation and performance), `worktree-agent-a34410e89db508f5e` (packed and Compile), and `worktree-agent-aaef5c17865d52b27` (Query, Dataset, JSON). None of them was finished or verified.
