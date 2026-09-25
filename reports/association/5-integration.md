## 5. Integration, tests & documentation

This section asks how Association fits with the rest of Mathilda: the packed-array/NDArray substrate, `Compile[]`, the pattern matcher and rule engine, and the other subsystems that return associations. It also covers the test suites and the reference documentation. Every claim below was checked against the shipped binary (`./mathilda`, built Sep 24 16:40) and the prebuilt test binaries in `build/`. The probe scripts are in `/tmp/assocreport/p*.m`. Nothing in the repo was edited or rebuilt.

Source footprint: `src/assoc.c` (1624 lines, 25 registered builtins), `src/assoc_index.c` (83 lines, the persistent key→position hash), `src/compile/compile_assoc.c` (621 lines), `tests/test_association.c` (1689 lines), `tests/test_compile_assoc.c` (493 lines), and `tests/bench_assoc.c`.

### 5.1 Packed arrays and NDArray

**Design position.** `src/pack.c` keeps Association off the `AWARE` list on purpose. It is part of the *no-nesting invariant*:

```
 * DELIBERATELY ABSENT, and each for a reason:
 *   List, Rule, RuleDelayed, Association, Hold, HoldForm
 *       These enforce the NO-NESTING INVARIANT. ... a packed node can never sit
 *       inside a plain EXPR_FUNCTION tree
```

The only association-producing head on `AWARE` is `Counts`, which "relabels Tally's machine-word count" (the structural batch that names `src/assoc.c`). `tools/nd_fastpath_sweep.py` records the rest as known materialisers in its `OFF_BUFFER` ratchet: `"AssociationMap", "AssociationThread", "KeyExistsQ", "KeyFreeQ", "KeyMemberQ"`, `"CountsBy", "GroupBy"`, `"AssociationQ"` and `"Lookup"`. None of `check_packed_aware.py`, `nd_surface_audit.py` or `check_array_exactness.py` has an Association-specific probe. `nd_surface_audit.py` explicitly treats "returns ... an Association" as correctly unpacked.

**What happens in practice** (`p1.m`, `p2.m`, `p3.m`). Note that `Developer`PackedArrayQ` and `Developer`ToPackedArray` are *not* defined. The context-free `PackedArrayQ` and `ToPackedArray` (an alias of `ToNDArray`) are:

```
$ ./mathilda /tmp/assocreport/p1.m      (a = AssociationThread[Range[1000], N[Range[1000]]])
PackedArrayQ[Values[a]] = False
PackedArrayQ[Keys[a]] = False
Dev PackedArrayQ Values: Developer`PackedArrayQ[{1.0, 2.0, ...}]      <- undefined, returned unevaluated
PackedArrayQ[v] = True                                                <- v = Range[1000]*1.0
PackedArrayQ[b["x"]] = False                                          <- b = <|"x" -> v|>
PackedArrayQ[Lookup[b,"x"]] = False
Total[Values[a]] = 500500.0  Mean = 500.5
Counts vals packed: False len 10 total 100000
Total[a] = 500500.0 Mean[a] = 500.5
Map packed: False
Keys packed after KeySort: False

$ ./mathilda /tmp/assocreport/p2.m      (10^6 entries)
vv=Values[a]; PackedArrayQ[vv] = False
time Total[Values[a]] 0.034502
time Total[packed] 0.000207               <- same data as a packed list: ~170x faster
time Values[a] 0.014016
time AssociationThread[1e6] 0.388316
time Counts packed 1e6 0.003869
NDArray AssociationThread: AssociationThread[{a, b, c}, NDArray[{1.0, 2.0, 3.0}]]   <- unevaluated
NDArray AssociationMap:    AssociationMap[#1^2 &, NDArray[{1.0, 2.0, 3.0}]]         <- unevaluated
Counts NDArray: <|1 -> 2, 2 -> 1|>
GroupBy NDArray: GroupBy[NDArray[{1.0, 2.0, 3.0}], #1 > 1 &]                        <- unevaluated

$ ./mathilda /tmp/assocreport/p3.m
PackedArrayQ[ToPackedArray[Values[a]]] = True
PackedArrayQ[Values[a]+0.] = False
```

The gate diagnostics confirm where materialisation happens (`MATHILDA_PACK_DIAG=gate`):

```
r = AssociationThread[Range[10^5], v];   ->  AssociationThread   2   200000
r = GroupBy[v, EvenQ];                   ->  GroupBy             1   100000
r = <|"x" -> v|>;                        ->  Rule                1   100000
r = Counts[v];                           ->  (no materialisation)
```

In summary:
- Values and Keys never come back packed, even at 10^6 machine reals.
- A packed list stored as a value is boxed at `Rule` construction.
- `Total` and `Mean` over `Values` give correct answers but run on boxed `Expr`s, about 170× slower than the same reduction on a packed buffer.
- Nothing re-packs the result automatically, not even arithmetic (`Values[a]+0.`).
- A *visible* `NDArray` passed to `AssociationThread`, `AssociationMap` or `GroupBy` is left unevaluated. CLAUDE.md calls this outcome a wrong answer, not merely a slow one. `Counts` is the one head that handles both surfaces.
- A packed or NDArray list used as a **key** hashes identically to its boxed twin. `tests/test_packed_list.c:573-574` covers this, and it reproduces by hand (`<|ToNDArray[{1., 2.}] -> "x"|>[{1., 2.}]` gives `x`).

### 5.2 Compile

`compile_assoc.c` lowers these heads (`try_emit_assoc`, line 259): `Lookup` (2 or 3 args), `KeyExistsQ`/`KeyMemberQ`, `KeyFreeQ`, `Length`, `Values`, `KeyDrop`, `KeyTake`, `Counts` (on a rank-1 machine array), `Map[f, assoc]`, `Select[assoc, pred]` and `Append[assoc, k -> v]`. They work over a declared `{p, _Association[, _Integer|_Real|_Complex]}` argument, a literal `<|...|>`, or a folded global. Keys must be compile-time constants, except for an `_Integer` runtime key into `Lookup` (B2). Only `KeyDrop`, `KeyTake` and `Counts` count as *producer* heads (`assoc_producer_head`, line 94), so only their results can feed another association op.

```
$ ./mathilda /tmp/assocreport/p4.m ; ./mathilda /tmp/assocreport/p5.m      (abridged)
{{p,_Association},{x,_Real}}  Lookup[p,"a"] + x                 -> Compiled -> True, ResultType -> Real
{{p,_Association,_Real}}      Total[Values[p]]                  -> Compiled -> True, ResultType -> Real
{{p,_Association,_Real}}      Total[Values[KeyDrop[p,"b"]]]     -> Compiled -> True
{{v,_Integer,1}}              Length[Counts[v]]                 -> Compiled -> True, ResultType -> Integer
{{v,_Integer,1}} / {{v,_Real,1}}  Counts[v]                     -> Compiled -> True, ResultType -> Unknown
{{p,_Association,_Real}}      Map[#^2 &, p]                     -> Compiled -> True, ResultType -> Unknown
{{p,_Association,_Real}}      Table[Lookup[p,"a"] i, {i,3}]     -> Compiled -> True, ResultType -> Array
{{p,_Association,_Real}}      Total[Values[Map[#^2 &, p]]]      -> Compiled -> False  (Subexpression -> Values[Map[#1^2 &, p]])
{{p,_Association,_Real}}      Length[Select[p, # > 1. &]]       -> Compiled -> False
{{p,_Association,_Real},{x,_Real}}  Lookup[Append[p,"z"->x],"z"] -> Compiled -> False
{{p,_Association,_Real}}      Keys[p] / p["a"] / p[["a"]] / KeySort[p] / Merge[{p,p},Total]  -> Compiled -> False
{{v,_Real,1}}   Total[Values[AssociationThread[Range[Length[v]], v]]]            -> Compiled -> False
f = Compile[{{p,_Association,_Real}}, Total[Map[#^2&,p]]]; f[<|"a"->1.,"b"->2.|>] = 5.0   (falls back to the interpreter)
g = Compile[{{p,_Association,_Real}}, Lookup[p,"q",-1.]];  g[<|"q" -> "s"|>] = s             (graceful type fallback)
h = Compile[{{p,_Association,_Real}}, Lookup[p,"a"]+1.];   h[{1.,2.}] = 1.0 + Missing[KeyAbsent, a]
interpreter:  Lookup[{1.,2.},"a",-1.] = -1.0
```

The subset is coherent for read-only parameter bags, and it degrades cleanly: a wrong value type falls back instead of crashing, and there are no leaks (5.5). Four gaps stand out:
- `Map`, `Select` and `Append` lower *at the top level* but are not producers, so they cannot be composed (`Values[Map[...]]`, `Length[Select[...]]` and `Lookup[Append[...]]` all fall off the cliff).
- `Keys`, the `p[k]` accessor and `Part[p, k]` do not lower at all, although these are the most idiomatic read forms.
- A non-association list passed where `_Association` was declared is accepted silently rather than rejected.
- `tools/compile_coverage.py:218` still lists `"Counts": "returns an Association, which has no compiled type"` in `EXEMPT`, and `COMPILE_MISSING.md` repeats that reason. Both are now stale, because `Counts[v]` compiles.

### 5.3 Pattern matcher and rules

```
$ ./mathilda /tmp/assocreport/p6.m ; ./mathilda /tmp/assocreport/p7.m      (a = <|"a"->1,"b"->2,"c"->3|>)
1 ReplaceAll values: <|a -> 1, b -> 20, c -> 3|>
2 ReplaceAll rule key: <|a -> 1, b -> 20, c -> 3|>
4 AssociationQ after replace: True Lookup: 20            <- rebuilt association is indexed and usable
5 key rewrite: <|z -> 1, b -> 2, c -> 3|>   5b Lookup: 1
7 MatchQ KVP: True True False
8 Cases on assoc (values): {1, 3}
9 Select: <|b -> 2|>
10 Cases KVP cond: {<|n -> x, v -> 5|>}
11 Cases KVP extract: {x5, y1}
12 downvalue _Association: 6 no
13 literal assoc pattern g[<|"a" -> x_, ___|>]: 1 7
14 Condition on assoc: big small
18 Position: {{Key[b]}}
19 FreeQ/MemberQ: False True False                         <- MemberQ looks at values, not keys (correct)
21 KVP on nested: True        22 ReplaceAll inside nested: <|p -> <|q -> 2|>|>
23 ReplaceRepeated: <|a -> 5|>   24 downvalue KVP: 42   25 KVP orderless match: True
--- divergences ---
17 a /. ("a" -> 1) -> 5 :  Association[5, b -> 2, c -> 3]
   AssociationQ: True  Length: 3  Keys: Keys[Association[5, b -> 2, c -> 3]]
Level[a,{1}]: {a -> 1, b -> 2, c -> 3}                     (WL: {1, 2, 3})
Replace[a, 1 -> 100, {1}]: <|a -> 1, b -> 2, c -> 3|>      (WL: <|a -> 100, ...|>)
Map[f, a, {1}]: Association[f[a -> 1], f[b -> 2], f[c -> 3]]   (Map[f, a] is fine: <|a -> f[1]|>)
Apply[g, <|"a" -> {1,2}|>, {1}]: Association[g[a, {1, 2}]]
{"a","b"} /. a : {a, b}       ReplaceAll[x, <|x -> 1|>]: x   (WL uses an association as a rule set: {1, 2})
```

The core matcher integration is solid. `_Association`, `KeyValuePattern` (with a condition, nested, orderless, or bound in a DownValue), `Condition`, `Cases`/`Select`/`Position`, and `ReplaceAll` into both keys and values all behave correctly. Two classes of defect remain:

1. **Level-spec traversal sees the internal `Rule` entries** rather than the values. `Level`, `Replace` at `{1}`, and `Map`/`Apply` with an explicit level are affected. The `Map` and `Apply` cases build malformed `Association[...]` nodes.
2. **`AssociationQ` does not check well-formedness.** `is_association` is a bare head test (`assoc.h:32`), so it returns `True` for those malformed nodes while `Keys` stays unevaluated on them.

An association also cannot be used as a replacement-rule set.

### 5.4 Other subsystems returning associations

Code outside `assoc.c` that builds `Association` nodes directly includes `src/graph/galg_isoheads.c:310` (`gi_assoc`), `src/funcprog.c:536`, `1055` and `3452` (Map/MapThread rebuilds), and `src/list/setops.c:956`. The key index is lazy: `assoc_index.h` says it is "built LAZILY, by the first single-key read", so every producer is indexed on first use as long as its keys are unique.

```
$ ./mathilda /tmp/assocreport/p8.m     (abridged; each row checks AssociationQ, Keys, and a Lookup of the first key)
Counts / CountsBy / GroupBy / GroupBy+reducer / PositionIndex / AssociationMap / Merge   AssociationQ=True, Keys & Lookup OK
FindGraphIsomorphism: <|1 -> a, 2 -> b, 3 -> c|>  AssociationQ=True Keys={1, 2, 3} Lookup1st=a
MapThread over assocs: <|a -> 11, b -> 22|>   Join: <|1 -> 2, 3 -> 4|>
Union: <|1 -> 5, 3 -> 4|>   Intersection: <|3 -> 4|>   Association[Options[Plot]]: well-formed
LetterCounts / CharacterCounts / WordCounts / Dataset / DateValue[..., Association] / ComponentMeasurements: unevaluated

$ ./mathilda /tmp/assocreport/p9.m
undefined (Names[] empty): {LetterCounts, CharacterCounts, WordCounts, Dataset, ComponentMeasurements, KeyDropFrom,
  KeyComplement, KeyIntersection, JoinAcross, Query, Key, Missing, ToAssociations, TakeDrop, ClassifierMeasurements,
  DateObject, DateValue, GraphAssortativity, Dispatch, SubsetMap}
GroupBy-built 1e5: 1e5 Lookups: 0.017053 s
FGI 2000: 2000 Lookups: 0.000158 s  x20: 0.003024 s      <- linear in probes: index built once, O(1) reads
literal 1e5: 1e5 Lookups: 0.01607 s
```

Every association result produced elsewhere in the system is well-formed and gets O(1) lookups. The surface is narrow, though. Strings, datetime and ML return no associations at all. `WordCounts`, `LetterCounts` and `CharacterCounts` are missing, as are `Dataset`/`Query`, and so are `KeyComplement`, `KeyIntersection` and `KeyDropFrom` within the Association family itself. `src/calculus/integrate_goursat.c:16` still says the WL appendix "relies on Association, which Mathilda lacks". That comment is stale.

### 5.5 Tests and leak check

```
$ build/association_tests  > at.log; echo $?   -> 0     236 "Running test" lines, "All Association tests passed."
$ build/compile_assoc_tests > ct.log; echo $?  -> 0     19 tests, "All Compile[] Association (B1-B5) tests passed!"
$ build/bench_assoc; echo $?                   -> 0
  PASS: all operations scaled linearly (ratio < 3.3)
  PASS: single-key lookup is O(1) (ratio < 1.6)
  PASS: interpreter repeated Lookup is O(1)
  PASS: Do-loop over a loop-invariant association is O(1)
  PASS: no operation exceeded 2.5x its baseline cost    (Counts 0.09x, CountsBy 1.27x, GroupBy 0.83x, Merge 0.59x ...)
$ build/mapthread_tests  -> exit 0      $ build/graph_algos_tests -> exit 0
$ build/packed_list_tests -> exit 134 (abort in test_declines_what_it_cannot_represent: ToNDArray[{1, 2.5}] now gives
  {1.0, 2.5}), unrelated to Association, but the abort means the packed-key Association cases at line 573 never run.
```

Coverage map. The table counts `Head[` occurrences in `tests/test_association.c` and `tests/test_compile_assoc.c`:

| Well covered (≥ 4) | Thin (1–3) | Untested in the association suites |
|---|---|---|
| Lookup 17/41, GroupBy 17, KeyValuePattern 16, Association 9, Counts 7/11, GatherBy 7, KeyDrop 5/28, Merge 5, DeleteMissing 5, Keys 4, Values 4/12, KeyUnion 4, AssociationQ 4 | KeyExistsQ/KeyMemberQ/KeyFreeQ 3, KeyTake 3, KeyValueMap 2, AssociationThread 2, AssociateTo 2, KeySort 2, CountsBy 2, KeySortBy/KeyMap/KeySelect/PositionIndex/AssociationMap/Join 1 | Union/Intersection/Complement on associations, MapThread (covered in `test_mapthread.c`), FindGraphIsomorphism's association result, level-spec forms (`Level`, `Replace[..,{1}]`, `Map[f,a,{1}]`, `Apply` at level 1), associations as rule sets, **any packed or NDArray interplay** (0 `PackedArrayQ`/`NDArray` references in `test_association.c`), and non-canonical `AssociationQ` |

All 25 builtins in `assoc.c` appear in at least one test. The untested forms are exactly where section 5.3 found defects.

**Leaks.** `leak.m` exercises construction (literal, list of rules, `AssociationThread` of 1000 entries), mutation (`a[k]=`, `AssociateTo`, `KeyDrop`, nested `a["c","x"]=`, `a[["c","y"]] += 1`, `MapAt`, a `Do` loop of 50 key assignments), `Merge` (Total and Identity), `GroupBy` (with a reducer, and over records), `Counts`, `CountsBy`, `PositionIndex`, `KeySort`, `KeyMap`, `KeySelect`, `KeyUnion`, `KeyValueMap`, `Select`, `Map`, `Total`, bulk `Lookup`, `KeyValuePattern` in `Cases`, `ReplaceAll`, `Normal`, `Append`, and a compiled `KeyDrop`/`Values` call:

```
$ ./mathilda /tmp/assocreport/leak.m
4 <|a -> 112, b -> 20, c -> <|x -> {9, 9, 9}, y -> 4|>, new -> 5, d -> 4|> 7 51 499225.0
$ leaks --atExit -- ./mathilda /tmp/assocreport/leak.m
Process 16130: 33494 nodes malloced for 2459 KB
Process 16130: 0 leaks for 0 total leaked bytes.
$ leaks --atExit -- ./mathilda /tmp/assocreport/p7.m   -> 0 leaks   (malformed-association paths)
$ leaks --atExit -- ./mathilda /tmp/assocreport/p4.m   -> 0 leaks   (Compile/CompileDiagnostics paths)
```

### 5.6 Documentation and refpages

```
$ python3 - <<'EOF'
import sys, re; sys.path.insert(0,'site'); import generate as g
f=g.discover_builtins(); c,s=g.parse_spec_files()
names=sorted(n for n,v in f.items() if v["module"].startswith("src/assoc"))+["KeyValuePattern","DeleteMissing","Key","Missing"]
for n in names: b=s.get(n,{}).get("body",""); print(n, len(re.findall(r"^\s*In\[\d*\]",b,re.M)), len(re.findall(r"^\s*Out\[\d*\]",b,re.M)), len(f.get(n,{}).get("doc","")))
EOF
assoc.c builtins: 25     -> every one has an H2 in docs/spec/builtins/data-structures.md with >= 1 In/Out pair
                            (Association 5, GroupBy 5, KeyDrop 3, Lookup 3, most others 1-2)
KeyValuePattern 5/5, DeleteMissing 1/1   (src/patterns.c)
Key 0/0 doc=0      Missing 0/0 doc=0      <- no H2, no docstring
```

- **Generated pages:** all 27 have one, at `site/docs/documentation/data-structures/<Name>.md` (for example `.../data-structures/Association.md` and `.../KeyValuePattern.md`). `Key` and `Missing` have none.
- **Docstrings:** all 27 are present. Lengths run from 57 characters (KeyValueMap) to 264 (KeyValuePattern). `Information["Key"]` and `Information["Missing"]` both print `No information available for symbol ...`, and `Attributes[Key]` and `Attributes[Missing]` are `{}`.
- **Example accuracy:** I extracted all 49 distinct `In[]/Out[]` pairs from the association H2 sections and ran them through `./mathilda`. 48 matched character-for-character after ignoring whitespace and quotes. The one mismatch (KeyUnion) came from my extractor reading only the first line of a two-line `Out[]`; the full output agrees. The docs are faithful to the binary.

### Key findings

- [strength] **The core suites are green and leak-free.** `association_tests` passed 236 tests and `compile_assoc_tests` passed 19. All five `bench_assoc` scaling and O(1) gates pass. `leaks --atExit` reports 0 leaks over construction, mutation, nested updates, Merge, GroupBy and compiled calls.
- [strength] **Pattern-matcher integration is solid.** `_Association`, `KeyValuePattern` (condition, nested, DownValue binding), `Condition`, `Cases`, `Select`, `Position` and `ReplaceAll` into keys and values all behave. Associations produced elsewhere (FindGraphIsomorphism, MapThread, Union, GroupBy) are well-formed and get O(1) lookups through the lazy index.
- [strength] **The reference docs are complete and correct for the 27 association builtins.** Each has an H2 with In/Out examples, a generated refpage and a docstring, and all 49 doc examples reproduce.
- [risk] **Level-spec operations traverse the internal `Rule` entries.** `Level[a,{1}]` returns rules, `Replace[a, r, {1}]` is a no-op, and `Map[f, a, {1}]`, `Apply[g, a, {1}]` and `a /. (k -> v) -> x` build malformed `Association[...]` nodes. `AssociationQ`, being a bare head test, still reports `True` for them, and none of these forms is tested.
- [risk] **Associations and packed arrays don't mix.** Association is deliberately off `AWARE`, so `Values` and `Keys` never come back packed and a packed value is boxed at `Rule`. `Total[Values[a]]` runs about 170× slower than on the packed list. A visible `NDArray` passed to `AssociationThread`, `AssociationMap` or `GroupBy` is left unevaluated, which CLAUDE.md classes as a wrong answer. There are no association probes in the ND audit tools, and 0 packed-array tests in `test_association.c`.
- [gap] **The Compile subset stops at a cliff.** `Map`, `Select` and `Append` lower only at the top level (they are not producers), and `Keys`, `p[k]` and `Part[p,k]` never lower. A plain list passed as `_Association` is accepted silently. `tools/compile_coverage.py`'s `EXEMPT` entry for Counts and the matching text in `COMPILE_MISSING.md` are stale now that `Counts[v]` compiles.
- [gap] **Missing surface.** There are no `KeyComplement`, `KeyIntersection`, `KeyDropFrom`, `Dataset`/`Query`, `WordCounts`, `LetterCounts`, `CharacterCounts` or `DateValue[..., Association]`, and an association cannot be used as a rule set (`{"a"} /. <|"a" -> 1|>` stays `{"a"}`). `Key` and `Missing` have no docstring, spec section, refpage or attributes, and the `integrate_goursat.c` comment still says Mathilda lacks Association.
- [risk] **One related test binary aborts.** `build/packed_list_tests` exits with status 134 on a `ToNDArray[{1, 2.5}]` regression, which is not Association-specific. The abort stops the run before its packed-key Association cases (`test_packed_list.c:573-574`) execute. They pass when run by hand.
