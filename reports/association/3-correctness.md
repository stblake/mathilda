## 3. Correctness (differential vs Mathematica 15)

**Method.** A battery of 271 deterministic expressions (`/tmp/assocreport/diff/battery.m`) was run unchanged through `./mathilda battery.m` and `wolframscript -file battery.m` (Mathematica 15). Each case prints `id | ToString[expr, InputForm]` through a `HoldRest` helper, so both systems serialise results the same way. Stateful cases run inside `Module` so that no state leaks from one case to the next. Lines were joined by id and compared as exact strings (`mathilda.raw`, `wolfram.raw`, `mismatches.tsv`, `classification.tsv`). Each mismatch was then classified by hand. The WRONG-VALUE cases were reduced to minimal reproducers (`repro.m`, `repro.mathilda`, `repro.wolfram`), and a direct `Print[InputForm[...]]` / `Print[FullForm[...]]` check is in `print.m`.

Coverage, by case-id prefix: construction `con` (24), key types `key` (18), lookup `lk` (34), mutation and copy semantics `mu` (24), functional operations `fn` (99), structure and patterns `st` (60), printing `pr` (12).

### Overall numbers

| | Cases |
|---|---|
| Total cases | **271** |
| Exact match | **211 (77.9%)** |
| Mismatch | **60** |
| WRONG-VALUE | **32** |
| UNEVALUATED | **16** |
| EXTRA (Mathilda evaluates, Mathematica does not) | **5** |
| FORMAT-ONLY | **7** |

Match rate by area: construction 21/24, key types 16/18, lookup 27/34, mutation 17/24, functional 77/99, structure and patterns 45/60, printing 8/12. Differences in message text are not counted. For example, Mathilda prints no `Part::partw` or `AssociationThread::idim` messages, and it prints `AddTo::rvalue` where Mathematica succeeds.

### WRONG-VALUE and UNEVALUATED cases

| ID | Input | Mathilda | Mathematica 15 | Severity |
|---|---|---|---|---|
| con15 | `<\|"a" -> 1, "a" :> 2\|>` | `<\|"a" -> 2\|>` | `<\|"a" :> 2\|>` | WRONG-VALUE |
| con21 | `<\|x -> 1, y -> 2\|> /. x -> z` | `<\|z -> 1, y -> 2\|>` | `<\|x -> 1, y -> 2\|>` | WRONG-VALUE |
| lk07 | `<\|"a" -> 1, "b" -> 2, "c" -> 3\|>[[2 ;; 3]]` | `Missing["KeyAbsent", Span[2, 3]]` | `<\|"b" -> 2, "c" -> 3\|>` | WRONG-VALUE |
| lk08 | `<\|"a" -> 1, "b" -> 2, "c" -> 3\|>[[{1, 3}]]` | `{1, 3}` | `<\|"a" -> 1, "c" -> 3\|>` | WRONG-VALUE |
| lk09 | `<\|"a" -> 1, "b" -> 2, "c" -> 3\|>[[{"a", "c"}]]` | `{1, 3}` | `<\|"a" -> 1, "c" -> 3\|>` | WRONG-VALUE |
| lk23 | `<\|"a" -> 1, "b" -> 2\|>[[All]]` | `Missing["KeyAbsent", All]` | `<\|"a" -> 1, "b" -> 2\|>` | WRONG-VALUE |
| lk30 | `<\|"a" -> 1\|>[Key["a"]]` | `1` | `Missing["KeyAbsent", Key["a"]]` | WRONG-VALUE |
| mu08 | `Module[{a = <\|"a" -> 1, "b" -> 2\|>}, KeyDropFrom[a, "a"]; a]` | `<\|"a" -> 1, "b" -> 2\|>` | `<\|"b" -> 2\|>` | WRONG-VALUE |
| mu09 | `Module[{a = <\|"a" -> 1, "b" -> 2\|>}, KeyDropFrom[a, {"a", "b"}]; a]` | `<\|"a" -> 1, "b" -> 2\|>` | `<\|\|>` | WRONG-VALUE |
| mu11 | `Module[{a = <\|"a" -> 1, "b" -> 2\|>}, a["a"] =.; a]` | `<\|"a" -> 1, "b" -> 2\|>` | `<\|"b" -> 2\|>` | WRONG-VALUE |
| mu14 | `Module[{a = <\|"x" -> <\|"y" -> 1\|>\|>}, a["x"]["z"] = 3; a]` | `<\|"x" -> <\|"y" -> 1\|>\|>` | `<\|"x" -> <\|"y" -> 1, "z" -> 3\|>\|>` | WRONG-VALUE |
| mu17 | `Module[{a = <\|"a" -> 1\|>}, a["a"] += 5; a]` | `<\|"a" -> 1\|>` | `<\|"a" -> 6\|>` | WRONG-VALUE |
| mu18 | `Module[{a = <\|"a" -> 1\|>}, a["a"]++; a]` | `<\|"a" -> 1\|>` | `<\|"a" -> 2\|>` | WRONG-VALUE |
| fn12 | `SortBy[<\|"a" -> 3, "b" -> 1, "c" -> 2\|>, Minus]` | `<\|"b" -> 1, "c" -> 2, "a" -> 3\|>` | `<\|"a" -> 3, "c" -> 2, "b" -> 1\|>` | WRONG-VALUE |
| fn14 | `KeySortBy[<\|"c" -> 1, "a" -> 2, "b" -> 3\|>, Minus @* ToCharacterCode]` | `<\|"a" -> 2, "b" -> 3, "c" -> 1\|>` | `<\|"c" -> 1, "b" -> 3, "a" -> 2\|>` | WRONG-VALUE |
| fn53 | `KeyTake[<\|"a" -> 1, "b" -> 2, "c" -> 3\|>, {"c", "a", "z"}]` | `<\|"a" -> 1, "c" -> 3\|>` | `<\|"c" -> 3, "a" -> 1\|>` | WRONG-VALUE |
| fn59 | `Prepend[<\|"a" -> 1, "b" -> 2\|>, "b" -> 9]` | `<\|"b" -> 2, "a" -> 1\|>` | `<\|"b" -> 9, "a" -> 1\|>` | WRONG-VALUE |
| fn60 | `Delete[<\|"a" -> 1, "b" -> 2\|>, Key["a"]]` | `<\|"a" -> 1, "b" -> 2\|>` | `<\|"b" -> 2\|>` | WRONG-VALUE |
| fn75 | `Catenate[{<\|"a" -> 1\|>, <\|"b" -> 2\|>}]` | `<\|"a" -> 1, "b" -> 2\|>` | `{1, 2}` | WRONG-VALUE |
| fn87 | `Map[f, <\|"a" -> 1\|>, {2}]` | `<\|f["a"] -> f[1]\|>` | `<\|"a" -> 1\|>` | WRONG-VALUE |
| fn90 | `Level[<\|"a" -> 1, "b" -> 2\|>, {1}]` | `{"a" -> 1, "b" -> 2}` | `{1, 2}` | WRONG-VALUE |
| st17 | `<\|"a" -> 1, "b" -> 2\|> /. "a" -> "q"` | `<\|"q" -> 1, "b" -> 2\|>` | `<\|"a" -> 1, "b" -> 2\|>` | WRONG-VALUE |
| st18 | `<\|"a" -> 1, "b" -> 2\|> /. ("a" -> 1) -> ("c" -> 3)` | `<\|"c" -> 3, "b" -> 2\|>` | `<\|"a" -> 1, "b" -> 2\|>` | WRONG-VALUE |
| st30 | `OrderedQ[<\|"a" -> 2, "b" -> 1\|>]` | `True` | `False` | WRONG-VALUE |
| st33 | `Depth[<\|"a" -> 1\|>]` | `3` | `2` | WRONG-VALUE |
| st34 | `Depth[<\|"a" -> <\|"b" -> 1\|>\|>]` | `5` | `3` | WRONG-VALUE |
| st35 | `AtomQ[<\|"a" -> 1\|>]` | `False` | `True` | WRONG-VALUE |
| st43 | `Keys[<\|"a" -> 1, "b" -> 2\|> /. "b" -> "a"]` | `{"a"}` | `{"a", "b"}` | WRONG-VALUE |
| st47 | `Select[{<\|"a" -> 1\|>, <\|"a" -> 5\|>}, #a > 2 &]` | `{}` | `{<\|"a" -> 5\|>}` | WRONG-VALUE |
| st48 | `#a + #b &[<\|"a" -> 1, "b" -> 2\|>]` | `<\|"a" -> 1, "b" -> 2\|> a + <\|"a" -> 1, "b" -> 2\|> b` | `3` | WRONG-VALUE |
| st51 | `#a &[<\|"b" -> 3\|>]` | `<\|"b" -> 3\|> a` | `#a` | WRONG-VALUE |
| st56 | `FreeQ[<\|"a" -> 1\|>, "a"]` | `False` | `True` | WRONG-VALUE |
| con06 | `Association[{{"a" -> 1}, {"b" -> 2}}]` | `Association[{{"a" -> 1}, {"b" -> 2}}]` | `<\|"a" -> 1, "b" -> 2\|>` | UNEVALUATED |
| lk16 | `Lookup["a"][<\|"a" -> 7\|>]` | `Lookup["a"][<\|"a" -> 7\|>]` | `7` | UNEVALUATED |
| mu10 | `Module[{a = <\|"a" -> 1, "b" -> 2\|>}, KeyDropFrom[a, "zz"]]` | `KeyDropFrom[<\|"a" -> 1, "b" -> 2\|>, "zz"]` | `<\|"a" -> 1, "b" -> 2\|>` | UNEVALUATED |
| fn29 | `Merge[{"a" -> 1, "a" -> 5}, Max]` | `Merge[{"a" -> 1, "a" -> 5}, Max]` | `<\|"a" -> 5\|>` | UNEVALUATED |
| fn44 | `AssociationMap[Reverse, <\|"a" -> 1, "b" -> 2\|>]` | `AssociationMap[Reverse, <\|"a" -> 1, "b" -> 2\|>]` | `<\|1 -> "a", 2 -> "b"\|>` | UNEVALUATED |
| fn46 | `PositionIndex[<\|"x" -> 1, "y" -> 2, "z" -> 1\|>]` | `PositionIndex[<\|"x" -> 1, "y" -> 2, "z" -> 1\|>]` | `<\|1 -> {"x", "z"}, 2 -> {"y"}\|>` | UNEVALUATED |
| fn50 | `<\|"a" -> 1, "b" -> 2\|> + 1` | `1 + <\|"a" -> 1, "b" -> 2\|>` | `<\|"a" -> 2, "b" -> 3\|>` | UNEVALUATED |
| fn51 | `<\|"a" -> 1, "b" -> 2\|> * <\|"a" -> 3, "b" -> 4\|>` | `<\|"a" -> 1, "b" -> 2\|> <\|"a" -> 3, "b" -> 4\|>` | `<\|"a" -> 3, "b" -> 8\|>` | UNEVALUATED |
| fn52 | `Sqrt[<\|"a" -> 4\|>]` | `Sqrt[<\|"a" -> 4\|>]` | `<\|"a" -> 2\|>` | UNEVALUATED |
| fn70 | `KeyIntersection[{<\|"a" -> 1, "b" -> 2\|>, <\|"b" -> 3\|>}]` | `KeyIntersection[{<\|"a" -> 1, "b" -> 2\|>, <\|"b" -> 3\|>}]` | `{<\|"b" -> 2\|>, <\|"b" -> 3\|>}` | UNEVALUATED |
| fn71 | `KeyComplement[{<\|"a" -> 1, "b" -> 2\|>, <\|"b" -> 3\|>}]` | `KeyComplement[{<\|"a" -> 1, "b" -> 2\|>, <\|"b" -> 3\|>}]` | `<\|"a" -> 1\|>` | UNEVALUATED |
| fn72 | `Query["a"][<\|"a" -> 1\|>]` | `Query["a"][<\|"a" -> 1\|>]` | `1` | UNEVALUATED |
| fn73 | `Dataset`Dummy; Transpose[{<\|"a" -> 1, "b" -> 2\|>, <\|"a" -> 3, "b" -> 4\|>}]` | `Transpose[{<\|"a" -> 1, "b" -> 2\|>, <\|"a" -> 3, "b" -> 4\|>}]` | `{<\|"a" -> 1, "b" -> 2\|>, <\|"a" -> 3, "b" -> 4\|>}` | UNEVALUATED |
| st21 | `<\|"a" -> 1, "b" -> 2\|> == <\|"b" -> 2, "a" -> 1\|>` | `<\|"a" -> 1, "b" -> 2\|> == <\|"b" -> 2, "a" -> 1\|>` | `False` | UNEVALUATED |
| st24 | `<\|"a" -> 1\|> == <\|"a" -> 1.\|>` | `<\|"a" -> 1\|> == <\|"a" -> 1.0\|>` | `True` | UNEVALUATED |
| st50 | `Slot["a"] &[<\|"a" -> 3\|>]` | `#"a"` | `3` | UNEVALUATED |

### Minimal reproducers for WRONG-VALUE, grouped by root cause

Mathematica 15's answer is on the right of each arrow.

1. **ReplaceAll rewrites keys.** Mathematica treats an Association as atomic, so `/.` reaches only the values.
   - `<|1 -> 2|> /. 1 -> 0` gives `<|0 -> 2|>`; Mathematica gives `<|1 -> 2|>`.
   - `Keys[<|a -> 1, b -> 2|> /. b -> a]` gives `{a}`; Mathematica gives `{a, b}`. Silent data loss: two keys collapse into one.
   - The same pattern drives con21, st17, st18 and st43. A related case is `FreeQ[<|a -> 1|>, a]`, which gives `False`; Mathematica gives `True` (st56).
2. **Part with Span, a list of indices, or All does not return a sub-Association.**
   - `<|a -> 1, b -> 2|>[[1 ;; 1]]` gives `Missing["KeyAbsent", Span[1, 1]]`; Mathematica gives `<|a -> 1|>`.
   - `<|a -> 1, b -> 2|>[[{1}]]` gives `{1}`; Mathematica gives `<|a -> 1|>`.
   - `<|a -> 1|>[[All]]` gives `Missing["KeyAbsent", All]`.
3. **In-place mutation other than `a[k] = v` is a silent no-op.**
   - `s = <|a -> 1|>; KeyDropFrom[s, a]; s` gives `<|a -> 1|>`.
   - `s[a] =.` has no effect.
   - `s[a] += 1` fails with `AddTo::rvalue`.
   - `s[a][b] = 1` has no effect.
   - `a[k] = v`, `a[[k]] = v`, `a[k1, k2] = v`, `AssociateTo` and `AppendTo[a[k], ...]` all work.
4. **Minus does not evaluate on numbers.** `Minus[3]` stays `Minus[3]`; Mathematica gives `-3`. This is a core arithmetic bug, not an Association bug. `SortBy[{1, 2}, Minus]` returns `{1, 2}`, and the same failure causes the Association SortBy and KeySortBy mismatches (fn12, fn14). `SortBy[..., -# &]` is correct.
5. **Named slots are misparsed.** `#a` parses as `Times[Slot[1], a]` (`FullForm[Hold[#a]]` shows this).
   - `#a &[<|"a" -> 1|>]` gives `<|"a" -> 1|> a`; Mathematica gives `1`.
   - `Slot["a"]` is not applied to associations either (st50).
   - These cases break the common `Select[list, #key > x &]` idiom.
6. **Associations are not atomic in structural functions.**
   - `AtomQ[<|a -> 1|>]` gives `False`; Mathematica gives `True`.
   - `Depth[<|a -> 1|>]` gives `3`; Mathematica gives `2`.
   - `Level[<|a -> 1|>, {1}]` gives `{a -> 1}`; Mathematica gives `{1}`.
   - `Map[f, <|a -> 1|>, {2}]` gives `<|f[a] -> f[1]|>`; Mathematica leaves it unchanged.
   - `OrderedQ[<|a -> 2, b -> 1|>]` gives `True`; Mathematica gives `False` because it compares values. Mathilda appears to compare the rules, and so the keys.
7. **Key-order and duplicate-key handling.**
   - `KeyTake[<|a -> 1, b -> 2|>, {b, a}]` keeps the source order; Mathematica uses the requested order.
   - `Prepend[<|a -> 1|>, a -> 2]` keeps the old value; Mathematica gives `<|a -> 2|>`.
   - `<|a -> 1, a :> 2|>` drops the RuleDelayed; Mathematica gives `<|a :> 2|>`.
8. **Other operations.**
   - `Delete[<|a -> 1|>, Key[a]]` is a no-op.
   - `Catenate[{<|a -> 1|>}]` returns an Association; Mathematica returns the list of values `{1}`.
   - `<|a -> 1|>[Key[a]]` gives `1`; Mathematica gives `Missing["KeyAbsent", Key[a]]`.

### EXTRA (Mathilda evaluates where Mathematica does not)

| ID | Input | Mathilda | Mathematica 15 | Severity |
|---|---|---|---|---|
| fn81 | `Fold[Plus, 0, <\|"a" -> 1, "b" -> 2\|>]` | `3` | `Fold[Plus, 0, <\|"a" -> 1, "b" -> 2\|>]` | EXTRA |
| fn82 | `Accumulate[<\|"a" -> 1, "b" -> 2, "c" -> 3\|>]` | `<\|"a" -> 1, "b" -> 3, "c" -> 6\|>` | `Accumulate[<\|"a" -> 1, "b" -> 2, "c" -> 3\|>]` | EXTRA |
| fn83 | `Differences[<\|"a" -> 1, "b" -> 4\|>]` | `<\|"b" -> 3\|>` | `Differences[<\|"a" -> 1, "b" -> 4\|>]` | EXTRA |
| fn92 | `Partition[<\|"a" -> 1, "b" -> 2, "c" -> 3, "d" -> 4\|>, 2]` | `<\|"a" -> 1, "b" -> 2, "c" -> 3, "d" -> 4\|>` | `Partition[<\|"a" -> 1, "b" -> 2, "c" -> 3, "d" -> 4\|>, 2]` | EXTRA |
| st19 | `<\|"a" -> x, "b" -> x^2\|> /. x -> 3` | `<\|"a" -> 3, "b" -> 9\|>` | `<\|"a" -> 3, "b" -> 3^2\|>` | EXTRA |

`Fold`, `Accumulate` and `Differences` on an Association are reasonable extensions. `Partition` silently returning its input is not; it should stay unevaluated. For st19, Mathematica does not re-evaluate values after `/.` and gives `3^2`. Mathilda gives `9`, which is arguably more useful but diverges from Mathematica.

### Format-only differences

- **Real formatting.** Machine reals print as `1.0` where Mathematica prints `1.` (key01, key05, and the `<|"a" -> 1.0|>` in st24).
- **Part.** An unevaluated Part prints as `Part[a, 3]` rather than `a[[3]]` (lk33). Mathilda also emits no `Part::partw` message.
- **FullForm.** `ToString[FullForm[x], InputForm]` strips the `FullForm` wrapper in Mathilda; Mathematica keeps it (pr02, pr03). `FullForm` output itself is otherwise identical: `Association[Rule["a", 1], RuleDelayed["b", 2]]`.
- **Plain ToString.** `ToString[<|"a" -> 1|>]` (OutputForm) keeps the string quotes in Mathilda; Mathematica drops them (pr05).
- **Nested strings.** In `ToString[ToString[..., InputForm], InputForm]`, the inner quotes are not escaped as `\"` (pr06). The same quoting bug applies to any string, not only associations.
- **Direct Print.** `Print[InputForm[<|"a" -> 1|>]]` in Mathilda prints `<|a -> 1|>` without quotes. The output is therefore not valid InputForm (`print.m`).

### Key findings

- [strength] **Core construction and lookup are sound.** Construction follows Mathematica's rules: last key wins, `<||>`, rules mixed with lists of rules, RuleDelayed values held and evaluated on lookup, and non-rule arguments left unevaluated. Key identity is correct (`1`, `1.`, `"1"` and `1/2` are distinct keys), as are `Missing["KeyAbsent", k]`, `Lookup` with defaults and key lists, nested `a[["x", "y"]]` and `a["x", "y"]`, and positional `[[1]]`. Mathilda matches on 211 of 271 cases (77.9%).
- [strength] **Copy semantics and basic mutation match Mathematica.** After `b = a`, mutating `a` leaves `b` unchanged. `a[k] = v`, `a[[k]] = v`, nested `a[k1, k2] = v`, `AssociateTo` and `AppendTo[a[k], ...]` all behave identically. So do `Merge`, `GroupBy` (all forms, including the 3-argument form), `Counts`, `CountsBy`, `AssociationThread`, `KeyValuePattern`, `SameQ` and `Hash`.
- [risk] **ReplaceAll rewrites keys and can silently merge them.** `<|a -> 1, b -> 2|> /. b -> a` returns `<|a -> 2|>`. Mathilda does not treat Associations as atomic: AtomQ, Depth, Level, Map at level 2, FreeQ and OrderedQ all look inside key/value rules. This one root cause produces 11 of the 32 WRONG-VALUE cases.
- [risk] **Some mutations are silent no-ops.** KeyDropFrom, `a[k] =.`, `Delete[a, Key[k]]`, `a[k1][k2] = v` and compound assignment (`+=`, `++`) leave the association unchanged, mostly without any message, so programs keep running on stale state.
- [risk] **Named slots are misparsed.** `#a` parses as `Slot[1]*a`, so `#key &` idioms such as `Select[rows, #a > 2 &]` return wrong answers rather than failing. Separately, `Minus[3]` does not evaluate, which breaks `SortBy`/`KeySortBy` with `Minus` on lists and associations alike.
- [gap] **Part does not return sub-Associations.** Span, lists of indices or keys, and `All` all fail: `a[[2 ;; 3]]` yields `Missing[...]` and `a[[{1, 3}]]` yields a bare list of values.
- [gap] **Listable arithmetic and several heads do not thread over Associations.** `<|...|> + 1`, `Sqrt[<|...|>]` and `==` between associations stay unevaluated. So do `KeyIntersection`, `KeyComplement`, `Query`, `Lookup[k]` (operator form), `Merge` on a list of rules, `AssociationMap` and `PositionIndex` on an Association, and `Association` of a list of lists of rules. That is 16 UNEVALUATED cases in all.
- [gap] **Printing is mostly right, with string-quoting gaps.** InputForm and FullForm of associations match Mathematica. The differences are the `1.0` vs `1.` real format, unescaped quotes in nested InputForm strings, and `Print[InputForm[...]]` dropping string quotes.
