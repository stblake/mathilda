# Data Structures

Associations (`<| key -> value, ... |>`) are Mathilda's key–value data structure,
modelled on the Wolfram Language. Keys are unique and insertion-ordered (the
first occurrence fixes a key's position, the last occurrence its value), and
bulk operations are backed by a hash index for amortised `O(n)` construction,
grouping and lookup.

**Single-key access is amortised `O(1)`.** A persistent key→position hash index
is cached on the association (see `src/assoc_index.h`) and built on first
single-key read, so `Lookup[a, k]`, `KeyExistsQ`/`KeyMemberQ`/`KeyFreeQ`,
`a[[key]]`, `a[key]` and `Key[k][a]` do not scan. The index reuses spare space in
the expression node, so it costs no extra memory per non-association value, and it
is transparent — `SameQ`, `Union` and `Sort` see the association exactly as
before. Repeated access to a bound association stays `O(1)` per probe — including
inside `Do`/`Table`/`Fold`: a loop-invariant association (or any value built from
literals under `List`/`Association`/`Rule`/`Complex`/`Rational`) is a GROUND fixed
point that survives the iterator rebinding's eval-clock churn, so it is re-checked
in `O(1)` per iteration rather than re-canonicalised in `O(n)` (see the
evaluator's rule-epoch / GROUND-flag mechanism in `src/eval.c`).

### The association toolchain at a glance

| Category | Operations |
|----------|------------|
| **Construct** | `Association` / `<\|…\|>`, `AssociationThread`, `AssociationMap`, `Association[{rules}]` (splice) |
| **Access** | `assoc[key]`, `assoc[k1, k2]` (nested), `assoc[[key]]` / `[[Key[k]]]` / `[[i]]`, `Lookup`, `Normal` |
| **Write** | `assoc[[key]] = v`, `AssociateTo`, `Append`, `Prepend`, `Join` (merge) |
| **Keys & values** | `Keys`, `Values`, `KeyMap`, `KeyValueMap`, `KeyUnion`, `KeyIntersection`, `KeyComplement` |
| **Presence** | `KeyExistsQ`, `KeyMemberQ`, `KeyFreeQ`, `MemberQ` |
| **Transform** | `Map`, `Select`, `Discard`, `KeySelect`, `KeyTake`, `KeyDrop`, `DeleteMissing`, `Transpose`, arithmetic (Listable threading) |
| **Combine** | `JoinAcross`, `Merge`, `Join`, `SubsetQ`, `CountDistinct`, `CountDistinctBy` |
| **Aggregate** | `Total`, `Min`, `Max`, `Mean`, `Counts`, `CountsBy`, `GroupBy` (+reducer, `key->val`), `Gather`, `GatherBy`, `Merge`, `PositionIndex` |
| **Order / rank** | `Sort`, `SortBy` (multi-key), `ReverseSort`, `ReverseSortBy`, `KeySort`, `KeySortBy`, `MaximalBy`, `MinimalBy`, `TakeLargest`(`By`), `TakeSmallest`(`By`), `Reverse` |
| **Iterate / reduce** | `Table`/`Do`/`Sum`/`Product` (`{v, assoc}`), `Fold`, `FoldList`, `Scan`, `Cases`, `Count`, `DeleteCases`, `Position`, `FirstPosition`, `SelectFirst`, `FirstCase`, `AllTrue`, `AnyTrue`, `NoneTrue` |
| **Structural** | `First`, `Last`, `Rest`, `Most`, `Take`, `Drop`, `Length` |
| **Patterns** | `KeyValuePattern` (destructure/match, incl. in function definitions) |

Functions that consume or reduce a collection operate on an association's
**values**, keeping keys aligned (see *Design notes: value threading* below).

## Association
Represents a mapping from keys to values, written `<|k1 -> v1, k2 -> v2, ...|>`
or `Association[k1 -> v1, ...]`.
- Arguments may be rules, lists of rules, or other associations (which are spliced).
  Lists may nest to any depth (`Association[{{a -> 1}, {b -> 2}}]`), and
  `Splice[{rules}]` splices into a literal (see [Splice](#splice)).
- Duplicate keys collapse with last-value-wins, preserving first-occurrence order.
- `assoc[[key]]`, `assoc[[Key[key]]]` and `assoc[[i]]` extract values (a missing
  key gives `Missing["KeyAbsent", key]`); `assoc[[{k1, k2, ...}]]` gives the list
  of those values.
- **`assoc[key]`** — an association applied as a function looks the key up (the
  idiomatic accessor), giving the value or `Missing["KeyAbsent", key]`; `assoc[Key[k]]`
  is the explicit form.
- **`Key[k][assoc]`** — the curried complement: `Key[k]` applied to an association
  extracts the value at `k`. This is the record-field extractor for pipelines —
  `GroupBy[records, Key["field"]]`, `SortBy[records, Key["field"]]`,
  `Map[Key["field"], records]`.

```mathematica
In[1]:= <|"a" -> 1, "b" -> 2|>
Out[1]= <|"a" -> 1, "b" -> 2|>

In[2]:= <|"a" -> 1, "b" -> 2, "a" -> 99|>
Out[2]= <|"a" -> 99, "b" -> 2|>

In[3]:= <|"a" -> 10, "b" -> 20|>[["b"]]
Out[3]= 20

In[4]:= <|"a" -> 10, "b" -> 20|>["a"]
Out[4]= 10

In[5]:= <|"a" -> <|"b" -> 5|>|>["a", "b"]
Out[5]= 5

In[6]:= Association[{{"a" -> 1}, {"b" -> 2, {"c" -> 3}}}]
Out[6]= <|"a" -> 1, "b" -> 2, "c" -> 3|>
```

## AssociationQ
Tests whether an expression is an association.

```mathematica
In[1]:= AssociationQ[<|"a" -> 1|>]
Out[1]= True

In[2]:= AssociationQ[{1, 2, 3}]
Out[2]= False
```

## Keys
Gives the list of keys of an association (or a list of rules).
- `Keys[assoc]` gives `{k1, k2, ...}`; also accepts a rule or a list of rules.
- `Keys[{assoc1, assoc2, ...}]` threads over a list (to any depth) of
  associations and lists of rules.
- `Keys[assoc, f]` wraps each key: `{f[k1], f[k2], ...}`.

**Features**: an invalid list element gives `Keys::invrl` and leaves the call
unevaluated. Mathematica has no level argument (`Keys[a, 2]` is the `f` form
with `f = 2`), and neither does Mathilda.

```mathematica
In[1]:= Keys[<|"a" -> 1, "b" -> 2|>]
Out[1]= {"a", "b"}

In[2]:= Keys[<|"a" -> 1, "b" -> 2|>, f]
Out[2]= {f["a"], f["b"]}

In[3]:= Keys[{<|"a" -> 1|>, <|"b" -> 2, "c" -> 3|>}]
Out[3]= {{"a"}, {"b", "c"}}
```

## Values
Gives the list of values of an association (or a list of rules).
- `Values[assoc]` gives `{v1, v2, ...}`; also accepts a rule or a list of rules.
- `Values[{assoc1, assoc2, ...}]` threads like `Keys`.
- `Values[assoc, f]` wraps each value: `{f[v1], f[v2], ...}`.

**Features**: same threading and error rules as `Keys`.

```mathematica
In[1]:= Values[<|"a" -> 1, "b" -> 2|>]
Out[1]= {1, 2}

In[2]:= Values[<|"a" -> 1, "b" -> 2|>, f]
Out[2]= {f[1], f[2]}

In[3]:= Values[{<|"a" -> 1|>, {"b" -> 2}}]
Out[3]= {{1}, {2}}
```

## Lookup
Looks up the value stored under a key.
- `Lookup[assoc, key]` gives the value, or `Missing["KeyAbsent", key]`.
- `Lookup[assoc, key, default]` uses `default` when the key is absent.
- `Lookup[assoc, {k1, k2, ...}]` looks up several keys with a single hash-index
  build (`O(n + m)`).
- `Lookup[{a1, a2, ...}, key]` threads over a list of associations, extracting
  the key from each (a key/default thread through) — handy for pulling one field
  out of a column of records. The list may mix associations and lists of rules:
  `Lookup[{<|"a" -> 1|>, {"a" -> 3}, {}}, "a", 0]` gives `{1, 3, 0}`.
- Also accepts a bare list of rules (like `Keys`/`Values`).

```mathematica
In[1]:= Lookup[<|"a" -> 1, "b" -> 2|>, "b"]
Out[1]= 2

In[2]:= Lookup[<|"a" -> 1|>, "z", 0]
Out[2]= 0

In[3]:= Lookup[{<|"a" -> 1, "b" -> 2|>, <|"a" -> 3|>}, "a", 0]
Out[3]= {1, 3}
```

## KeyExistsQ, KeyMemberQ, KeyFreeQ
`KeyExistsQ`/`KeyMemberQ` test whether a key is present; `KeyFreeQ` is the
complement (True when the key is absent). All three accept an association or a
bare list of rules.

```mathematica
In[1]:= KeyExistsQ[<|"a" -> 1|>, "a"]
Out[1]= True

In[2]:= KeyFreeQ[<|"a" -> 1|>, "b"]
Out[2]= True
```

## KeyDrop
Gives an association with the specified key or keys removed (order preserved).

```mathematica
In[1]:= KeyDrop[<|"a" -> 1, "b" -> 2, "c" -> 3|>, "b"]
Out[1]= <|"a" -> 1, "c" -> 3|>

In[2]:= KeyDrop[<|"a" -> 1, "b" -> 2, "c" -> 3|>, {"a", "c"}]
Out[2]= <|"b" -> 2|>
```

`KeyDrop` (like `KeyTake`) also threads over a list of associations — dropping
the keys from each record:

```mathematica
In[3]:= KeyDrop[{<|"a" -> 1, "b" -> 2|>, <|"a" -> 3, "b" -> 4|>}, "b"]
Out[3]= {<|"a" -> 1|>, <|"a" -> 3|>}
```

## KeyTake
Gives the association of only the specified keys (association order preserved).
Over a list of associations it threads, keeping the keys in each record.

```mathematica
In[1]:= KeyTake[<|"a" -> 1, "b" -> 2, "c" -> 3|>, {"c", "a"}]
Out[1]= <|"a" -> 1, "c" -> 3|>

In[2]:= KeyTake[{<|"a" -> 1, "b" -> 2|>, <|"a" -> 3, "b" -> 4|>}, {"a"}]
Out[2]= {<|"a" -> 1|>, <|"a" -> 3|>}
```

## KeyUnion
`KeyUnion[{assoc1, assoc2, ...}]` pads every association to the union of all
their keys (in first-appearance order), filling a key absent from an
association with `Missing["KeyAbsent", key]`. The result is the list of
equalised associations — a common preparation step before row-wise or tabular
processing. Hash-indexed: the union and every rebuild are `O(n)`.

```mathematica
In[1]:= KeyUnion[{<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "c" -> 4|>}]
Out[1]= {<|"a" -> 1, "b" -> 2, "c" -> Missing["KeyAbsent", "c"]|>,
         <|"a" -> Missing["KeyAbsent", "a"], "b" -> 3, "c" -> 4|>}

In[2]:= KeyUnion[{<|"a" -> 1|>, <|"b" -> 2|>}, 0 &]
Out[2]= {<|"a" -> 1, "b" -> 0|>, <|"a" -> 0, "b" -> 2|>}
```

- `KeyUnion[{...}, f]` fills an absent key `k` with `f[k]` instead of `Missing`.
- The list may also hold rules and lists of rules.

## KeyValueMap
Applies `f` to each key–value pair, giving `{f[k1, v1], f[k2, v2], ...}`.

```mathematica
In[1]:= KeyValueMap[f, <|a -> 1, b -> 2|>]
Out[1]= {f[a, 1], f[b, 2]}

In[2]:= KeyValueMap[Plus, <|1 -> 10, 2 -> 20|>]
Out[2]= {11, 22}
```

## AssociationThread
Builds an association from parallel key and value lists.

```mathematica
In[1]:= AssociationThread[{"a", "b"}, {1, 2}]
Out[1]= <|"a" -> 1, "b" -> 2|>

In[2]:= AssociationThread[{"a", "b"} -> {1, 2}]
Out[2]= <|"a" -> 1, "b" -> 2|>
```

## Counts
Tallies each distinct element of a list, giving `<|element -> count, ...|>`.
Hash-indexed: a single `O(n)` pass. Over an association it tallies the values.

```mathematica
In[1]:= Counts[{1, 2, 2, 3, 3, 3}]
Out[1]= <|1 -> 1, 2 -> 2, 3 -> 3|>

In[2]:= Counts[<|"a" -> 1, "b" -> 1, "c" -> 2|>]
Out[2]= <|1 -> 2, 2 -> 1|>
```

> **Packed arrays.** Over a packed list or an `NDArray`, `Counts` keys on the
> raw machine words — the same count `Tally` does — and relabels the result as
> an `Association`, at one pass over the *distinct* values.

## GroupBy
Groups the elements of a list by the value of `f` applied to each element.
Hash-indexed grouping in `O(n)` plus the cost of `f`.
- `GroupBy[list, f, g]` applies the reducer `g` to each group (split-apply-combine).
- `GroupBy[list, keyfn -> valfn]` groups by `keyfn[x]` but collects `valfn[x]` in
  each group; `GroupBy[list, keyfn -> valfn, g]` then reduces — the full
  group / extract / reduce pipeline in one call.
- `GroupBy[assoc, f]` groups an association's **entries** by `f[value]` into
  sub-associations (keys preserved); `GroupBy[assoc, f, g]` reduces each
  sub-association (so `GroupBy[assoc, f, Total]` composes with value-threading
  `Total`). Distinct from `GroupBy[records, Key["field"]]`, which groups a
  *list* of record-associations by a field.

```mathematica
In[1]:= GroupBy[{1, 2, 3, 4, 5, 6}, EvenQ]
Out[1]= <|False -> {1, 3, 5}, True -> {2, 4, 6}|>

In[2]:= GroupBy[Range[10], EvenQ, Total]
Out[2]= <|False -> 25, True -> 30|>

In[3]:= GroupBy[{{"x", 1}, {"y", 2}, {"x", 3}}, First -> Last, Total]
Out[3]= <|"x" -> 4, "y" -> 2|>

In[4]:= GroupBy[<|"a" -> 1, "b" -> 2, "c" -> 3, "d" -> 4|>, EvenQ]
Out[4]= <|False -> <|"a" -> 1, "c" -> 3|>, True -> <|"b" -> 2, "d" -> 4|>|>

In[5]:= GroupBy[<|"a" -> 1, "b" -> 2, "c" -> 3, "d" -> 4|>, EvenQ, Total]
Out[5]= <|False -> 4, True -> 6|>
```

## Gather
Gathers identical elements into sublists, in order of each element's first
occurrence; within a sublist elements keep their input order. Equal elements are
collected from anywhere in the list, not only from adjacent runs — this is the
difference from [`Split`](structural-manipulation.md#split). Equivalent to
`GatherBy[list, Identity]`, and implemented on the same grouping engine
(hash-indexed, O(n)) with the identity key applied directly rather than as `n`
`Identity[x]` applications.

```mathematica
In[1]:= Gather[{1, 7, 3, 7, 2, 3, 9}]
Out[1]= {{1}, {7, 7}, {3, 3}, {2}, {9}}

In[2]:= Gather[{a, b, a}]
Out[2]= {{a, a}, {b}}

In[3]:= Gather[{}]
Out[3]= {}
```

Notes:
- Group order is first-occurrence order, not sorted order.
- Elements are grouped by structural identity (`expr_eq`), so `1` and `1.0` are
  distinct: `Gather[{1, 1.0}]` gives `{{1}, {1.0}}`.
- Over an association, the entries are gathered by value into sub-associations
  (keys preserved), matching `GatherBy[assoc, Identity]`.
- Anything other than a list or association, and any arity other than 1, is
  returned unevaluated.

## GatherBy
Gathers the elements with equal `f[element]` into sublists, in first-appearance
order (like `GroupBy` but returning the groups as a plain list of lists). Over an
association it gathers the **entries** by `f[value]` into sub-associations (keys
preserved), returned as a list — `GroupBy[assoc, f]` without the outer keys.

```mathematica
In[1]:= GatherBy[{1, 2, 3, 4, 5, 6}, EvenQ]
Out[1]= {{1, 3, 5}, {2, 4, 6}}

In[2]:= GatherBy[<|"a" -> 1, "b" -> 2, "c" -> 3, "d" -> 4|>, EvenQ]
Out[2]= {<|"a" -> 1, "c" -> 3|>, <|"b" -> 2, "d" -> 4|>}
```

## Merge
Combines several associations, applying `f` to the list of values collected for
each key (in first-seen key order). The first argument may be a list *or* an
association of associations. The list may also hold rules and lists of rules;
each rule contributes one value (`Merge[{a -> 1, a -> 2}, Total]` is `<|a -> 3|>`).

```mathematica
In[1]:= Merge[{<|"a" -> 1|>, <|"a" -> 2, "b" -> 3|>}, Total]
Out[1]= <|"a" -> 3, "b" -> 3|>

In[2]:= Merge[<|"g1" -> <|"a" -> 1|>, "g2" -> <|"a" -> 2, "b" -> 3|>|>, Total]
Out[2]= <|"a" -> 3, "b" -> 3|>

In[3]:= Merge[{"a" -> 1, "b" -> 2, "a" -> 3}, Total]
Out[3]= <|"a" -> 4, "b" -> 2|>
```

## AssociateTo
Adds or updates key–value pairs in the association held by a symbol, modifying
the symbol in place (like `AppendTo`). Has attribute `HoldFirst`.

```mathematica
In[1]:= asc = <|"a" -> 1|>; AssociateTo[asc, "b" -> 2]; asc
Out[1]= <|"a" -> 1, "b" -> 2|>
```

## Part assignment
Associations support in-place update through `Part`: `a[[key]] = val` updates an
existing key or adds a new one, `a[[Key[k]]] = val` targets a key explicitly,
and `a[[i]] = val` updates the i-th value positionally. `a[[All]] = val`,
`a[[i ;; j]] = val` and `a[[{k1, ...}]] = val` assign into every / spanned /
listed entry's value. Multi-index assignment descends into nested associations
and lists (`a[[k1, k2]] = val`).

```mathematica
In[1]:= a = <|"x" -> 1, "y" -> 2|>; a[["x"]] = 99; a
Out[1]= <|"x" -> 99, "y" -> 2|>

In[2]:= a[["z"]] = 5; a
Out[2]= <|"x" -> 99, "y" -> 2, "z" -> 5|>

In[3]:= n = <|"p" -> <|"q" -> 1|>|>; n[["p", "q"]] = 42; n
Out[3]= <|"p" -> <|"q" -> 42|>|>
```

## Min, Max, MinMax
`Min[assoc]` / `Max[assoc]` give the extreme values; `MinMax[assoc]` gives
`{min, max}` in one shot (the value-range, e.g. for plot bounds).

```mathematica
In[1]:= MinMax[<|"a" -> 3, "b" -> 1, "c" -> 9|>]
Out[1]= {1, 9}
```

## Apply (@@)
`f @@ assoc` uses the association's **values** as `f`'s arguments —
`f @@ <|k1 -> v1, ...|>` is `f[v1, ...]` (so `Total[assoc]` is `Plus @@ assoc`,
and `Apply[List, assoc]` is `Values[assoc]`).

```mathematica
In[1]:= Plus @@ <|"a" -> 1, "b" -> 2, "c" -> 3|>
Out[1]= 6

In[2]:= Apply[List, <|"a" -> 1, "b" -> 2|>]
Out[2]= {1, 2}
```

## Mean
`Mean[assoc]` averages the association's values (like `Total`/`Min`/`Max`).

```mathematica
In[1]:= Mean[<|"a" -> 2, "b" -> 4, "c" -> 6|>]
Out[1]= 4
```

## Median, Variance, StandardDeviation
The rest of the descriptive-statistics family reduces over an association's
values too, matching `Mean`.

```mathematica
In[1]:= Median[<|"a" -> 1, "b" -> 3, "c" -> 5|>]
Out[1]= 3

In[2]:= Variance[<|"a" -> 2, "b" -> 4, "c" -> 6|>]
Out[2]= 4

In[3]:= StandardDeviation[<|"a" -> 2, "b" -> 4, "c" -> 6|>]
Out[3]= 2
```

## Accumulate, Differences, Ratios, FoldList
These windowed transforms act on the values but — unlike `Total` or `Fold`,
which collapse — keep the association shape. `Accumulate[assoc]` gives the
running totals with every key retained; `Differences[assoc]` and `Ratios[assoc]`
give the successive differences / ratios, keyed by the trailing key of each pair
(so the leading key drops, `n` entries → `n - 1`); `FoldList[f, assoc]` pairs
each key with the running result of `f` (`n → n`), the general form of
`Accumulate`. All share one `assoc_rekey_over_values` mechanism.

```mathematica
In[1]:= Accumulate[<|"a" -> 1, "b" -> 2, "c" -> 3|>]
Out[1]= <|"a" -> 1, "b" -> 3, "c" -> 6|>

In[2]:= Differences[<|"a" -> 1, "b" -> 4, "c" -> 9|>]
Out[2]= <|"b" -> 3, "c" -> 5|>

In[3]:= Ratios[<|"a" -> 1, "b" -> 2, "c" -> 6|>]
Out[3]= <|"b" -> 2, "c" -> 3|>

In[4]:= FoldList[Times, <|"a" -> 2, "b" -> 3, "c" -> 4|>]
Out[4]= <|"a" -> 2, "b" -> 6, "c" -> 24|>
```

## Tally, Commonest
`Tally[assoc]` tallies the association's values as `{value, count}` pairs, and
`Commonest[assoc]` returns the most frequent value(s) — the value-oriented
siblings of `Counts` / `CountsBy`.

```mathematica
In[1]:= Tally[<|"a" -> 1, "b" -> 1, "c" -> 2|>]
Out[1]= {{1, 2}, {2, 1}}

In[2]:= Commonest[<|"a" -> 1, "b" -> 1, "c" -> 2|>]
Out[2]= {1}
```

## TakeWhile, LengthWhile
Act on the leading run of **values**. `TakeWhile[assoc, crit]` keeps the leading
entries whose value satisfies `crit` (keys preserved); `LengthWhile[assoc, crit]`
counts them. (See functional-programming for the list forms.)

```mathematica
In[1]:= TakeWhile[<|"a" -> 1, "b" -> 2, "c" -> 5, "d" -> 1|>, # < 3 &]
Out[1]= <|"a" -> 1, "b" -> 2|>

In[2]:= LengthWhile[<|"a" -> 1, "b" -> 2, "c" -> 5|>, # < 3 &]
Out[2]= 2
```

## Cases, Count, DeleteCases
Pattern operations act on the **values** of an association. `Cases` returns the
list of matching values, `Count` their number, and `DeleteCases` returns the
association with the matching-value entries removed.

```mathematica
In[1]:= Cases[<|"a" -> 1, "b" -> 2, "c" -> 3|>, x_ /; x > 1]
Out[1]= {2, 3}

In[2]:= Count[<|"a" -> 1, "b" -> 2, "c" -> 3|>, x_ /; x > 1]
Out[2]= 2

In[3]:= DeleteCases[<|"a" -> 1, "b" -> 2, "c" -> 3|>, x_ /; x > 1]
Out[3]= <|"a" -> 1|>
```

## Top-N by value (TakeLargest, TakeSmallest)
`TakeLargest[assoc, n]` / `TakeSmallest[assoc, n]` return the `n` entries with
the largest / smallest values (as an association, ranked). `TakeLargestBy` /
`TakeSmallestBy` rank by `f` of each value.

```mathematica
In[1]:= TakeLargest[<|"a" -> 3, "b" -> 9, "c" -> 1, "d" -> 6|>, 2]
Out[1]= <|"b" -> 9, "d" -> 6|>
```

## Extremes by value (MaximalBy, MinimalBy)
`MaximalBy[assoc, f]` / `MinimalBy[assoc, f]` return the entries whose value
maximises / minimises `f` (all ties), as an association.

```mathematica
In[1]:= MaximalBy[<|"a" -> 1, "b" -> 3, "c" -> 3|>, Identity]
Out[1]= <|"b" -> 3, "c" -> 3|>
```

## DeleteMissing
Removes `Missing[...]` elements — the natural cleanup after a multi-key
`Lookup`. Over an association it drops entries whose value is `Missing[...]`.

- `DeleteMissing[expr, n]` removes `Missing[...]` at levels 1 through `n`
  (`n` a positive integer or `Infinity`); association values count as one level
  below the association.
- `DeleteMissing[expr, n, d]` removes the elements at levels 1..`n` that contain
  a `Missing[...]` at depth `d` or less (`d = 0`: the element itself). Levels are
  processed deepest first, so a sublist emptied at a deeper level is kept.
- A first argument that is not a list or an association gives
  `DeleteMissing::invrp`.

```mathematica
In[1]:= DeleteMissing[Lookup[<|"a" -> 1, "b" -> 2|>, {"a", "z", "b"}]]
Out[1]= {1, 2}

In[2]:= DeleteMissing[{1, {Missing[], 2}}, 2]
Out[2]= {1, {2}}

In[3]:= DeleteMissing[{{1, Missing[]}, {2}}, 1, 1]
Out[3]= {{2}}
```

## DeleteDuplicates
Over an association, keeps the first entry for each distinct value and returns
an association (values are compared, keys along for the ride). Hash-indexed: a
single `O(n)` pass.

```mathematica
In[1]:= DeleteDuplicates[<|"a" -> 1, "b" -> 1, "c" -> 2, "d" -> 2, "e" -> 3|>]
Out[1]= <|"a" -> 1, "c" -> 2, "e" -> 3|>
```

## DeleteDuplicatesBy
`DeleteDuplicatesBy[expr, f]` keeps the first element for each distinct
`f[element]`, preserving order. Over an association, `f` is applied to the
values and the surviving entries keep their keys.

```mathematica
In[1]:= DeleteDuplicatesBy[{1, 12, 3, 14, 5}, EvenQ]
Out[1]= {1, 12}

In[2]:= DeleteDuplicatesBy[<|"a" -> 1, "b" -> 12, "c" -> 3|>, EvenQ]
Out[2]= <|"a" -> 1, "b" -> 12|>
```

## ReplacePart
`ReplacePart[assoc, {Key[k]} -> v]` (or bare `Key[k]`, or a positional `i`)
replaces the value at that key; several positions may be given at once. It is
replace-only — an absent key is left unchanged (unlike Part assignment, which
would add it).

```mathematica
In[1]:= ReplacePart[<|"a" -> 1, "b" -> 2, "c" -> 3|>, {{Key["a"]} -> 10, {Key["c"]} -> 30}]
Out[1]= <|"a" -> 10, "b" -> 2, "c" -> 30|>
```

## Delete
`Delete[assoc, {Key[k]}]` removes an entry by key position; a nested position
(`{Key[k1], Key[k2]}`) descends into inner associations, and a list of positions
(`{{Key[k1]}, {Key[k2]}}`) removes several. (`KeyDrop` is the key-oriented
equivalent.)

```mathematica
In[1]:= Delete[<|"a" -> 1, "b" -> 2, "c" -> 3|>, {Key["b"]}]
Out[1]= <|"a" -> 1, "c" -> 3|>

In[2]:= Delete[<|"a" -> <|"x" -> 5, "y" -> 6|>|>, {Key["a"], Key["x"]}]
Out[2]= <|"a" -> <|"y" -> 6|>|>
```

## MapAt
`MapAt[f, assoc, key]` (or `Key[k]`, a positional index, or a nested position
`{Key[k1], Key[k2]}`) applies `f` to the value at that position — composing with
the `{Key[k]}` positions `Position` returns. `All`, `Span` and the head index `0`
address the entry list, as they do for a `List`. An absent key or an
out-of-range index leaves `MapAt` unevaluated. `ReplaceAt` accepts exactly the
same association positions.

```mathematica
In[1]:= MapAt[#^2 &, <|"a" -> 3, "b" -> 4|>, "b"]
Out[1]= <|"a" -> 3, "b" -> 16|>

In[2]:= p = <|"a" -> 1, "b" -> 9|>; MapAt[-# &, p, First[Position[p, 9]]]
Out[2]= <|"a" -> 1, "b" -> -9|>

In[3]:= MapAt[f, <|"a" -> 1, "b" -> 2|>, All]
Out[3]= <|"a" -> f[1], "b" -> f[2]|>

In[4]:= MapAt[f, <|"a" -> 1, "b" -> 2, "c" -> 3|>, 1 ;; 2]
Out[4]= <|"a" -> f[1], "b" -> f[2], "c" -> 3|>

In[5]:= ReplaceAt[<|"a" -> 1, "b" -> 2|>, 1 -> 9, Key["a"]]
Out[5]= <|"a" -> 9, "b" -> 2|>
```

## Position
`Position[assoc, patt]` gives the positions of matches inside the **values** as
`{Key[k], subpos...}` (Wolfram semantics), descending into nested values; the
positions round-trip through `Part`/`Extract`.

```mathematica
In[1]:= Position[<|"a" -> 1, "b" -> 2, "c" -> 1|>, 1]
Out[1]= {{Key["a"]}, {Key["c"]}}

In[2]:= Position[<|"a" -> {1, 2}, "b" -> 3, "c" -> 1|>, 1]
Out[2]= {{Key["a"], 1}, {Key["c"]}}

In[3]:= Extract[<|"a" -> {10, 20}|>, {Key["a"], 2}]
Out[3]= 20
```

`FirstPosition[assoc, patt]` returns the first such `{Key[k], subpos...}` (or
`Missing["NotFound"]` / a supplied default), sharing this Key-remapping with
`Position`.

## Predicate tests on values (AllTrue, AnyTrue, NoneTrue, MemberQ)
Predicate tests apply to the **values** of an association (and to the elements
of a list). `AllTrue`/`AnyTrue`/`NoneTrue` are left unevaluated if a test
result is neither `True` nor `False`.

```mathematica
In[1]:= AllTrue[<|"a" -> 2, "b" -> 4|>, EvenQ]
Out[1]= True

In[2]:= AnyTrue[{1, 3, 4}, EvenQ]
Out[2]= True

In[3]:= NoneTrue[{1, 3, 5}, EvenQ]
Out[3]= True

In[4]:= MemberQ[<|"a" -> 1, "b" -> 2|>, 2]
Out[4]= True
```

## Map, Select
`Map` and `Select` thread over the **values** of an association, preserving keys
(matching Wolfram semantics) — `Map[f, <|k -> v|>]` gives `<|k -> f[v]|>`, and
`Select` keeps the entries whose value satisfies the predicate.

```mathematica
In[1]:= Map[#^2 &, <|"x" -> 3, "y" -> 4|>]
Out[1]= <|"x" -> 9, "y" -> 16|>

In[2]:= Select[<|"a" -> 1, "b" -> 2, "c" -> 3|>, # > 1 &]
Out[2]= <|"b" -> 2, "c" -> 3|>
```

## Sort, SortBy, Total, Min, Max, Join
Ordering and aggregation act on the **values** of an association: `Sort` orders
the entries by value and `SortBy[assoc, f]` by `f` of each value (keys follow),
while `Total`/`Min`/`Max` reduce over the values. `Join` merges associations
(later values win).

```mathematica
In[1]:= Sort[<|"a" -> 3, "b" -> 1, "c" -> 2|>]
Out[1]= <|"b" -> 1, "c" -> 2, "a" -> 3|>

In[2]:= SortBy[<|"a" -> {9}, "b" -> {1}|>, First]
Out[2]= <|"b" -> {1}, "a" -> {9}|>

In[3]:= Total[<|"a" -> 3, "b" -> 1, "c" -> 2|>]
Out[3]= 6

In[4]:= Join[<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "c" -> 4|>]
Out[4]= <|"a" -> 1, "b" -> 3, "c" -> 4|>
```

## KeySort
Sorts an association into canonical key order.

```mathematica
In[1]:= KeySort[<|"c" -> 3, "a" -> 1, "b" -> 2|>]
Out[1]= <|"a" -> 1, "b" -> 2, "c" -> 3|>
```

## KeySortBy
Sorts an association by `f` applied to each key (stable on ties).

```mathematica
In[1]:= KeySortBy[<|"bbb" -> 1, "a" -> 2, "cc" -> 3|>, StringLength]
Out[1]= <|"a" -> 2, "cc" -> 3, "bbb" -> 1|>
```

## KeyMap
Applies `f` to every key, keeping values. Keys that collide collapse with
last-value-wins.

```mathematica
In[1]:= KeyMap[f, <|1 -> 10, 2 -> 20|>]
Out[1]= <|f[1] -> 10, f[2] -> 20|>
```

## KeySelect
Keeps the entries whose **key** satisfies the predicate.

```mathematica
In[1]:= KeySelect[<|1 -> 10, 2 -> 20, 3 -> 30|>, EvenQ]
Out[1]= <|2 -> 20|>
```

## CountsBy
Tallies elements of a list by `f[element]`, giving `<|f[x] -> count, ...|>`.
Hash-indexed, `O(n)`.

```mathematica
In[1]:= CountsBy[Range[10], EvenQ]
Out[1]= <|False -> 5, True -> 5|>
```

## PositionIndex
Maps each distinct element of a list to the list of 1-based positions where it
occurs — `<|value -> {positions}|>`. Hash-indexed, `O(n)`.

`PositionIndex[assoc]` maps each distinct value to the list of keys at which
it occurs.

```mathematica
In[1]:= PositionIndex[{a, b, a, c, a, b}]
Out[1]= <|a -> {1, 3, 5}, b -> {2, 6}, c -> {4}|>

In[2]:= PositionIndex[<|"a" -> x, "b" -> y, "c" -> x|>]
Out[2]= <|x -> {"a", "c"}, y -> {"b"}|>
```

## AssociationMap
Builds `<|k1 -> f[k1], k2 -> f[k2], ...|>` from a list of keys.
- `AssociationMap[f, assoc]` applies `f` to each rule `k -> v`; each result may
  be a rule, a list of rules, an association, or `{}`/`Nothing` (dropped), and
  together they form the new association (last value wins on a repeated key).
  Any other result gives `AssociationMap::invrlf`.

```mathematica
In[1]:= AssociationMap[#^2 &, {1, 2, 3, 4}]
Out[1]= <|1 -> 1, 2 -> 4, 3 -> 9, 4 -> 16|>

In[2]:= AssociationMap[Reverse, <|"a" -> 1, "b" -> 2|>]
Out[2]= <|1 -> "a", 2 -> "b"|>
```

## Iterating an association
Iterator-driven builtins (`Table`, `Do`, `Sum`, `Product`) walk an association's
**values** when given `{var, assoc}`.

```mathematica
In[1]:= Table[v^2, {v, <|"a" -> 2, "b" -> 3|>}]
Out[1]= {4, 9}

In[2]:= Sum[v, {v, <|"a" -> 10, "b" -> 20|>}]
Out[2]= 30
```

## Design notes: value threading

A unifying principle runs through the association operations above: **functions
that consume or reduce a collection operate on an association's values, keeping
keys aligned**, matching the Wolfram Language.

- Element-wise (`Map`, `Select`) return an association: `Map[f, <|k -> v|>]` is `<|k -> f[v]|>`.
- Reductions (`Total`, `Min`, `Max`, `Mean`) return a scalar over the values — and are
  defined to be exactly the list reduction over `Values[assoc]`, so empty-collection and
  edge behaviour is identical to the list case (`Total[<||>]` behaves as `Total[{}]`).
- Ordering (`Sort`, `SortBy`, `ReverseSort`, `MaximalBy`, `TakeLargest`) reorders/selects
  entries by value, keeping the entries intact.
- Pattern/predicate ops (`Cases`, `Count`, `DeleteCases`, `MemberQ`, `AllTrue`, `SelectFirst`,
  `FirstCase`) test the values; `Cases`/`SelectFirst`/`FirstCase` return the matching values,
  while `DeleteCases`/`Select` return the surviving entries as an association.

Because these all thread through a single shared helper (`assoc_apply_over_values`)
or the same key-preserving rebuild, they compose predictably into pipelines —
e.g. `ReverseSort[GroupBy[txns, First, Total[#[[All, 2]]] &]]` groups, reduces
and ranks in one expression (see [`../../../examples/association-showcase.md`](../../../examples/association-showcase.md)).

## Append, Prepend
`Append`/`Prepend` extend an association with new entries (non-mutating siblings
of `AssociateTo`); an existing key is updated in place, preserving order.

```mathematica
In[1]:= Append[<|"a" -> 1, "b" -> 2|>, "c" -> 3]
Out[1]= <|"a" -> 1, "b" -> 2, "c" -> 3|>

In[2]:= Append[<|"a" -> 1|>, "a" -> 99]
Out[2]= <|"a" -> 99|>
```

## KeyValuePattern
`KeyValuePattern[{k1 -> p1, ...}]` is a **pattern** that matches an association
(or a list of rules) containing keys matching `k1, ...` with values matching
`p1, ...`. Value patterns may bind, so associations can be destructured in rules
and used to filter records.

```mathematica
In[1]:= MatchQ[<|"a" -> 1, "b" -> 2|>, KeyValuePattern[{"a" -> _}]]
Out[1]= True

In[2]:= Replace[<|"a" -> 5, "b" -> 2|>, KeyValuePattern[{"a" -> v_}] :> v]
Out[2]= 5

In[3]:= Cases[{<|"t" -> 1|>, <|"t" -> 2|>, <|"x" -> 3|>}, KeyValuePattern[{"t" -> _}]]
Out[3]= {<|"t" -> 1|>, <|"t" -> 2|>}
```

Requirements are matched with backtracking (so shared bound variables resolve),
and the pattern composes with a `/;` condition over its bindings:

```mathematica
In[4]:= Cases[{<|"p" -> 3|>, <|"p" -> 9|>}, KeyValuePattern[{"p" -> v_}] /; v > 5 :> v]
Out[4]= {9}
```

`KeyValuePattern` also works in function definitions, so associations can be
destructured directly in a rule's left-hand side:

```mathematica
In[4]:= area[KeyValuePattern[{"w" -> w_, "h" -> h_}]] := w h; area[<|"w" -> 3, "h" -> 4|>]
Out[4]= 12
```

## First, Last, Rest, Most, Take, Drop
Structural extractors follow Wolfram semantics on associations: `First`/`Last`
give the first/last **value**, while `Rest`/`Most`/`Take`/`Drop` slice **entries**
and return an association (order preserved). `First[expr, default]` and
`Last[expr, default]` return `default` when `expr` has no elements (an empty
association, empty list, or an atom) — this works for lists too.

```mathematica
In[1]:= First[<|"a" -> 10, "b" -> 20|>]
Out[1]= 10

In[1b]:= First[<||>, 0]
Out[1b]= 0

In[2]:= Rest[<|"a" -> 10, "b" -> 20, "c" -> 30|>]
Out[2]= <|"b" -> 20, "c" -> 30|>

In[3]:= Take[<|"a" -> 1, "b" -> 2, "c" -> 3|>, 2]
Out[3]= <|"a" -> 1, "b" -> 2|>
```

## Arithmetic threading over associations
Listable functions (`Plus`, `Times`, `Power`, `Sqrt`, `Sin`, `Abs`, `N`,
`StringLength`, any user symbol with the `Listable` attribute, ...) map over
the values of association arguments, keeping the keys.
- Non-association arguments are repeated for every key: `<|a -> 1|> + 1` is
  `<|a -> 2|>`.
- Several association arguments must have the same keys **in the same order**;
  otherwise `Association::incmp` is issued and the call stays unevaluated (as in
  Mathematica 15, where `<|a -> 1, b -> 2|> + <|b -> 1, a -> 2|>` is also
  rejected).
- When a `List` argument is present too, List threading runs first, so the list
  is the outer level: `<|a -> 1, b -> 2|> + {10, 20}` is a list of two
  associations.
- `Total`, `Mean` and the other reductions over a *list* of associations
  follow from this (they reduce with `Plus`), while `Total[assoc]` and
  `Mean[assoc]` still reduce the values of a single association.
- A `RuleDelayed` entry of the first association stays delayed.

**Features**: the rule lives in the evaluator's Listable step
(`assoc_thread_listable`, `src/assoc_ops.c`). One pass classifies the
arguments, so calls with a `List` argument, or with neither kind, take exactly
the path they took before.

```mathematica
In[1]:= <|"a" -> 1, "b" -> 2|> + 10
Out[1]= <|"a" -> 11, "b" -> 12|>

In[2]:= 2 <|"a" -> 1, "b" -> 2|> + <|"a" -> 100, "b" -> 200|>
Out[2]= <|"a" -> 102, "b" -> 204|>

In[3]:= Sqrt[<|"a" -> 4, "b" -> 9|>]
Out[3]= <|"a" -> 2, "b" -> 3|>

In[4]:= <|"a" -> 1, "b" -> 2|> + {10, 20}
Out[4]= {<|"a" -> 11, "b" -> 12|>, <|"a" -> 21, "b" -> 22|>}

In[5]:= Total[{<|"a" -> 1, "b" -> 2|>, <|"a" -> 3, "b" -> 4|>}]
Out[5]= <|"a" -> 4, "b" -> 6|>

In[6]:= Mean[{<|"a" -> 1, "b" -> 2|>, <|"a" -> 3, "b" -> 4|>}]
Out[6]= <|"a" -> 2, "b" -> 3|>

In[7]:= <|"a" -> 1, "b" -> 2|> + <|"b" -> 1, "a" -> 2|>
Association::incmp: The arguments <|"a" -> 1, "b" -> 2|> and <|"b" -> 1, "a" -> 2|> in <|"a" -> 1, "b" -> 2|> + <|"b" -> 1, "a" -> 2|> are incompatible.
Out[7]= <|"a" -> 1, "b" -> 2|> + <|"b" -> 1, "a" -> 2|>
```

## Key
`Key[k]` names the key `k` of an association.
- `Key[k][assoc]` gives the value at `k`, or `Missing["KeyAbsent", k]`.
- `assoc[[Key[k]]]`, `assoc[Key[k]]` and `Lookup[assoc, Key[k]]` read the key
  literally, which matters when `k` is an integer (not a position) or a list
  (not a list of keys).
- Key specs in `GroupBy`, `SortBy`, `JoinAcross`, `MapAt`,
  `ReplacePart` and `Extract` accept `Key[k]`.

**Features**: attributes `{Protected}`.

```mathematica
In[1]:= Key["b"][<|"a" -> 1, "b" -> 2|>]
Out[1]= 2

In[2]:= Lookup[<|1 -> "one", {1, 2} -> "pair"|>, Key[{1, 2}]]
Out[2]= "pair"

In[3]:= Key["z"][<|"a" -> 1|>]
Out[3]= Missing["KeyAbsent", "z"]
```

## Missing
`Missing[]`, `Missing["reason"]` and `Missing["reason", data]` represent missing
data. Key access on an absent key gives `Missing["KeyAbsent", k]`; `KeyUnion`
fills gaps with `Missing["KeyAbsent", k]` and `JoinAcross` with
`Missing["Unmatched"]` / `Missing["NotAvailable"]`.
- Test with [`MissingQ`](#missingq); remove with [`DeleteMissing`](#deletemissing).

**Features**: an inert head; attributes `{Protected}` (Mathematica also sets
`ReadProtected`, which Mathilda does not implement).

```mathematica
In[1]:= <|"a" -> 1|>["q"]
Out[1]= Missing["KeyAbsent", "q"]

In[2]:= Missing["NotAvailable"]
Out[2]= Missing["NotAvailable"]

In[3]:= DeleteMissing[{1, Missing["NotAvailable"], 3}]
Out[3]= {1, 3}
```

## MissingQ
`MissingQ[expr]` gives `True` if `expr` has head `Missing`, and `False`
otherwise.

**Features**: attributes `{Protected}`; any other argument count gives
`MissingQ::argx`.

```mathematica
In[1]:= MissingQ[Missing["KeyAbsent", "z"]]
Out[1]= True

In[2]:= MissingQ[<|"a" -> 1|>["z"]]
Out[2]= True

In[3]:= Select[{1, Missing[], 3}, Not @* MissingQ]
Out[3]= {1, 3}
```

## KeyIntersection
`KeyIntersection[{assoc1, assoc2, ...}]` restricts every association to the
keys common to all of them, in the key order of `assoc1`.
- Elements may also be rules or lists of rules.
- `KeyIntersection[{}]` is `{}`; anything else that is not a list of
  associations or rules gives `KeyIntersection::invar`.

**Features**: attributes `{Protected}`; `O(total size)` via the persistent key
index.

```mathematica
In[1]:= KeyIntersection[{<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "c" -> 4|>}]
Out[1]= {<|"b" -> 2|>, <|"b" -> 3|>}

In[2]:= KeyIntersection[{<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "a" -> 4|>}]
Out[2]= {<|"a" -> 1, "b" -> 2|>, <|"a" -> 4, "b" -> 3|>}
```

## KeyComplement
`KeyComplement[{assoc1, assoc2, ...}]` gives the entries of `assoc1` whose keys
occur in none of the other associations.
- Elements may also be rules or lists of rules; an empty list gives
  `KeyComplement::empt`.

**Features**: attributes `{Protected}`; `O(total size)`.

```mathematica
In[1]:= KeyComplement[{<|"a" -> 1, "b" -> 2, "c" -> 3|>, <|"b" -> 0|>, <|"c" -> 0|>}]
Out[1]= <|"a" -> 1|>
```

## Discard
The complement of `Select`.
- `Discard[expr, crit]` drops the elements `e` for which `crit[e]` is `True`.
- `Discard[expr, crit, n]` drops only the first `n` of them (`n` a non-negative
  integer or `Infinity`).
- On an association `crit` tests the values and the surviving entries keep
  their keys; any other non-atomic expression keeps its head.

**Features**: attributes `{Protected}`. An atomic first argument gives
`Discard::normal`, a bad `n` gives `Discard::innf`. The operator form
`Discard[crit]` is provided by the generic curried-form mechanism.

```mathematica
In[1]:= Discard[{1, 2, 3, 4, 5, 6}, EvenQ]
Out[1]= {1, 3, 5}

In[2]:= Discard[<|"a" -> 1, "b" -> 2, "c" -> 4|>, EvenQ]
Out[2]= <|"a" -> 1|>

In[3]:= Discard[{1, 2, 3, 4, 5, 6}, EvenQ, 2]
Out[3]= {1, 3, 5, 6}
```

## CountDistinct
`CountDistinct[expr]` gives the number of distinct elements of `expr` (of its
values, for an association).

**Features**: attributes `{Protected}`; one hash pass, `O(n)`; elements are
distinct under `SameQ`, so `1` and `1.` count separately. An atomic argument
gives `CountDistinct::normal`.

```mathematica
In[1]:= CountDistinct[{1, 2, 1, 3, 2}]
Out[1]= 3

In[2]:= CountDistinct[<|"a" -> 1, "b" -> 1, "c" -> 2|>]
Out[2]= 2
```

## CountDistinctBy
`CountDistinctBy[expr, f]` gives the number of distinct values of `f[e]` over the
elements `e` of `expr` (the values, for an association).

**Features**: attributes `{Protected}`; `f` is applied once per element.

```mathematica
In[1]:= CountDistinctBy[{1, 2, 3, 4, 5}, EvenQ]
Out[1]= 2

In[2]:= CountDistinctBy[{"apple", "avocado", "banana"}, StringTake[#, 1] &]
Out[2]= 2
```

## SubsetQ
`SubsetQ[a, b]` gives `True` if every element of `b` occurs in `a`.
- Multiplicity is ignored: `SubsetQ[{1, 1, 2}, {1, 1, 1}]` is `True`.
- Lists and associations (compared by value) may be mixed; any other two
  expressions must share their head (`SubsetQ[f[1, 2], f[1]]`), otherwise
  `SubsetQ::heads`.

**Features**: attributes `{Protected}`; `O(|a| + |b|)` with a hash set.

```mathematica
In[1]:= SubsetQ[{1, 2, 3}, {3, 1}]
Out[1]= True

In[2]:= SubsetQ[{1, 2}, {1, 4}]
Out[2]= False

In[3]:= SubsetQ[<|"a" -> 1, "b" -> 2|>, {2}]
Out[3]= True
```

## Splice
`Splice[{e1, e2, ...}]` is replaced by the sequence `e1, e2, ...` when it appears
as an argument of a `List` or an `Association`.
- `Splice[list, h]` splices into any head matching the pattern `h`
  (`Splice[list, _]` everywhere). As in Mathematica, an explicit `h` never
  splices into an `Association`.
- Anywhere else it stays inert, so `f[Splice[{1, 2}]]` is unchanged.

**Features**: attributes `{Protected}`. Handled in the evaluator's
Sequence-flattening step, only once a `Splice` argument has been seen, so no
other call pays for it; heads with `SequenceHold`/`HoldAllComplete` are left
alone.

```mathematica
In[1]:= {1, Splice[{2, 3}], 4}
Out[1]= {1, 2, 3, 4}

In[2]:= <|"a" -> 1, Splice[{"b" -> 2, "c" -> 3}]|>
Out[2]= <|"a" -> 1, "b" -> 2, "c" -> 3|>

In[3]:= Table[Splice[{i, -i}], {i, 3}]
Out[3]= {1, -1, 2, -2, 3, -3}

In[4]:= f[1, Splice[{2, 3}]]
Out[4]= f[1, Splice[{2, 3}]]

In[5]:= f[1, Splice[{2, 3}, _]]
Out[5]= f[1, 2, 3]
```

## JoinAcross
Joins two lists of associations on shared key values (a relational join).
- `JoinAcross[{a1, ...}, {b1, ...}, spec]` merges each `ai` with every `bj`
  whose join keys agree (inner join). `spec` is a key `k`, `Key[k]`, `k1 -> k2`
  (keys named differently on the two sides), or a list of these (join on
  several keys).
- `JoinAcross[..., spec, "Inner" | "Left" | "Right" | "Outer"]` keeps unmatched
  left and/or right rows; an unmatched row gets `Missing["Unmatched"]` for the
  other side's keys.
- Option `KeyCollisionFunction -> Left` (default) | `Right` | `f` resolves a
  non-join key present on both sides; `f[k]` gives the pair of keys under which
  the left and right values are kept.

**Features**: attributes `{Protected}`; `Options[JoinAcross]` is
`{KeyCollisionFunction -> Left}`. The right rows are hash-indexed on the
join-key tuple, so a join costs `O(|left| + |right| + |output|)`. Output order
matches Mathematica 15: matched pairs (left-major), then unmatched left rows,
then unmatched right rows; rows whose key sets then differ are padded to the
union with `Missing["NotAvailable"]`.

```mathematica
In[1]:= JoinAcross[{<|"id" -> 1, "name" -> "Ada"|>, <|"id" -> 2, "name" -> "Bob"|>}, {<|"id" -> 1, "dept" -> "R&D"|>, <|"id" -> 3, "dept" -> "Ops"|>}, Key["id"]]
Out[1]= {<|"id" -> 1, "name" -> "Ada", "dept" -> "R&D"|>}

In[2]:= JoinAcross[{<|"id" -> 1, "name" -> "Ada"|>, <|"id" -> 2, "name" -> "Bob"|>}, {<|"id" -> 1, "dept" -> "R&D"|>, <|"id" -> 3, "dept" -> "Ops"|>}, "id", "Outer"]
Out[2]= {<|"id" -> 1, "name" -> "Ada", "dept" -> "R&D"|>, <|"id" -> 2, "name" -> "Bob", "dept" -> Missing["Unmatched"]|>, <|"id" -> 3, "name" -> Missing["Unmatched"], "dept" -> "Ops"|>}

In[3]:= JoinAcross[{<|"k" -> 1, "x" -> 10|>}, {<|"key" -> 1, "y" -> 20|>}, "k" -> "key"]
Out[3]= {<|"k" -> 1, "x" -> 10, "key" -> 1, "y" -> 20|>}

In[4]:= JoinAcross[{<|"a" -> 1, "v" -> "L"|>}, {<|"a" -> 1, "v" -> "R"|>}, "a", KeyCollisionFunction -> Right]
Out[4]= {<|"a" -> 1, "v" -> "R"|>}
```

## ApplyTo
`ApplyTo[x, f]` sets `x` to `f[x]` and returns the new value (`x //= f` in
Mathematica's notation).
- `x` may be a symbol with a value, a part `s[[i]]`, or an association entry
  `s[key]`; the write-back goes through `Set`.
- A target without a value gives `ApplyTo::rvalue` and stays unevaluated.

**Features**: attributes `{HoldFirst, Protected}`.

```mathematica
In[1]:= x = 5; ApplyTo[x, #^2 &]; x
Out[1]= 25

In[2]:= r = <|"n" -> 1|>; ApplyTo[r["n"], # + 10 &]; r
Out[2]= <|"n" -> 11|>
```

## AssociationComap
`AssociationComap[{f1, f2, ...}, x]` gives `<|f1 -> f1[x], f2 -> f2[x], ...|>`.

**Features**: attributes `{Protected}`; a first argument that is not a list
gives `AssociationComap::invl`. The operator form `AssociationComap[{f, ...}]`
is provided by the generic curried-form mechanism.

```mathematica
In[1]:= AssociationComap[{Min, Max, Length}, {3, 1, 2}]
Out[1]= <|Min -> 1, Max -> 3, Length -> 3|>
```

## Transpose (associations)
`Transpose[<|k1 -> <|j1 -> v11, ...|>, ...|>]` swaps the two key levels of an
association of associations.
- The inner associations must share their keys in the same order, otherwise
  `Transpose::nmtx`.
- A list of associations is a vector (associations are atomic at that level),
  which `Transpose` returns unchanged.

```mathematica
In[1]:= Transpose[<|"a" -> <|"x" -> 1, "y" -> 2|>, "b" -> <|"x" -> 3, "y" -> 4|>|>]
Out[1]= <|"x" -> <|"a" -> 1, "b" -> 3|>, "y" -> <|"a" -> 2, "b" -> 4|>|>

In[2]:= Transpose[{<|"x" -> 1|>, <|"x" -> 2|>}]
Out[2]= {<|"x" -> 1|>, <|"x" -> 2|>}
```

## Normal (associations)
`Normal[expr]` converts every association in `expr` to its list of rules, at
any depth, but does not descend into the values of a converted association
(Mathematica's rule). `Normal[expr, Association]` converts only associations.
Held expressions (`Hold`, `HoldComplete`, ...) are left alone.

```mathematica
In[1]:= Normal[{<|"a" -> 1|>, f[<|"b" -> 2|>]}]
Out[1]= {{"a" -> 1}, f[{"b" -> 2}]}

In[2]:= Normal[<|"a" -> <|"b" -> 1|>|>]
Out[2]= {"a" -> <|"b" -> 1|>}
```
