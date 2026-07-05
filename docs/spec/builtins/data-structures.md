# Data Structures

Associations (`<| key -> value, ... |>`) are Mathilda's key–value data structure,
modelled on the Wolfram Language. Keys are unique and insertion-ordered (the
first occurrence fixes a key's position, the last occurrence its value), and
bulk operations are backed by a hash index for amortised `O(n)` construction,
grouping and lookup.

## Association
Represents a mapping from keys to values, written `<|k1 -> v1, k2 -> v2, ...|>`
or `Association[k1 -> v1, ...]`.
- Arguments may be rules, lists of rules, or other associations (which are spliced).
- Duplicate keys collapse with last-value-wins, preserving first-occurrence order.
- `assoc[[key]]`, `assoc[[Key[key]]]` and `assoc[[i]]` extract values (a missing
  key gives `Missing["KeyAbsent", key]`).

```mathematica
In[1]:= <|"a" -> 1, "b" -> 2|>
Out[1]= <|"a" -> 1, "b" -> 2|>

In[2]:= <|"a" -> 1, "b" -> 2, "a" -> 99|>
Out[2]= <|"a" -> 99, "b" -> 2|>

In[3]:= <|"a" -> 10, "b" -> 20|>[["b"]]
Out[3]= 20
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

```mathematica
In[1]:= Keys[<|"a" -> 1, "b" -> 2|>]
Out[1]= {"a", "b"}
```

## Values
Gives the list of values of an association (or a list of rules).

```mathematica
In[1]:= Values[<|"a" -> 1, "b" -> 2|>]
Out[1]= {1, 2}
```

## Lookup
Looks up the value stored under a key.
- `Lookup[assoc, key]` gives the value, or `Missing["KeyAbsent", key]`.
- `Lookup[assoc, key, default]` uses `default` when the key is absent.
- `Lookup[assoc, {k1, k2, ...}]` looks up several keys with a single hash-index
  build (`O(n + m)`).

```mathematica
In[1]:= Lookup[<|"a" -> 1, "b" -> 2|>, "b"]
Out[1]= 2

In[2]:= Lookup[<|"a" -> 1|>, "z", 0]
Out[2]= 0
```

## KeyExistsQ
Tests whether a key is present in an association.

```mathematica
In[1]:= KeyExistsQ[<|"a" -> 1|>, "a"]
Out[1]= True

In[2]:= KeyExistsQ[<|"a" -> 1|>, "b"]
Out[2]= False
```

## KeyDrop
Gives an association with the specified key or keys removed (order preserved).

```mathematica
In[1]:= KeyDrop[<|"a" -> 1, "b" -> 2, "c" -> 3|>, "b"]
Out[1]= <|"a" -> 1, "c" -> 3|>

In[2]:= KeyDrop[<|"a" -> 1, "b" -> 2, "c" -> 3|>, {"a", "c"}]
Out[2]= <|"b" -> 2|>
```

## KeyTake
Gives the association of only the specified keys (association order preserved).

```mathematica
In[1]:= KeyTake[<|"a" -> 1, "b" -> 2, "c" -> 3|>, {"c", "a"}]
Out[1]= <|"a" -> 1, "c" -> 3|>
```

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
Hash-indexed: a single `O(n)` pass.

```mathematica
In[1]:= Counts[{1, 2, 2, 3, 3, 3}]
Out[1]= <|1 -> 1, 2 -> 2, 3 -> 3|>
```

## GroupBy
Groups the elements of a list by the value of `f` applied to each element.
Hash-indexed grouping in `O(n)` plus the cost of `f`.

```mathematica
In[1]:= GroupBy[{1, 2, 3, 4, 5, 6}, EvenQ]
Out[1]= <|False -> {1, 3, 5}, True -> {2, 4, 6}|>
```

## Merge
Combines several associations, applying `f` to the list of values collected for
each key (in first-seen key order).

```mathematica
In[1]:= Merge[{<|"a" -> 1|>, <|"a" -> 2, "b" -> 3|>}, Total]
Out[1]= <|"a" -> 3, "b" -> 3|>
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
and `a[[i]] = val` updates the i-th value positionally. Multi-index assignment
descends into nested associations and lists (`a[[k1, k2]] = val`).

```mathematica
In[1]:= a = <|"x" -> 1, "y" -> 2|>; a[["x"]] = 99; a
Out[1]= <|"x" -> 99, "y" -> 2|>

In[2]:= a[["z"]] = 5; a
Out[2]= <|"x" -> 99, "y" -> 2, "z" -> 5|>

In[3]:= n = <|"p" -> <|"q" -> 1|>|>; n[["p", "q"]] = 42; n
Out[3]= <|"p" -> <|"q" -> 42|>|>
```

## Mean
`Mean[assoc]` averages the association's values (like `Total`/`Min`/`Max`).

```mathematica
In[1]:= Mean[<|"a" -> 2, "b" -> 4, "c" -> 6|>]
Out[1]= 4
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

```mathematica
In[1]:= PositionIndex[{a, b, a, c, a, b}]
Out[1]= <|a -> {1, 3, 5}, b -> {2, 6}, c -> {4}|>
```

## AssociationMap
Builds `<|k1 -> f[k1], k2 -> f[k2], ...|>` from a list of keys.

```mathematica
In[1]:= AssociationMap[#^2 &, {1, 2, 3, 4}]
Out[1]= <|1 -> 1, 2 -> 4, 3 -> 9, 4 -> 16|>
```
