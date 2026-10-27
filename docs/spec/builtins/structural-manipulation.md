# Structural Manipulation

## Part
Extracts subparts of an expression.
- `expr[[i]]` or `Part[expr, i]`: Extracts the $i$-th element.
- `expr[[i, j, ...]]` or `Part[expr, i, j, ...]]`: Extracts nested parts.
- `expr[[{i, j, ...}]]`: Extracts a list of specific parts.
- `expr[[All]]`: Represents all elements at a given level.

**Features**: 
- Supports negative indices to count from the end (`-1` is the last element).
- `expr[[0]]` returns the `Head` of the expression. This is permitted even for atomic expressions.
- Mapping `All` across a dimension allows column extraction from matrices.

```mathematica
In[1]:= {a, b, c, d}[[2]]
Out[1]= b

In[2]:= {a, b, c, d}[[-1]]
Out[2]= d

In[3]:= 123[[0]]
Out[3]= Integer
```

## Extract
Extracts the part of an expression at the position specified by `pos`.
- `Extract[expr, pos]`
- `Extract[expr, {pos1, pos2, ...}]`
- `Extract[expr, pos, h]`
- `Extract[pos]`

**Features**:
- Position specifications have the same form as those returned by `Position`.
- `Extract[expr, {i, j, ...}]` is equivalent to `Part[expr, i, j, ...]`.
- `pos` can be of the more general form `{part1, part2, ...}` where `parti` are `Part` specifications such as an integer `i`, `All` or `Span`.
- You can use `Extract[expr, ..., Hold]` to extract parts without evaluation.
- `Extract[expr, {}]` is `{}` (an empty list of positions), and `Extract[expr, {{}}]`
  is `{expr}` (the one empty path) (v0.586).
- Reads a [packed list](packed-arrays.md) through the buffer, the same gather
  `Part` uses, so pulling three elements out of a 10^6-element array costs three
  reads rather than a materialisation of the whole array. An element from an
  integer buffer comes back as an exact `Integer`. `Extract[expr, {0}]` is head
  extraction and is unaffected.

## Span
- `i;;j`: Represents a span of elements `i` through `j`.
- `i;;`: Represents a span from `i` to the end.
- `;;j`: Represents a span from the beginning to `j`.
- `;;`: Represents a span that includes all elements.
- `i;;j;;k`: Represents a span from `i` through `j` in steps of `k`.
- `i;;;;k`: Represents a span from `i` to the end in steps of `k`.
- `;;j;;k`: Represents a span from the beginning to `j` in steps of `k`.
- `;;;;k`: Represents a span from the beginning to the end in steps of `k`.

**Features**:
- `m[[i;;j;;k]]` is equivalent to `Take[m, {i, j, k}]` but evaluated natively within `Part`.
- `m[[i;;j]] = v` can be used to assign `v` iteratively over a span of elements. If `v` is a list, elements are assigned sequentially. If `v` is a non-list expression, it is assigned uniformly to all elements in the span.
- When used in `Part`, negative `i` and `j` count from the end.
- `i` and `j` can be of the form `UpTo[n]` to restrict endpoints to the actual length of the list.
- Any argument of `Span[...]` can be `All`.

```mathematica
In[1]:= {a, b, c, d, e, f, g, h}[[2;;5]]
Out[1]= {b, c, d, e}

In[2]:= {a, b, c, d, e, f, g, h}[[1;;-1;;3]]
Out[2]= {a, d, g}

In[3]:= t = {a, b, c, d, e, f, g, h}; t[[2;;5]] = x; t
Out[3]= {a, x, x, x, x, f, g, h}

In[4]:= t = {a, b, c, d, e, f, g, h}; t[[2;;5]] = {p, q, r, s}; t
Out[4]= {a, p, q, r, s, f, g, h}

In[5]:= Range[10][[3;;All]]
Out[5]= {3, 4, 5, 6, 7, 8, 9, 10}
```

## Head
Returns the top-level wrapper of an expression.
- `Head[expr]`: gives the head of `expr`.
- `Head[expr, h]`: wraps the result with `h`, i.e. returns `h[Head[expr]]`.

**Features**:
- For functions, returns the symbol or expression acting as the head.
- For atoms, returns the symbolic type name: `Integer`, `Real`, `Rational`, `Complex`, `Symbol`, or `String`.
- The raw objects return their own heads: `SparseArray`, `NDArray` (`List` for
  a packed list), and `CompiledFunction` for the object `Compile` returns
  (v0.6060; before, `Head` stayed unevaluated on it).

```mathematica
In[1]:= Head[f[x]]
Out[1]= f

In[2]:= Head[3/4]
Out[2]= Rational

In[3]:= Head[a + b, f]
Out[3]= f[Plus]
```

## Length
Returns the number of arguments in an expression.
- `Length[expr]`

**Features**: 
- Returns the count of top-level arguments for functions.
- Returns `0` for all atomic expressions.

```mathematica
In[1]:= Length[{a, b, c}]
Out[1]= 3
```

## Dimensions
Returns the list of lengths of levels in a nested rectangular structure.
- `Dimensions[expr]`: gives the dimensions of `expr` at every level where
  it is not ragged.
- `Dimensions[expr, n]`: limits the result to the first `n` levels.

`expr` is treated as a full array only at levels where every sub-piece
shares the same head and length; ragged levels are not counted.
`Dimensions` always returns a `List` (the empty `List {}` for atoms).

```mathematica
In[1]:= Dimensions[{{1, 2}, {3, 4}}]
Out[1]= {2, 2}

In[2]:= Dimensions[{{a, b, c}, {d, e}, {f}}]
Out[2]= {3}

In[3]:= Dimensions[{{{{a, b}}}}]
Out[3]= {1, 1, 1, 2}

In[4]:= Dimensions[{{{{a, b}}}}, 2]
Out[4]= {1, 1}

In[5]:= Dimensions[1]
Out[5]= {}
```

## First, Last, Most, Rest
Convenience functions for accessing parts of a sequence.
- `First[expr]`: Equivalent to `expr[[1]]`.
- `Last[expr]`: Equivalent to `expr[[-1]]`.
- `Most[expr]`: Returns all elements except the last.
- `Rest[expr]`: Returns all elements except the first.

## Reverse
Reverses the order of elements.
- `Reverse[expr]`: Reverses top-level elements.
- `Reverse[expr, n]`: Reverses elements at level `n`.
- `Reverse[expr, {n1, n2, ...}]`: Reverses at levels `n1`, `n2`, etc.
- Numbers are atoms: a level past the array never reaches inside a `Rational`
  or a `Complex`, so `Reverse[{5/4, 7/2}, 2]` is `{5/4, 7/2}` (non-atomic
  expressions of any head are reversed: `Reverse[{f[1, 2]}, 2]` is `{f[2, 1]}`).

## RotateLeft, RotateRight
Cycles elements.
- `RotateLeft[expr, n]`: Cycles `n` positions to the left.
- `RotateLeft[expr]`: Cycles 1 position to the left.
- `RotateLeft[expr, {n1, n2, ...}]`: Cycles at successive levels.
- `RotateRight[...]`: Similar, but cycles to the right.

## PadLeft, PadRight
Pads a (possibly ragged, possibly nested) array to a requested length or set
of dimensions. `PadLeft` pads on the left; `PadRight` is its exact mirror,
padding on the right. The head of `list` need not be `List` (it is preserved).
- `PadLeft[list, n]`: length-`n` list, padded with `0` on the left.
- `PadLeft[list, n, x]`: pads by repeating the element `x`.
- `PadLeft[list, n, {x1, x2, ...}]`: pads by cyclically repeating the `xi`.
- `PadLeft[list, n, padding, m]`: leaves a margin of `m` padding elements on
  the right (`PadRight`: on the left); negative `m` truncates `m` elements.
- `PadLeft[list, {n1, n2, ...}]`: builds a full nested array with length `ni`
  at level `i`, padding ragged data; a nested `padding` block is tiled.
- `PadLeft[list]` / `PadLeft[list, Automatic, x]`: pads a ragged array to a
  full rectangular array (with `0`, resp. `x`).
- A negative length `n` pads on the opposite side.

```
In[1]:= PadLeft[{a, b, c}, 10, {x, y, z}]
Out[1]= {z, x, y, z, x, y, z, a, b, c}

In[2]:= PadRight[{{a, b}, {c}}, {3, 5}]
Out[2]= {{a, b, 0, 0, 0}, {c, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}
```

## Transpose
Transposes levels in a rectangular array.
- `Transpose[list]`: Transposes the first two levels.
- `Transpose[list, {n1, n2, ...}]`: Level `k` becomes level `nk` in result.
  A repeated index (e.g. `{1, 1}`) extracts the corresponding diagonal; a spec
  shorter than the depth leaves the remaining levels in place.

**Features**:
- `Protected`.
- Works only on rectangular arrays.
- A vector is its own transpose: `Transpose[{1, 2, 3}]` and
  `Transpose[{1, 2, 3}, {1}]` are `{1, 2, 3}` (a packed vector stays packed).
- `Transpose[m, {1, 1}]` extracts the diagonal of a square matrix; levels of
  unequal lengths emit `Transpose::diagnl` ("Level rearrangement {1, 1} would
  require collapsing dimensions of unequal lengths.") and stay unevaluated --
  for a `SparseArray` too.
- Spec messages (the call stays unevaluated): `::perm1` (an entry that is not
  a positive machine integer), `::perm2` (an entry beyond the spec's length),
  `::newdims` (a result level with no source, `{2, 2}`), `::tperm` (a
  non-identity spec longer than the depth), `::permspec` (a non-list spec),
  `::nmtx` (`Transpose[{a, {b}}]`: the first two levels are not rectangular).

```mathematica
In[1]:= Transpose[{{a, b}, {c, d}}]
Out[1]= {{a, c}, {b, d}}

In[2]:= Transpose[{{a, b}, {c, d}}, {1, 1}]
Out[2]= {a, d}
```

## ConjugateTranspose
Conjugate transpose (Hermitian transpose) of an array.
- `ConjugateTranspose[m]`: equivalent to `Conjugate[Transpose[m]]`.
- `ConjugateTranspose[m, spec]`: equivalent to `Conjugate[Transpose[m, spec]]`.

**Features**:
- `Protected`.
- On a 1-D vector, `ConjugateTranspose[vec]` conjugates the entries but
  does not change the shape of `vec` (matches the reference convention).
- For symbolic entries, `Conjugate[x]` is left wrapped around `x`.
- Works on higher-rank tensors with the same `spec` semantics as
  `Transpose`.
- Reads and returns a [packed list](packed-arrays.md) at every rank, including
  the rank-1 form. A **real** buffer skips the conjugation entirely — it is the
  identity there, and running it would be a second full pass over the data.

```mathematica
In[1]:= ConjugateTranspose[{{1, 2 I, 3}, {3 + 4 I, 5, I}}]
Out[1]= {{1, 3 - 4 I}, {-2 I, 5}, {3, -I}}

In[2]:= ConjugateTranspose[{{a + b I, c + d I}}]
Out[2]= {{Conjugate[a] - I Conjugate[b]}, {Conjugate[c] - I Conjugate[d]}}

In[3]:= ConjugateTranspose[{1, 2 I, 3 + 4 I}]
Out[3]= {1, -2 I, 3 - 4 I}
```

## Take, Drop

In[1]:= Reverse[{a, b, c}]
Out[1]= {c, b, a}

In[2]:= RotateLeft[{a, b, c}, 1]
Out[2]= {b, c, a}

In[3]:= RotateRight[{a, b, c}, 1]
Out[3]= {c, a, b}
```

## Take, Drop

In[1]:= Most[{a, b, c}]
Out[1]= {a, b}

In[2]:= Rest[f[x, y]]
Out[2]= f[y]
```

## Take, Drop
Extract or remove sequences of elements.
- `Take[expr, spec]`: Extracts elements according to `spec`.
- `Drop[expr, spec]`: Removes elements according to `spec`.

## Flatten
Flattens out nested lists.
- `Flatten[list]`: Flattens all levels.
- `Flatten[list, n]`: Flattens up to level `n`; `n` may be `Infinity`.
- `Flatten[list, n, h]`: Flattens subexpressions with head `h`.
- A negative level emits `Flatten::flev` and stays unevaluated, for a List, a
  packed list and a `SparseArray` alike.

## FlattenAt
Splices the subexpression at one or more positions into its parent, removing
that subexpression's head. Works on any head, not just `List`. Position
resolution (negative indices, single deep path, list of paths) is the shared
`MapAt`/`ReplaceAt` walker; an out-of-range position leaves the call
unevaluated. Attribute: `Protected`.
- `FlattenAt[list, n]`: Flattens the sublist at element `n`; a negative `n`
  counts from the end. `FlattenAt[{a,{b,c},{d,e},{f}}, 2]` → `{a,b,c,{d,e},{f}}`.
- `FlattenAt[expr, {i,j,...}]`: Flattens the part at a single deep position.
- `FlattenAt[expr, {{i1,...},{i2,...},...}]`: Flattens at several positions,
  each resolved against the original expression.
  `FlattenAt[{a,{b,c},{d,e},{f}}, {{2},{4}}]` → `{a,b,c,{d,e},f}`.

Because the head at the position is removed, `FlattenAt[{1,{{2},{3}},4}, {2}]`
gives `{1,{2},{3},4}`; use `MapAt[Flatten, …]` to flatten *within* a part
instead.

## ArrayFlatten
Creates a single flattened matrix from a matrix of matrices (block matrix).
- `ArrayFlatten[a]`: For a matrix of matrices, yields a matrix whose elements are in the same order as in `MatrixForm[a]`. Equivalent to `Flatten[a, {{1,3},{2,4}}]`.
- `ArrayFlatten[a, r]`: Flattens out `r` pairs of levels of a rank-`2r` array, giving a rank-`r` array. Equivalent to `Flatten[a, {{1,r+1},...,{r,2r}}]`.

**Features**:
- `Protected`. Default `r = 2`.
- Blocks must fit: matrices in the same grid row must share their first dimension, and matrices in the same column their second; the output size along an axis is the sum of the block sizes. Disagreeing blocks leave the call unevaluated.
- Elements whose array depth is less than `r` are treated as scalars and replicated to fill a rank-`r` block (e.g. `0` becomes a zero block).

```
In[1]:= m = {{1, 2}, {3, 4}}; ArrayFlatten[{{0, 0, m}, {m, m, 0}}]
Out[1]= {{0, 0, 0, 0, 1, 2}, {0, 0, 0, 0, 3, 4}, {1, 2, 1, 2, 0, 0}, {3, 4, 3, 4, 0, 0}}
```

## Join
Concatenates lists or other expressions that share the same head.
- `Join[list1, list2, ...]`: Concatenates the elements of all lists into a single expression.
- `Join[list1, list2, ..., n]`: Joins the objects at level `n` in each of the lists.

**Features**:
- `Protected`.
- All arguments must share the same head; returns unevaluated if heads differ.
- Works on any head, not just `List` (e.g., `Join[f[a], f[b]]` gives `f[a, b]`).
- `Join[list1, list2, ..., n]` handles ragged arrays by concatenating successive elements at level `n`.

```mathematica
In[1]:= Join[{a, b, c}, {x, y}, {u, v, w}]
Out[1]= {a, b, c, x, y, u, v, w}

In[2]:= Join[{1, 2}, {3, 4}]
Out[2]= {1, 2, 3, 4}

In[3]:= Join[f[a, b], f[c, d]]
Out[3]= f[a, b, c, d]

In[4]:= Join[{{a, b}, {c, d}}, {{1, 2}, {3, 4}}, 2]
Out[4]= {{a, b, 1, 2}, {c, d, 3, 4}}

In[5]:= Join[{{1}, {5, 6}}, {{2, 3}, {7}}, {{4}, {8}}, 2]
Out[5]= {{1, 2, 3, 4}, {5, 6, 7, 8}}

In[6]:= Join[{{x}}, {{1, 2}, {3, 4}}, 2]
Out[6]= {{x, 1, 2}, {3, 4}}
```

## Catenate
`Catenate[{e1, e2, ...}]` takes a single list and flattens one level: the `ei`
must share a head, and their elements are concatenated under it. A list of
associations merges into one association (later keys win) — so `Catenate`
composes with `GroupBy` / `Merge` pipelines that yield a list of associations.

**Features**:
- `Protected`.
- All elements must share a head; returns unevaluated for mixed heads.
- `Catenate[{}]` gives `{}`.

```mathematica
In[1]:= Catenate[{{1, 2}, {3, 4}, {5}}]
Out[1]= {1, 2, 3, 4, 5}

In[2]:= Catenate[{<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "c" -> 4|>}]
Out[2]= <|"a" -> 1, "b" -> 3, "c" -> 4|>

In[3]:= Catenate[Values[GroupBy[Range[6], EvenQ]]]
Out[3]= {1, 3, 5, 2, 4, 6}
```

## Partition
Partitions a list into sublists.
- `Partition[list, n]`: Non-overlapping sublists of length `n`.
- `Partition[list, n, d]`: Sublists with offset `d`.
- `Partition[list, {n1, n2, ...}]`: Multi-level partitioning.
- `Partition[list, spec, dspec]`: Multi-level partitioning with offsets.
- `Partition[list, UpTo[n]]`: Allows shorter final sublist.

**Features**:
- `Protected`.
- Works on any expression with arguments.
- `Partition[list, n, d]` only includes full sublists of length `n` unless `UpTo` is used.

```mathematica
In[1]:= Partition[{a, b, c, d, e}, 2]
Out[1]= {{a, b}, {c, d}}

In[2]:= Partition[{a, b, c, d, e}, 2, 1]
Out[2]= {{a, b}, {b, c}, {c, d}, {d, e}}

In[3]:= Partition[{a, b, c, d, e}, UpTo[2]]
Out[3]= {{a, b}, {c, d}, {e}}

In[4]:= Partition[{{1, 2, 3}, {4, 5, 6}}, {2, 2}]
Out[4]= {{{{1, 2}, {4, 5}}}}
```

## Split
Splits a list into sublists of identical adjacent elements.
- `Split[list]`
- `Split[list, test]`

**Features**:
- `Protected`.
- Uses `SameQ` as the default test.
- Result has the same head as the input.
- Only *adjacent* elements are grouped. To collect equal elements from anywhere
  in the list, use [`Gather`](data-structures.md#gather):
  `Split[{a, b, a}]` gives `{{a}, {b}, {a}}` where `Gather` gives `{{a, a}, {b}}`.

```mathematica
In[1]:= Split[{a, a, a, b, b, a, a, c}]
Out[1]= {{a, a, a}, {b, b}, {a, a}, {c}}

In[2]:= Split[{1, 2, 3, 4, 3, 2, 1}, Less]
Out[2]= {{1, 2, 3, 4}, {3}, {2}, {1}}
```

## Order
- `Order[e1, e2]`: Gives `1` if `e1` is before `e2` in canonical order, `-1` if `e1` is after `e2`, and `0` if `e1` is identical to `e2`.

**Features**:
- `Protected`.
- Uses the same internal canonical comparison (`expr_compare`) as `Sort` and `OrderedQ` — see the canonical-order rules under `Sort` below.
- Compares **structurally**, not by numerical value: `Order[6, Pi]` is `1` (the Integer `6` sorts before the symbol `Pi`), whereas `Order[6, N[Pi]]` is `-1` (two numeric atoms, compared by value).
- Compilable inside `Compile[]` and auto-compiled by `Plot`/`Table`/`NIntegrate`: over machine numbers it lowers to `Sign[e2 - e1]`, returning the Integer `{1, 0, -1}` (matching the interpreter's head). Complex/array arguments fall back to the interpreter.
- Requires exactly two arguments; otherwise it stays unevaluated.

```mathematica
In[1]:= {Order[a, a], Order[a, b], Order[b, a]}
Out[1]= {0, 1, -1}

In[2]:= {Order[6, Pi], Order[6, N[Pi]]}
Out[2]= {1, -1}

In[3]:= Order @@@ Tuples[{0, 1, 2}, 2]
Out[3]= {0, 1, 1, -1, 0, 1, -1, -1, 0}
```

## OrderedQ
- `OrderedQ[expr]`: Gives `True` if the elements of `expr` are in canonical order, and `False` otherwise.
- `OrderedQ[expr, p]`: Uses the ordering function `p` to determine whether each pair of elements is in order.

**Features**:
- `Protected`.
- Uses the same internal canonical comparison logic as `Sort` by default.
- Custom ordering function `p` may return `1`, `0`, `-1`, `True`, or `False`.
- `OrderedQ` works with any expression head, not just `List`.
- Automatically handles 0- and 1-element lists.

```mathematica
In[1]:= OrderedQ[{1, 4, 2}]
Out[1]= False

In[2]:= OrderedQ[{"cat", "catfish", "fish"}]
Out[2]= True

In[3]:= OrderedQ[{1, Sqrt[2], 2, E, 3, Pi}, Less]
Out[3]= True

In[4]:= OrderedQ[{{a, 2}, {c, 1}, {d, 3}}, #1[[2]] < #2[[2]] &]
Out[4]= False

In[5]:= OrderedQ[f[b, a, c]]
Out[5]= False
```

## Sort
Sorts elements of an expression into canonical order.
- `Sort[list]`
- `Sort[list, p]`

**Features**:
- `Protected`.
- Uses an efficient quicksort algorithm.
- Canonical order (v0.597; each rule measured against the reference):
    - **Numbers first, by value.** A complex number orders by its real part,
      then by the size of its imaginary part, negative first:
      `-1 < -1 - I < 0 < -I < I < 1 < 1 - I < 1 + I < 1 - 2 I`. Equal values
      put Integer before Real before Rational (`1 < 1.`, `0.5 < 1/2`).
    - **Products** compare their non-numeric factors from the *last* one
      backwards (any other expression counts as one factor), the shorter list
      winning a tie, and then their numeric coefficients:
      `-2 x < -x < x < 2 x < x^2`, `b < a x`, `y < x y < y^2 < x z`.
    - **Powers** compare bases, then exponents, and any other expression
      counts as `base^1`: `x^-1 < Sqrt[x] < x < x^2 < x^n < y`, and
      `E^x < x` since `E < x`.
    - **Sums** compare their terms the same way as products, with the
      numeric term compared last: `a < 2 + a < a + b < 1 + a + b < a + c`.
      A power is compared before a sum (`x^2 < a + x`, against `(a + x)^1`),
      and a product before both (`a b < a + b`).
    - **Atoms come before compound expressions**, and strings before
      symbols. Symbols use the string collation of their names below, so
      `a < A < b < B` and `x < x1 < xa < x$`.
    - **Other expressions** compare their heads, then their lengths, then
      their arguments from left to right: `f[b] < f[a, b] < g[a] < f[x][a]`,
      and `{b} < {a, b}`.
    - `Infinity` and `ComplexInfinity` order as the `DirectedInfinity[1]` and
      `DirectedInfinity[]` they stand for.
- The same comparison orders the arguments of `Orderless` heads, so
  `Plus` and `Times` print in the reference's order: `a b C d`, `x^-2 + x^-1 + Sqrt[x]`.
- Custom ordering function `p` can return `1`, `0`, `-1`, `True`, or `False`.

**String order.** This is the same collation as WL's, checked against its sort
of 5,233 random printable-ASCII strings. The three levels:

1. **Primary characters only.** Digits come before letters, letters are
   compared case-insensitively, then non-ASCII characters (by code point),
   then `$`. Punctuation and space are ignored at this level.
2. **Where the punctuation falls.** A string with more leading punctuation
   sorts first (`"!-0" < "!0" < "0"`). Among punctuation-only strings, the
   shorter comes first.
3. **Position by position.** Punctuation is compared by its rank in
   ``! @ # % ^ & * ( ) _ - + = ~ ` { [ } ] | \ < > . , ; " ? ' / :`` and then
   space. A lowercase letter comes before its uppercase form.

The result: `"a" < "A" < "a!" < "a1"`,
`"a_b" < "a-b" < "a.b" < "ab" < "aB" < "Ab" < "ab."`, and `".hid"` sorts
beside `"hid"`, after `"a.m"`. Only identical strings tie. Symbols use the
same collation on their names, without the context (v0.597).

```mathematica
In[1]:= Sort[{d, b, c, a}]
Out[1]= {a, b, c, d}

In[2]:= Sort[{Pi, E, 2, 3, 1, Sqrt[2]}, Less]
Out[2]= {1, Sqrt[2], 2, E, 3, Pi}
```

## Ordering
Gives the permutation of positions that sorts an expression.
- `Ordering[list]`: The positions in `list` at which each successive element of `Sort[list]` appears, so that `list[[Ordering[list]]]` is `Sort[list]`.
- `Ordering[list, n]`: The positions of the `n` smallest elements; `Ordering[list, -n]` gives the positions of the `n` largest.
- `Ordering[list, seq]`: Equivalent to `Take[Ordering[list], seq]`; `seq` may be an integer `n` or `-n`, a `{m, n}` or `{m, n, s}` span, `UpTo[k]`, or `All`.
- `Ordering[list, seq, p]` / `Ordering[list, All, p]`: Orders using the ordering function `p`, as in `Sort[list, p]`.

**Features**:
- `Protected`.
- Uses the same internal canonical comparison (`expr_compare`) as `Sort`, and the same custom-ordering-function convention (`p` may return `1`, `0`, `-1`, `True`, or `False`).
- **Stable**: ties are broken by original position, so `Ordering[list, 1]` gives the position of the *first* minimum and `Ordering[{2, 2, 1}]` is `{3, 1, 2}`.
- The result is always a `List` of integer positions, regardless of `list`'s head — `Ordering[f[3, 1, 2]]` is `{2, 3, 1}`.
- Over an `Association`, orders by the **values** and returns their positions.
- Packed-array fast path: on a machine-number vector it argsorts the buffer directly (int64 argsort past `2^53` is exact), returning a packed int64 permutation.
- Compilable inside `Compile[]` and auto-compiled: `Ordering[vector]` lowers to a delegated buffer argsort whose result element type is always integer. A complex dtype, rank ≥ 2, or a custom comparator fall back to the interpreter.

```mathematica
In[1]:= Ordering[{c, a, b}]
Out[1]= {2, 3, 1}

In[2]:= Ordering[{2, 6, 1, 9, 1, 2, 3}, 4]
Out[2]= {3, 5, 1, 6}

In[3]:= Ordering[{2, 6, 1, 9, 1, 2, 3}, -1]
Out[3]= {4}

In[4]:= Ordering[{2, 6, 1, 9, 1, 2, 3}, All, Greater]
Out[4]= {4, 2, 7, 6, 1, 5, 3}

In[5]:= Ordering[<|1 -> c, 2 -> a, 3 -> b|>]
Out[5]= {2, 3, 1}
```

## Union
Gives a sorted list of all distinct elements from one or more expressions.
- `Union[list]`
- `Union[list1, list2, ...]`
- `Union[..., SameTest -> test]`

**Features**:
- `Flat`, `OneIdentity`, `Protected`.
- Without `SameTest`, duplicates are removed by hashing first and only the
  distinct elements are sorted (v0.6052). This is O(total + unique log unique)
  rather than a sort of every element.
- All expressions must have the same head.
- Result has the same head as the inputs.

```mathematica
In[1]:= Union[{1, 2, 1, 3, 6, 2, 2}]
Out[1]= {1, 2, 3, 6}

In[2]:= Union[{a, b, a, c}, {d, a, e, b}, {c, a}]
Out[2]= {a, b, c, d, e}
```

`Union`, `Intersection`, `Complement` and `DeleteDuplicates` have a machine
fast path over a rank-1 buffer of exact integers — the domain where a set
operation is most often a graph traversal or a k-mer count. It is reached from
**either** array representation, the invisible packed `List` and an explicit
`NDArray[...]`, and the result keeps whichever it was given. Until 2026-08-01
only the packed form reached it, so the same call on the same values ran 145×
slower when the argument was written as an `NDArray`.

Reals are deliberately excluded: `0.` and `-0.` compare equal and print
differently, so which of two equal elements survived would depend on the
representation. A custom `SameTest`, a non-`List` head, or any other element
type takes the general path and gives the same answer.

## Intersection
Gives a sorted list of the elements common to all of the inputs.
- `Intersection[list]`
- `Intersection[list1, list2, ...]`
- `Intersection[..., SameTest -> test]`

**Features**:
- `Flat`, `OneIdentity`, `Protected`.
- All expressions must have the same head, which need not be `List`.
- Result has the same head as the inputs; the empty intersection is `{}`.
- With `SameTest -> f`, elements `a`, `b` are treated as equal when `f[a, b]`
  is `True`; the canonically-greatest member of each class is kept.

```mathematica
In[1]:= Intersection[{1, 1, 2, 3}, {3, 1, 4}, {4, 1, 3, 3}]
Out[1]= {1, 3}

In[2]:= Intersection[f[a, b], f[c, a], f[b, b, a]]
Out[2]= f[a]

In[3]:= Intersection[Divisors[45], Divisors[78]]
Out[3]= {1, 3}

In[4]:= Intersection[{1.1, 3.4, .5, 7.6, 7.1, 1.9}, {1.2, 3.3, 7.7, 1.3},
          SameTest -> (Floor[#1] == Floor[#2] &)]
Out[4]= {1.9, 3.4, 7.6}
```

## Complement
Gives the sorted list of the distinct elements of the first argument that are
not in any of the others (set difference).
- `Complement[eall, e1, e2, ...]`
- `Complement[eall, ..., SameTest -> test]`

**Features**:
- `Protected`. Unlike `Union`/`Intersection`, `Complement` is order-sensitive in
  its first argument, so it is *not* `Flat`/`OneIdentity`.
- All expressions must have the same head, which need not be `List`.
- Result has the same head as the inputs; deduplicated and sorted into standard
  order. If nothing survives the removals the result is `{}`.
- Default option `SameTest -> Automatic` (`Options[Complement]`). With
  `SameTest -> f`, elements `a`, `b` are treated as equal when `f[a, b]` is
  `True`; the canonically-smallest member of each class is kept.

```mathematica
In[1]:= Complement[{a, b, c, d, e}, {a, c}, {d}]
Out[1]= {b, e}

In[2]:= Complement[f[a, b, c, d], f[c, a], f[b, b, a]]
Out[2]= f[d]

In[3]:= Complement[{b, e, d, a, b, c, d}, {b, c}]
Out[3]= {a, d, e}

In[4]:= Complement[{1.1, 3.4, .5, 7.6, 7.1, 1.9}, {1.2, 3.3, 1.3},
          SameTest -> (Floor[#1] == Floor[#2] &)]
Out[4]= {0.5, 7.1}
```

## DeleteDuplicates
Removes duplicate elements while preserving order.
- `DeleteDuplicates[list]`
- `DeleteDuplicates[list, test]`

**Features**:
- `Protected`.
- Preserves the order of first occurrences.

```mathematica
In[1]:= DeleteDuplicates[{a, a, b, a, c, b, a}]
Out[1]= {a, b, c}
```

## IntersectingQ
`IntersectingQ[a, b]` gives `True` if `a` and `b` have at least one element in
common, and `False` otherwise.
- Elements are compared with `SameQ`, so `1` and `1.` differ.
- `SameTest -> f` counts `x` and `y` as common when `f[x, y] === True`.
- `a` and `b` may have any head, but it must be the same head for both. A
  mismatch emits `IntersectingQ::heads`; an atom emits `IntersectingQ::lsa`.
- `Protected`. Options: `SameTest -> Automatic`.
- Performance: the default test hashes `b` and probes `a`. Two packed int64
  lists intersect on their machine words, and any other packed or `NDArray`
  shape is delisted. `SameTest` is O(|a| |b|).
- `Compile`: `IntersectingQ[v, w]` over two `_Integer` or `_Real` vectors
  lowers to a Boolean, delegated to this builtin.

## DuplicateFreeQ
`DuplicateFreeQ[list]` gives `True` if no two elements of `list` are the same
(`SameQ`), and `False` otherwise.
- `DuplicateFreeQ[list, test]` counts `x` and `y` as the same when
  `test[x, y]` gives `True`.
- Any head is accepted. An atom emits `DuplicateFreeQ::normal`.
- `Protected`.
- Performance: one hash pass that stops at the first repeat. A packed int64
  list is checked on its machine words.
- `Compile`: `DuplicateFreeQ[v]` over an `_Integer` or `_Real` vector lowers to
  a Boolean.

## Tally
Counts occurrences of elements.
- `Tally[list]`
- `Tally[list, test]`

**Features**:
- `Protected`.
- Returns a list of `{element, count}` pairs.
- Elements appear in the order of their first occurrence.

```mathematica
In[1]:= Tally[{a, a, b, a, c, b, a}]
Out[1]= {{a, 4}, {b, 2}, {c, 1}}
```

## Commonest
Gives a list of the elements that are the most common in an expression.
- `Commonest[list]`
- `Commonest[list, n]`
- `Commonest[list, UpTo[n]]`

**Features**:
- `Protected`.
- When several elements occur with equal frequency, `Commonest` picks first the ones that occur first in `list`.
- `Commonest[list, n]` returns the `n` commonest elements in the order they appear in `list`.
- `Commonest[list, UpTo[n]]` returns the `n` commonest elements, or as many as are available.
- A message `Commonest::dstlms` is generated if there are fewer distinct elements than requested by an integer `n`.
- `n` zero or negative selects nothing.
- **Packed / NDArray fast path.** A rank-1 `int64` or `float64` buffer is
  counted on the machine words, sharing `Tally`'s counting routine so the two
  agree on how a count tie breaks. The selected elements come out of the buffer,
  so the result keeps the input's dtype and stays packed. Everything else (rank
  ≥ 2, a complex dtype, a count that is neither an `Integer` nor
  `UpTo[Integer]`) is answered by the ordinary list code.

```mathematica
In[1]:= Commonest[{b, a, c, 2, a, b, 1, 2}]
Out[1]= {b, a, 2}

In[2]:= Commonest[{b, a, c, 2, a, b, 1, 2}, 4]
Out[2]= {b, a, c, 2}

In[3]:= Commonest[{b, a, c, 2, a, b, 1, 2}, UpTo[6]]
Out[3]= {b, a, c, 2, 1}

In[4]:= Commonest[{1, 2, 2, 3, 3, 3, 4}]
Out[4]= {3}

In[5]:= Commonest[{a, E, Sin[y], E, a, 7}]
Out[5]= {a, E}
```

## Min, Max
Returns the numerically smallest or largest elements.
- `Min[x1, x2, ...]`
- `Max[x1, x2, ...]`
- `Min[{x1, x2, ...}, {y1, ...}, ...]`

**Features**:
- `Flat`, `NumericFunction`, `OneIdentity`, `Orderless`, `Protected`.
- Flattens `List` arguments.
- `Min[]` returns `Infinity`.
- `Max[]` returns `-Infinity`.
- Handles `Infinity` and `-Infinity`.
- Simplifies numeric arguments to a single value; a single element reduces to
  itself (`Min[5]` is `5`, `Max[x]` is `x`).
- Over an association, uses the values (`Min[<|"a" -> 3, "b" -> 1|>]` is `1`).
- Real numeric quantities order with the numbers and keep their exact form
  (v0.477): `Max[1, Pi] = Pi`, `Max[1., Pi] = Pi`, `Min[1/2, Sqrt[2]] = 1/2`,
  `Max[1., Pi, x] = Max[Pi, x]`. A complex quantity has no order and stays.
- Packed arrays and `NDArray`s flatten like Lists, also among other
  arguments: `Max[NDArray[{1., 5.}], 2]` is `5.`, `Max[Range[10], {20}]` is
  `20`. Each array is reduced over its buffer first, so `Max[v, 2]` over 10^6
  elements costs one buffer scan (~0.4 ms), not one Expr per element. An empty
  array contributes nothing: `Max[NDArray[{}]]` is `-Infinity`, as `Max[{}]`
  is.

```mathematica
In[1]:= Min[9, 2]
Out[1]= 2

In[2]:= Min[{4, 1, 7, 2}]
Out[2]= 1

In[3]:= Max[Infinity, 5]
Out[3]= Infinity
```

## RankedMin, RankedMax
Selects the *n*-th smallest or largest element of a list — the order statistics
between `Min`/`Max`.
- `RankedMin[list, n]` — the *n*-th **smallest** element; `RankedMin[list, -n]`
  the *n*-th **largest**.
- `RankedMax[list, n]` — the *n*-th **largest** element; `RankedMax[list, -n]`
  the *n*-th **smallest**.

**Features**:
- `Protected`.
- `RankedMax[list, k]` is `RankedMin[list, -k]`.
- `RankedMin[list, 1]` is `Min[list]`; `RankedMin[list, -1]` is `Max[list]`.
- Yields a definite result whenever every element is a real number, including
  symbolic real constants (`Pi`, `E`, `Sqrt[2]`, `Pi + E`), which order by value;
  `Infinity`/`-Infinity` rank as `±∞`. Returns the element in its exact form.
- Exact for arbitrary-precision integers and rationals; a symbolic non-real
  element (a free symbol or a non-real complex), an empty list, or `|n|` out of
  range leaves the call unevaluated.
- Packed-array fast path (int64 exact, real via O(*n*) quickselect) and a
  `Compile[]` lowering, so `RankedMin[v, k]`/`RankedMax[v, k]` compile and
  auto-compile.

```mathematica
In[1]:= RankedMin[{12, 13, 11}, 2]
Out[1]= 12

In[2]:= RankedMin[{Pi, Sqrt[2], E, 3}, 3]
Out[2]= 3

In[3]:= RankedMax[{2.5, E, 12, 15, 485}, -2]
Out[3]= E

In[4]:= RankedMax[{Infinity, 5, Infinity, -Infinity}, 2]
Out[4]= Infinity
```

## MinMax
`MinMax[list]` gives `{Min[list], Max[list]}` in one shot — the range of the
data, handy for plot bounds. Over an association it uses the values. Delegates
to `Min` / `Max`, so it inherits their numeric handling.

Packed-array aware, and `int64`-exact, since 2026-07-31: it accepts a buffer and
passes it to `Min` and `Max` intact, each of which takes its own buffer path
(`ndred_min` / `ndred_max`). Before that the head was absent from `pack.c`'s
`AWARE` list, so the transparency gate materialised one `Expr` per element and
both delegates ran the generic List code on the boxed copy — 307 ms on a 10⁶
vector against 483 µs now.

```mathematica
In[1]:= MinMax[{3, 1, 4, 1, 5, 9, 2}]
Out[1]= {1, 9}

In[2]:= MinMax[<|"a" -> 3, "b" -> 1, "c" -> 9|>]
Out[2]= {1, 9}
```

## Append, Prepend

- `Protected`.
- Default head to flatten is `List`.
- $n$ must be a non-negative integer.

```mathematica
In[1]:= Flatten[{{a, b}, {c, {d, e}}}]
Out[1]= {a, b, c, d, e}

In[2]:= Flatten[{{a, b}, {c, {d, e}}}, 1]
Out[2]= {a, b, c, {d, e}}

In[3]:= Flatten[f[f[a], b], -1, f]
Out[3]= f[a, b]
```

## Append, Prepend

- `n`: first `n` elements.
- `-n`: last `n` elements.
- `{m, n}`: elements from `m` to `n`.
- `{m, n, s}`: elements from `m` to `n` with step `s`.
- `UpTo[n]`: up to `n` elements as available.

```mathematica
In[1]:= Take[{a, b, c, d}, 2]
Out[1]= {a, b}

In[2]:= Drop[{a, b, c, d}, -1]
Out[2]= {a, b, c}
```

## Append, Prepend
Adds an element to the end or beginning of an expression.
- `Append[expr, elem]`
- `Prepend[expr, elem]`

```mathematica
In[1]:= Append[{a, b}, c]
Out[1]= {a, b, c}
```

## AppendTo, PrependTo
Updates a variable by appending or prepending an element.
- `AppendTo[symbol, elem]`
- `PrependTo[symbol, elem]`

**Features**:
- `HoldFirst`. The target may be a symbol or any assignable form with a value,
  such as `m[[1]]`, `f[1]`, `Options[f]`, or a `Module` local.
- The new value is `Append[x, elem]` / `Prepend[x, elem]`, so these accept every
  representation `Append` does. A packed array keeps its buffer
  (`a = Range[5]; AppendTo[a, 6]` leaves `a` packed), and an element whose
  type does not fit the buffer gives an ordinary List (`AppendTo[a, x]`).
  Before v0.522 any packed variable was silently left unchanged.
- An association appends by key; a non-rule element assigns the unevaluated
  `Append[assoc, elem]`, as in the reference.
- `X::rvalue` when the symbol has no value; `X::normal` when the value is an
  atom (`b = 5; AppendTo[b, 1]`). Both leave the call unevaluated.

```mathematica
In[1]:= m = Table[i + j, {i, 2}, {j, 2}]; AppendTo[m[[1]], 9]; m
Out[1]= {{2, 3, 9}, {3, 4}}
```

## Insert, Delete
Inserts or removes elements at specified positions.
- `Insert[expr, elem, pos]`
- `Delete[expr, pos]`

**Features**:
- `pos` can be a single index, a list (path), or a list of paths.
- `Delete[expr, 0]` replaces the head with `Sequence`.

```mathematica
In[1]:= Insert[{a, b, c}, x, 2]
Out[1]= {a, x, b, c}

In[2]:= Delete[{a, b, c}, 2]
Out[1]= {a, c}
```

## Depth
Gives the maximum number of indices needed to specify any part of expr, plus 1.
- `Depth[expr]`

**Features**:
- `Protected`.
- Default option: `Heads -> False`.
- Raw objects (atoms) have depth 1.
- Numbers, `Rational`, and `Complex` have depth 1.
- Symbolic constants like `Pi`, `E`, `I` have depth 1.
- Compound expressions have depth `1 + Max(depths of arguments)`; an empty one
  has depth 2 (`Depth[{}]` is 2, as in the reference).
- An association counts only its values: `Depth[<|a -> 1|>]` is 2 (see
  [Association atomicity](data-structures.md#association-atomicity-structural-functions)).
- With `Heads -> True`, it includes heads of expressions and their parts.

```mathematica
In[1]:= Depth[f[g[h[x]]]]
Out[1]= 4

In[2]:= Depth[1/2]
Out[2]= 1

In[3]:= Depth[h[{{{a}}}][x, y]]
Out[3]= 2

In[4]:= Depth[h[{{{a}}}][x, y], Heads -> True]
Out[4]= 6
```

## Level
Gives a list of all subexpressions of expr on levels specified by levelspec.
- `Level[expr, levelspec]`
- `Level[expr, levelspec, f]`

**Features**:
- `Protected`.
- Default option: `Heads -> False`.
- Standard level specifications:
  - `n`: levels 1 through `n`.
  - `Infinity`: levels 1 through `Infinity`.
  - `{n}`: level `n` only.
  - `{n1, n2}`: levels `n1` through `n2`.
- Positive level `n` refers to distance from the top (level 0 is the whole expression).
- Negative level `-n` refers to distance from the bottom (depth `n`).
- Level `-1` corresponds to atomic objects.
- An association's parts are its values: `Level[<|a -> 1, b -> 2|>, {1}]` is
  `{1, 2}` (see [Association atomicity](data-structures.md#association-atomicity-structural-functions)).
- Lists subexpressions in post-order (depth-first), resulting in lexicographic ordering of indices.
- A visible `NDArray` answers as its List does:
  `Level[NDArray[{{1, 2, 3}, {5, 6, 7}}], {2}]` is `NDArray[{1, 2, 3, 5, 6, 7}]`.

```mathematica
In[1]:= Level[a + f[x, y^n], {-1}]
Out[1]= {a, x, y, n}

In[2]:= Level[a + f[x, y^n], 2]
Out[2]= {a, x, y^n, f[x, y^n]}

In[3]:= Level[x^2 + y^3, 3, Heads -> True]
Out[3]= {Plus, Power, x, 2, x^2, Power, y, 3, y^3}
```

## Variables
Gives an ordered list of all independent variables in a polynomial.
- `Variables[poly]`

**Features**:
- `Protected`.
- Looks for variables only inside `Plus`, `Times`, and `Power` with rational exponents.
- Returns a sorted `List` of variables.
- Symbolic constants like `Pi`, `E`, and `I` are not treated as variables.
- **A power with a symbolic exponent splits** as the reference reads it:
  `b^(n + k t)` is `b^n (b^t)^k`, so `Variables[x^(1 + e) + x^(2 + e)]` is
  `{x, x^e}`. `x^(2 e)` and `x^(-e)` give `x^e`, while `x^(e/2)` and
  `x^(a e)` are their own variables. A numeric base is a constant:
  `Variables[2^(1 + x)]` is `{}`. The same generators drive `Factor` and
  `FactorList` (`src/poly/expsplit.c`).

```mathematica
In[1]:= Variables[(x + y)^2 + 3 z^2 - y z + 7]
Out[1]= {x, y, z}

In[2]:= Variables[Sin[x] + Cos[x]]
Out[2]= {Cos[x], Sin[x]}

In[3]:= Variables[E^x]
Out[3]= {}
```

## Expand
Expands out products and positive integer powers in an expression.
- `Expand[expr]`
- `Expand[expr, patt]`

**Features**:
- `Protected`.
- Works only on positive integer powers and distributes products over sums.
- Threads over equations, inequalities, and lists.
- **Polynomials over the rationals are expanded through FLINT** (packed
  `fmpq_mpoly` arithmetic), which distributes and collects like terms orders of
  magnitude faster than the generic tree multiplier on dense, high-degree
  inputs — e.g. `(1 + x + 5 x^3 + 8 x^17)^341` (degree 5797, 5758 terms) expands
  in a few hundredths of a second. Non-polynomial parts (transcendental heads,
  symbolic/fractional exponents, inexact coefficients) fall back to the generic
  binary-splitting distributor and repeated-squaring power expander.
- **Opaque kernels are FLINT generators too** (v0.6054). A polynomial in
  kernels such as `f[1]`, `g[x[i], x[j]]`, `Sin[x]` or `Log[y]` is expanded in
  packed arithmetic, each distinct kernel one variable (found through a hash
  table), and the result is evaluated, so products that `Times` rewrites
  (`Sin[x] Csc[x]` is `1`) come out exactly as term-by-term distribution gives
  them. Speed-ups: `(f[1] + ... + f[6])^8` takes 4.9 ms (was 78 ms), and
  `(x + Sin[x] + g[y] + 1)^12` takes 1.4 ms (was 26 ms).
  - **When it is used.** Only when sums are actually multiplied together (a
    power of a sum, or two or more sum factors). Distributing one monomial over
    an expanded sum stays on the classical path, which keeps canonical order.
  - **What stays classical.** Heads with their own multiplication semantics:
    System heads without `NumericFunction` (`List`, `Interval`, `SeriesData`,
    `ConditionalExpression`, ...), a head with UpValues, user rules for
    `Plus`/`Times`/`Power`, and inputs needing more than 128 generators.
- In a **mixed product** — one where a non-polynomial factor (e.g. `Log[x]`)
  would otherwise force the whole product onto the generic path — the
  polynomial-over-`Q` factors are still multiplied together through FLINT, and
  only the non-polynomial factors are distributed generically. So
  `Expand[Log[x] (1+y)^300 (1-y)^300]` collapses the two degree-300 factors in
  packed arithmetic rather than distributing 90 601 terms by hand.
- `ExpandNumerator` and `ExpandDenominator` inherit the same acceleration (they
  expand their numerator/denominator polynomial through the same path).
- **Never silently declines.** `Expand` expands anything that fits in memory; an
  expansion whose (estimated) size exceeds the memory ceiling returns
  `Overflow[]` rather than leaving the expression unevaluated. The size estimate
  is a Newton-box upper bound — the tighter of the multinomial term count and
  the per-variable degree box — so univariate and low-dimensional expansions of
  *any* degree are recognised as cheap and always run; only a genuinely
  high-dimensional combinatorial blow-up (e.g. `(a+b+c+d+e+f+g)^60`) overflows.
- `Expand[expr, patt]` leaves unexpanded any parts of `expr` that are free of the pattern `patt`. Inside a product, pattern-free factors are carried along as an unexpanded coefficient rather than being distributed.

```mathematica
In[1]:= Expand[(x+3)(x+2)]
Out[1]= 6 + 5 x + x^2

In[2]:= Expand[(x+y)^2 (x-y)^2]
Out[2]= x^4 - 2 x^2 y^2 + y^4

In[3]:= Expand[(x+1)^2 + (y+1)^2, x]
Out[3]= 1 + 2x + x^2 + (1+y)^2

In[4]:= Expand[(a+b)(x[1]+x[2])^2, x[_]]
Out[4]= (a+b) x[1]^2 + 2 (a+b) x[1] x[2] + (a+b) x[2]^2

In[5]:= Length[Expand[(1 + x + 5 x^3 + 8 x^17)^341]]
Out[5]= 5758

In[6]:= Expand[(a+b+c+d+e+f+g)^60]
Out[6]= Overflow[]
```

## ExpandNumerator
Expands out products and powers that appear in the numerator of an expression.
- `ExpandNumerator[expr]`

**Features**:
- `Protected`.
- Acts only on factors with positive integer exponents (the "numerator part" of `expr`).
- Applies only to the top level in `expr`; it does not descend into function bodies.
- Leaves the denominator factors (those with negative integer exponents) unchanged.
- Does not separate the fraction into a sum of fractions; only `Expand` does that.
- Threads over `List`, `Equal`, `Unequal`, `Less`, `LessEqual`, `Greater`, `GreaterEqual`, `And`, `Or`, `Not`, and `Plus` (so each summand of a sum-of-fractions is processed independently).

```mathematica
In[1]:= ExpandNumerator[(x-1)(x-2)/((x-3)(x-4))]
Out[1]= (2 - 3 x + x^2)/((x-3)(x-4))

In[2]:= ExpandNumerator[(a+b)^2/x + (c+d)(c-d)/y]
Out[2]= (a^2 + 2 a b + b^2)/x + (c^2 - d^2)/y

In[3]:= ExpandNumerator[x == (a+b)^2/c && y >= (a-b)^2/c]
Out[3]= x == (a^2 + 2 a b + b^2)/c && y >= (a^2 - 2 a b + b^2)/c
```

## ExpandDenominator
Expands out products and powers that appear as denominators in an expression.
- `ExpandDenominator[expr]`

**Features**:
- `Protected`.
- Acts only on factors with negative integer exponents (the "denominator part" of `expr`).
- Combines all denominator factors of a top-level `Times` into a single expanded polynomial wrapped in `Power[..., -1]`.
- Applies only to the top level in `expr`; it does not descend into function bodies.
- Leaves the numerator factors unchanged.
- Threads over `List`, `Equal`, `Unequal`, `Less`, `LessEqual`, `Greater`, `GreaterEqual`, `And`, `Or`, `Not`, and `Plus`.

```mathematica
In[1]:= ExpandDenominator[(x-1)(x-2)/((x-3)(x-4))]
Out[1]= ((x-1)(x-2))/(12 - 7 x + x^2)

In[2]:= ExpandDenominator[1/(x+1) + 2/(x+1)^2 + 3/(x+1)^3]
Out[2]= 1/(1 + x) + 2/(1 + 2 x + x^2) + 3/(1 + 3 x + 3 x^2 + x^3)

In[3]:= ExpandDenominator[(a+b)(a-b)/((c+d)(c-d))]
Out[3]= ((a+b)(a-b))/(c^2 - d^2)
```

## ExpandAll
Expands out all products and integer powers in *any* part of an expression.
- `ExpandAll[expr]`
- `ExpandAll[expr, patt]`

**Features**:
- `Protected`.
- A thin recursive driver over the accelerated `Expand`: it descends into every
  part of `expr` — function heads, arguments, exponents, and the bases of
  denominators — and applies `Expand` (and `ExpandDenominator`) at each node,
  where a top-level `Expand` reaches none of them.
- `ExpandAll[expr, patt]` avoids expanding parts of `expr` that do not contain
  terms matching the pattern `patt`.
- Held parts are left alone: an argument that its head holds (`HoldForm[x]`,
  `Hold[...]`, the branches of `If`) is not descended into, so a held
  abbreviation is never evaluated (v0.586).
- Threads over `List`, equations, inequalities, and logic functions.
- Called with a number of arguments other than 1 or 2, emits `ExpandAll::argt`
  and stays unevaluated.

```mathematica
In[1]:= ExpandAll[1/(1+x)^3 + Sin[(1+x)^3]]
Out[1]= 1/(1 + 3 x + 3 x^2 + x^3) + Sin[1 + 3 x + 3 x^2 + x^3]

In[2]:= ExpandAll[(x+z)^2/(x+y)^2]
Out[2]= x^2/(x^2 + 2 x y + y^2) + (2 x z)/(x^2 + 2 x y + y^2) + z^2/(x^2 + 2 x y + y^2)

In[3]:= ExpandAll[E^(I a (t-b))]
Out[3]= E^(-I a b + I a t)

In[4]:= ExpandAll[((1+a) (1+b))[x]]
Out[4]= (1 + a + b + a b)[x]

In[5]:= ExpandAll[(f[(x+y)^2] + g[(y+z)^2])^2, x]
Out[5]= f[x^2 + 2 x y + y^2]^2 + 2 f[x^2 + 2 x y + y^2] g[(y+z)^2] + g[(y+z)^2]^2
```

## PowerExpand
Expands powers of products and nested powers, and logarithms/arguments of
products and powers.
- `PowerExpand[expr]`
- `PowerExpand[expr, {x1, x2, ...}]`
- `PowerExpand[expr, Assumptions -> assum]`

**Features**:
- `Protected`.
- Applies `(a b)^c -> a^c b^c`, `(a^b)^c -> a^(b c)`, `Log[a b] -> Log[a] + Log[b]`,
  `Log[a^b] -> b Log[a]`, and `Arg[a b] -> Arg[a] + Arg[b]`. Because `Sqrt[x]` is
  `Power[x, 1/2]` and `Log[1/z]` is `Log[z^-1]`, these also give
  `Sqrt[x y] -> Sqrt[x] Sqrt[y]`, `Sqrt[z^2] -> z`, and `Log[1/z] -> -Log[z]`.
- Rules are applied top-down to a fixed point, so `Log[(a b)^c]` becomes
  `c (Log[a] + Log[b])` rather than expanding the inner power first.
- When a product raised to a non-integer (root) exponent carries a negative
  numeric coefficient together with a `Plus` factor, the sign is folded into
  that factor so the root stays real, e.g.
  `Sqrt[-4 Dt[u]^2 (-1 + u)] -> 2 Dt[u] Sqrt[1 - u]` rather than
  `2 I Dt[u] Sqrt[-1 + u]`. This is value-preserving since `(-1)^c p^c = (-p)^c`.
- `f^-1[f[x]] -> x` for the inverse-trig / inverse-hyperbolic pairs
  (`ArcTan[Tan[x]] -> x`, `ArcSin[Sin[x]] -> x`, …).
- The default `Assumptions -> Automatic` makes the textbook transforms, which
  are correct when the bases are positive reals or the exponents are integers;
  branch cuts are otherwise ignored.
- `Assumptions -> True` emits the universally-correct result, attaching a
  branch-correction term built from `Floor`, `Arg`, `Im`, `E`, `I`, and `Pi`.
- `Assumptions -> assum` emits the `True`-mode formula and then refines the
  correction terms under `assum` (resolving `Arg`/`Im` of known-sign reals and
  evaluating `Floor` over assumption-bounded intervals). This is faithful on
  the documented examples; for assumptions outside this reasoning it degrades
  to the symbolic `True`-mode form.
- `PowerExpand[expr, {x1, ...}]` expands only subexpressions that mention one
  of the listed variables.
- Threads over `List`, equations, inequalities, and logic functions.

Limitations: in `Assumptions -> True` mode the nested-power correction is left
in `Im[...]` form rather than rewritten to `Arg[...]` (mathematically equal).

```mathematica
In[1]:= PowerExpand[Sqrt[x y]]
Out[1]= Sqrt[x] Sqrt[y]

In[2]:= PowerExpand[Log[(a b)^c]]
Out[2]= c (Log[a] + Log[b])

In[3]:= PowerExpand[Sqrt[a b] + Sqrt[c d], {a, b}]
Out[3]= Sqrt[a] Sqrt[b] + Sqrt[c d]

In[4]:= PowerExpand[Sqrt[z^2], Assumptions -> z < 0]
Out[4]= -z

In[5]:= PowerExpand[Log[x y], Assumptions -> True]
Out[5]= Log[x] + Log[y] + 2 I Pi Floor[1/2 - Arg[x]/(2 Pi) - Arg[y]/(2 Pi)]
```

## Exponent
Gives the maximum power with which a form appears in the expanded form of an expression.
- `Exponent[expr, form]`
- `Exponent[expr, form, h]`

**Features**:
- `Listable`, `Protected`.
- The default aggregator is `h = Max`. `Exponent[expr, form, Min]` gives the lowest power; `Exponent[expr, form, List]` gives the sorted, de-duplicated set of exponents.
- `form` may be a symbol, a kernel (e.g. `Sin[x]`), or a product of terms.
- Works whether or not `expr` is explicitly given in expanded form (it expands internally).
- Purely syntactic: it does not attempt to recognise zero coefficients.
- Exponents may be rational numbers or symbolic expressions.
- A monomial given as a non-flat nested `Times` is handled: a term such as
  `Times[-1, Times[a, x^2]]` (which `-a x^2` evaluates to beside a complex
  sibling like `I x`) is recursed into, so `Exponent[-a x^2 + I x, x]` is `2`,
  not `1`.
- `Exponent[0, x]` is `-Infinity` (empty exponent set, `h = Max`).
- The `Listable` attribute makes `Exponent[expr, {form1, form2, ...}]` give the list of exponents for each `formi`.

```mathematica
In[1]:= Exponent[1 + x^2 + a x^3, x]
Out[1]= 3

In[2]:= Exponent[0, x]
Out[2]= -Infinity

In[3]:= Exponent[x^(n+1) + 2 Sqrt[x] + 1, x]
Out[3]= Max[1/2, 1 + n]

In[4]:= Exponent[(x^2+1)^3 - 1, x, Min]
Out[4]= 2

In[5]:= Exponent[1 + x^2 + a x^3, x, List]
Out[5]= {0, 2, 3}
```

## Coefficient
Gives the coefficient of a specific form in a polynomial.
- `Coefficient[expr, form]`
- `Coefficient[expr, form, n]`

**Features**:
- `Protected`, `Listable`.
- `Coefficient[expr, form, 0]` picks out terms that do NOT contain `form`.
- `n` may be negative: `Coefficient[a/x + b x, x, -1]` is `a`, and a term in
  `1/x` is not part of the `x^0` coefficient (`Coefficient[x1/e + x2, e, 0]`
  is `x2`).
- A term whose other factors put `form` in a denominator (`x/(1 + x)` for form
  `x`) has no coefficient of any order; factors that merely contain `form`
  elsewhere count (`Coefficient[x Sin[x] + x, x]` is `1 + Sin[x]`).
- `Coefficient[expr, form, n]` stays unevaluated when `n` is not an explicit
  integer (a symbol, `Real`, or `Rational`).
- Expands only in `form`, as the reference does, after taking out the factors
  free of it: `Coefficient[(a l^2 + b l^2 + c)/d, l, 2]` is `(a + b)/d` and
  `Coefficient[(a + b) x/c + d x, x]` is `(a + b)/c + d` (v0.592).
- Treats distinct transcendental powers as algebraically unrelated (e.g., `x^s` is treated as a separate base from `x`).

```mathematica
In[1]:= Coefficient[(x+1)^3, x, 2]
Out[1]= 3

In[2]:= Coefficient[(x+y)^4, x y^3]
Out[2]= 4

In[3]:= Coefficient[x^s x, x^s]
Out[1]= x
```

## CoefficientList
Gives a list of coefficients of powers of variables in a polynomial.
- `CoefficientList[poly, var]`
- `CoefficientList[poly, {var1, var2, ...}]`

**Features**:
- `Protected`.
- Gives an array of coefficients starting with power 0.
- Returns a full rectangular array for multiple variables. Combinations of powers that do not appear in `poly` give zeros in the array.
- Automatically expands the polynomial internally.

```mathematica
In[1]:= CoefficientList[1 + 6 x - x^4, x]
Out[1]= {1, 6, 0, 0, -1}

In[2]:= CoefficientList[(1 + x)^10, x]
Out[2]= {1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1}

In[3]:= CoefficientList[1 + a x^2 + b x y + c y^2, {x, y}]
Out[3]= {{1, 0, c}, {0, b, 0}, {a, 0, 0}}
```

## MonomialList
Gives the list of monomials of a polynomial, sorted by a monomial order.
- `MonomialList[poly]`
- `MonomialList[poly, {var1, var2, ...}]`
- `MonomialList[poly, vars, order]`

**Features**:
- `Protected`.
- `MonomialList[poly]` is equivalent to `MonomialList[poly, Variables[poly]]`.
- Works whether or not `poly` is given in expanded form (it is expanded internally).
- `order` is `"Lexicographic"` (default), `"DegreeLexicographic"`,
  `"DegreeReverseLexicographic"`, `"NegativeLexicographic"`,
  `"NegativeDegreeLexicographic"`, `"NegativeDegreeReverseLexicographic"`, or an
  explicit weight matrix `w` (monomials ranked by the lexicographic order of `w.v`).
- `"NegativeLexicographic"` is `Sort` of the exponent vectors; `"Lexicographic"` is
  its reverse.
- `Modulus -> m` reduces coefficients modulo `m` (dropping monomials whose
  coefficient vanishes).
- `vars` may be `All` (equivalent to `Variables[poly]`).
- `Plus @@ MonomialList[poly, vars]` reconstructs the expanded polynomial.

```mathematica
In[1]:= MonomialList[(x + y)^3]
Out[1]= {x^3, 3 x^2 y, 3 x y^2, y^3}

In[2]:= MonomialList[x^2 y^2 + x^3, {x, y}, "DegreeLexicographic"]
Out[2]= {x^2 y^2, x^3}

In[3]:= MonomialList[(x + 1)^5, x, Modulus -> 2]
Out[3]= {x^5, x^4, x, 1}
```

## CoefficientRules
Gives the exponent vectors and coefficients of the monomials of a polynomial as rules.
- `CoefficientRules[poly]`
- `CoefficientRules[poly, {var1, var2, ...}]`
- `CoefficientRules[poly, vars, order]`

**Features**:
- `Protected`.
- Returns `{expvec -> coeff, ...}`, one rule per monomial; the exponent vector lists
  the powers of `vars` in order.
- `CoefficientRules[poly]` is equivalent to `CoefficientRules[poly, Variables[poly]]`.
- Same `order` settings and `Modulus -> m` option as `MonomialList`; `vars` may be
  `All`. Works whether or not `poly` is expanded.
- `FromCoefficientRules` is the inverse.
- Note: the no-variable form uses `Variables[poly]`, which Mathilda returns in
  canonical (sorted) order — so e.g. `CoefficientRules[y + x z]` uses the variable
  order `{x, y, z}`.
- Polynomials whose coefficients live in one number field `Q(θ)`
  (`AlgebraicNumber[θ, {..}]`, all sharing one θ — the ParallelMixedTower assembly's
  representation) take a native FLINT read-off (θ→a fresh variable, group over the
  field, read each coefficient back mod θ's minimal polynomial), byte-identical to
  the generic path but without the per-term evaluator work. `MonomialList` shares it.

```mathematica
In[1]:= CoefficientRules[(x + y)^3]
Out[1]= {{3, 0} -> 1, {2, 1} -> 3, {1, 2} -> 3, {0, 3} -> 1}

In[2]:= CoefficientRules[a x y^2 + b x^2 z, {x, y, z}, "DegreeReverseLexicographic"]
Out[2]= {{1, 2, 0} -> a, {2, 0, 1} -> b}

In[3]:= CoefficientRules[(x + 1)^5, x, Modulus -> 2]
Out[3]= {{5} -> 1, {4} -> 1, {1} -> 1, {0} -> 1}
```

## CoefficientArrays
Gives the coefficient arrays of a polynomial or polynomial system, grouped by total degree.
- `CoefficientArrays[polys, vars]`
- `CoefficientArrays[poly, vars]` (a non-List first argument is one polynomial)

**Features**:
- `Protected`. `vars` may be `All` or omitted (then `Variables` of the whole input).
- Element `d + 1` holds every degree-`d` coefficient at the index tuple naming its
  variables in non-decreasing order (upper-triangular, not symmetric); a system
  carries the equation index as the leading axis, so `{b, m} = CoefficientArrays[eqs, vars]`
  is the constant vector and the coefficient matrix of a linear system.
- **Every array element is a `SparseArray`** with background `0`;
  the constant term of a single polynomial is a scalar. Storage follows the number of
  terms, not `k^d`. `Normal` gives the dense arrays; `Part`, `Dot`, `Map`, `Transpose`
  and the linear-algebra heads accept the result (v0.586).
- An element of rank above 32 (degree 33 or more) has no SparseArray and is a dense
  nested List; a dense total over 2·10⁷ entries is declined with `CoefficientArrays::dense`.
- `Modulus -> p` reduces the coefficients (a Mathilda extension).

```mathematica
In[1]:= CoefficientArrays[{x + 2 y - 1, 3 x - y}, {x, y}]
Out[1]= {SparseArray[<1>, {2}], SparseArray[<4>, {2, 2}]}

In[2]:= Normal[%]
Out[2]= {{-1, 0}, {{1, 2}, {3, -1}}}

In[3]:= CoefficientArrays[x y + c y^2 + a x, {x, y}] // Normal
Out[3]= {0, {a, 0}, {{0, 1}, {0, c}}}
```

## GroebnerBasis`DistributedTermsList
The distributed form of a polynomial: its terms as exponent vectors with coefficients.
- ``GroebnerBasis`DistributedTermsList[poly, vars]``
- ``GroebnerBasis`DistributedTermsList[{poly1, ...}, vars]``

**Features**:
- `Protected`. An internal helper that symbolic packages call directly.
- Returns `{{{e1, ..., ek}, c}, ...}, vars}` in descending monomial order; a list of
  polynomials gives one term list per polynomial. `{{}, vars}` for `0`.
- Option `MonomialOrder -> o` (default `Lexicographic`); `o` may be a string or a
  symbol, as in `CoefficientRules` (which now also accepts a symbol order).
- A generalised variable that would have to be adjoined (`y^(1/2)` in `vars` `{x, y}`)
  is declined.

```mathematica
In[1]:= GroebnerBasis`DistributedTermsList[a x^2 + 3 x y + 5 y^3, {x, y}]
Out[1]= {{{{2, 0}, a}, {{1, 1}, 3}, {{0, 3}, 5}}, {x, y}}
```

## FromCoefficientRules
Reconstructs a polynomial from its `CoefficientRules` representation.
- `FromCoefficientRules[rules, {var1, var2, ...}]`

**Features**:
- `Protected`.
- Builds `Sum[coeff * var1^e1 * var2^e2 * ..., over rules]`; the exponent vectors must
  match the number of variables.
- Inverse of `CoefficientRules`: `FromCoefficientRules[CoefficientRules[poly, vars], vars]`
  is `Expand[poly]`.

```mathematica
In[1]:= FromCoefficientRules[{{2, 0} -> a, {1, 1} -> b, {0, 2} -> c}, {x, y}]
Out[1]= a x^2 + b x y + c y^2
```

## Collect
Collects together terms involving the same powers of objects matching a variable or variables.
- `Collect[expr, x]`
- `Collect[expr, {x1, x2, ...}]`
- `Collect[expr, var, h]`

**Features**:
- `Protected`.
- Automatically threads over lists, equations, inequalities, and logic functions.
- Effectively writes `expr` as a polynomial in `x` or a fractional power of `x`.
- `Collect[expr, var, h]` applies `h` to the expression that forms the coefficient of each term obtained.
- `var` may be a pattern (`Power[_, 1/2]`, `f[_]`, `x^_`, `p | q`, `f[n_] /; n > 0`). Every distinct subexpression of `expr` that matches it, at levels 1 to Infinity, becomes a collection variable of its own, in canonical order (v0.516).
- Numeric coefficients are never powers of the key: `Collect[2 x + 2 y, Sqrt[2]]` stays `2 x + 2 y` (v0.516).

```mathematica
In[1]:= Collect[a x + b y + c x + d y, x]
Out[1]= b y + d y + (a + c) x

In[2]:= Collect[a x^4 + b x^4 + 2 a^2 x - 3 b x + x - 7, x]
Out[2]= -7 + (1 + 2 a^2 - 3 b) x + (a + b) x^4

In[3]:= Collect[a Sqrt[x] + Sqrt[x] + x^(2/3) - c x + 3x - 2b x^(2/3) + 5, x]
Out[3]= 5 + (1 + a) Sqrt[x] + (1 - 2b) x^(2/3) + (3 - c) x

In[4]:= Collect[1 + 3 Sqrt[1 + x^2] + x Sqrt[1 + x^2], Power[_, 1/2]]
Out[4]= 1 + (3 + x) Sqrt[1 + x^2]

In[5]:= Collect[a f[1] + b f[1] + c f[2] + d f[2], f[_]]
Out[5]= (a + b) f[1] + (c + d) f[2]
```

## FactorSquareFree
Pulls out any multiple factors in a polynomial.
- `FactorSquareFree[poly]`
- `FactorSquareFree[poly, Extension -> None]` -- factor over Q (the default).
- `FactorSquareFree[poly, Extension -> alpha]` -- decompose over Q(alpha) (`alpha` = `Sqrt[c]`, `c^(1/n)`, `I`, or a list of generators).
- `FactorSquareFree[poly, Modulus -> p]` -- decompose over GF(p).

**Features**:
- `Listable`, `Protected`.  `Options[FactorSquareFree]` is `{Extension -> None, Modulus -> 0, Trig -> False}`; any other option is `FactorSquareFree::optx` (in particular `GaussianIntegers`).
- **FLINT** (v0.6061): a polynomial over Q, and each side of a rational function, decomposes through `fmpq_mpoly_factor_squarefree`, 65-340x faster than the classical Yun loop on expression trees (a 4-variable, 10920-term input: 102.7 s -> 0.57 s).  The numerical content stands in front, monomial factors are split out and each factor has a positive leading coefficient: `FactorSquareFree[x + 7 x^2 + 12 x^3]` is `x (1 + 7 x + 12 x^2)`, `FactorSquareFree[(x^2 - 1)^2/(x + 1)^3]` is `(-1 + x)^2/(1 + x)`.  The classical path below remains for builds without FLINT and for inputs FLINT declines.
- `Modulus -> p` decomposes over GF(p) (`fmpz_mod_mpoly_factor_squarefree`): `FactorSquareFree[x^3 + x^2, Modulus -> 2]` is `x^2 (1 + x)`.
- Automatically threads over lists, as well as equations, inequalities and logic functions.
- Works on both univariate and multivariate polynomials.
- Bigint coefficients (above 2^63) are coefficients, not variables: `FactorSquareFree[Expand[3^45 (1 + t)^2]]` gives `2954312706550833698643 (1 + t)^2` (v0.458; it used to come back unfactored).
- Self-checking: if the computed content and square-free factors do not reproduce the input exactly, which can happen when the "variables" are transcendental atoms over a coefficient ring the classical PRS does not model, the input is returned unfactored, never as an unequal product (v0.458).
- `Extension` runs the decomposition over an algebraic number field Q(alpha).  It reduces to `Factor[poly, Extension -> alpha]` (Trager engine) and groups the irreducible-over-Q(alpha) factors by multiplicity, so it also accepts inputs with algebraic-number coefficients, e.g. `FactorSquareFree[x^2 - 2 Sqrt[2] x + 2, Extension -> Sqrt[2]]` → `(-Sqrt[2] + x)^2`.  Because square-free structure is field-independent for rational inputs, an extension leaves a rational-coefficient answer unchanged; the default `Extension -> None` uses the fast over-Q Yun path.
- Multivariate inputs use a cheap squarefree pre-check (F4 Stage 1): after content extraction in the main variable, `sqfree_cheap_check` substitutes integer values from `{1, -1, 2, -2, 3, -3, 4}` for the other variables and tests `gcd(image, image')` over `Z[x]`.  If any image is squarefree at an alpha that preserves the leading-x degree, the pre-check proves `pp` is squarefree in x and the expensive multivariate `gcd(pp, pp')` is skipped.  Soundness comes from content extraction guaranteeing any repeated factor of `pp` involves the main variable nontrivially.  Measured 6.6× speedup on 4-variable squarefree inputs (6.27 s → 0.95 s); non-squarefree inputs fall through to the original Yun loop with negligible overhead.
- The cheap pre-check's univariate `gcd(image, image')` runs through `zupoly_gcd` (subresultant PRS, GMP `mpz_t` coefficients).  The previous implementation used `poly_gcd_internal` (Knuth-style primitive PRS at the Expr level) which suffers exponential coefficient growth on the intermediate pseudo-remainders; on a degree-31 univariate image (e.g. `Factor[Expand[x^2 (z^13 - x^12)(z^4 + 3 x^9 - y^13)(17 - 5 y - z^14)]]`) it ran for >120 s.  Routing the same gcd through subresultant PRS keeps coefficient sizes polynomially bounded and runs in sub-millisecond time on the same input, bringing the full Factor call to under 1 s.

```mathematica
In[1]:= FactorSquareFree[x^5 - x^3 - x^2 + 1]
Out[1]= (-1 + x)^2 (1 + 2 x + 2 x^2 + x^3)

In[2]:= FactorSquareFree[x^4 - 9x^3 + 29x^2 - 39x + 18]
Out[2]= (-3 + x)^2 (2 - 3 x + x^2)

In[3]:= FactorSquareFree[x^5 - x^3 y^2 - x^2 y^3 + y^5]
Out[3]= (x - y)^2 (x^3 + 2 x^2 y + 2 x y^2 + y^3)

In[4]:= FactorSquareFree[{(x^2 - 1)(x - 1), (x^4 - 1)(x^2 - 1)}]
Out[4]= {(1 + x) (-1 + x)^2, (1 + x^2) (-1 + x^2)^2}
```

## FactorSquareFreeList
Gives a list of the square-free factors of a polynomial together with their multiplicities.
- `FactorSquareFreeList[poly]`
- `FactorSquareFreeList[poly, Extension -> alpha]`
- `FactorSquareFreeList[poly, Modulus -> p]`

**Features**:
- `Protected` (not `Listable`, v0.6066: a list or an equation is one factor of itself, `{{1, 1}, {list, 1}}`).
- A thin wrapper over `FactorSquareFree`: it decomposes via `FactorSquareFree[poly, opts...]` and splits the product into `{factor, exponent}` pairs.  Options are validated under its own name (`FactorSquareFreeList::optx`, `::modp`).  The relationship to `FactorSquareFree` is exactly `FactorList`'s to `Factor`.
- The first element is the overall numerical factor `{c, 1}` — `{1, 1}` when there is none; a rational factor `p/q` is listed as `{p, 1}, {q, -1}`, as in `FactorList`.
- Each square-free factor keeps its multiplicity; unlike `FactorList`, the factors are *not* further reduced to irreducibles (`x^2 - 1` stays whole rather than splitting into `(-1 + x)(1 + x)`).
- `Times @@ Power @@@ FactorSquareFreeList[poly]` reconstructs `poly` (up to `FactorSquareFree`'s normal form).

```mathematica
In[1]:= FactorSquareFreeList[x^5 - x^3 - x^2 + 1]
Out[1]= {{1, 1}, {-1 + x, 2}, {1 + 2 x + 2 x^2 + x^3, 1}}

In[2]:= FactorSquareFreeList[2 x^3 + 2 x^2 - 2 x - 2]
Out[2]= {{2, 1}, {-1 + x, 1}, {1 + x, 2}}

In[3]:= FactorSquareFreeList[x^2 - 2 Sqrt[2] x + 2, Extension -> Sqrt[2]]
Out[3]= {{1, 1}, {-Sqrt[2] + x, 2}}
```

## Factor
Factors a polynomial over the rationals, or over an algebraic extension or a finite field.
- `Factor[poly]` -- over Q.
- `Factor[poly, Extension -> alpha]` / `Extension -> {a1, ...}` / `Extension -> Automatic` -- over Q(alpha).
- `Factor[poly, GaussianIntegers -> True]` -- over the Gaussian rationals Q(I).
- `Factor[poly, Modulus -> p]` -- over GF(p), p prime.
- `Factor[poly, Trig -> True]` -- accepted; no trigonometric rewriting is done.

**Features**:
- `Listable`, `Protected`.  `Options[Factor]` is `{Extension -> None, GaussianIntegers -> False, Modulus -> 0, Trig -> False}`; `SetOptions` applies.
- **Option handling** (v0.6062). Options are parsed once: every trailing argument must be a rule (`Factor::nonopt`) naming one of the four options (`Factor::optx`, e.g. `Method`, which `Factor` does not have); the first occurrence of an option wins.  A call whose options are all at their defaults (`Extension -> None`, `Modulus -> 0`, ...) is the plain form and takes every fast path; before, *any* option pushed a multivariate input onto the classical path (`Factor[x^99 - y^99, Extension -> None]`: over 100 s, now 4 ms).
- **`Modulus -> p`** (v0.6062) factors over GF(p) through `fmpz_mod_mpoly_factor`, for a prime of any size and any number of variables: monic factors with coefficients in `[0, p)`, the leading coefficient in front.  `Factor[3 x^2 + 3, Modulus -> 5]` is `3 (2 + x) (3 + x)`, `Factor[x^4 - y^4, Modulus -> 3]` is `(x + y) (x + 2 y) (x^2 + y^2)`, `Factor[x^4 + 1, Modulus -> 18446744073709551629]` works past a machine word.  A rational coefficient `a/b` is `a b^-1 mod p`; a rational function factors its numerator and denominator; lists and relations thread; a product or power factors base by base, so `(1 + x)^32768 (1 + x + x^2)` is never expanded over Q.  `Modulus -> 0` or a symbolic value means Q and `-p` means `p`.  Messages: `Factor::modp` (not a prime, e.g. 9 or 1), `::bmod` (not an integer), `::mext` (with an extension), `::facim` (with `GaussianIntegers -> True` or a complex coefficient), `::nmod` (a coefficient denominator divisible by p), `Power::infy` (a denominator that vanishes mod p).  A non-polynomial input is returned unchanged.
- **Fast entry** (v0.6063). A `Times` or an integer power factors base by base and multiplies back — factoring is multiplicative and equal factors merge — so `Factor[(x + 1)^500 (x - 1)^500]` takes 0.06 ms (1.44 s before, through the expansion) and `Factor[(x^2 - 1)/(x + 1)]` cancels to `-1 + x`.  A polynomial over Q in plain symbols goes straight to FLINT ahead of the pre-stage cascade.
- **Cyclotomic splitting.** `c x^k (x^m ± 1)` factors into cyclotomic polynomials (`fmpz_poly_cyclotomic`) instead of a general factorisation: `x^2310 - 1` (33 factors) went from 2.9 s to milliseconds.
- **Exponent bound.** A variable exponent of `2^20` or more is refused with `Factor::lrgexp`; it was read through an `int`, so `Factor[x^(2^32 + 2) + x^2 - 1]` returned `(-1 + x) (1 + x)`.  A large power of a sum is fine: `Factor[(x^2 - 1)^(2^31)]` is `(-1 + x)^2147483648 (1 + x)^2147483648`.
- **Kernels as generators** (v0.6065). A numeric base under a symbolic exponent (`E^x`, `2^x`, `Pi^x`, `(-1)^n`) and a positive fractional power (`x^(1/3)`, `(1 + x^2)^(1/4)`, `Sqrt[2]`) are generators of their own: `Factor[E^x y - E^x]` is `E^x (-1 + y)`, `Factor[(Sqrt[2] x^2 - Sqrt[2])/(x - 1)]` is `Sqrt[2] (1 + x)`.  Treating a radical as independent of the other variables can miss a factorisation, never invent one.
- **Rational functions** (v0.6061): each side of `Together`'s fraction goes to FLINT whatever its variable count.  A univariate side used to take the int64 Berlekamp-Zassenhaus path, which hung on `Factor[1/(x^43 - 1)]`, crashed on `Factor[1/(x^6001 + x + 1)]` and left wide-coefficient sides unfactored.
- **Inexact input** is made exact, factored, and numericalised.  A coefficient below `1e-4` was rationalised to `0` (`Factor[x - 0.000001]` gave `x`), and a real of `2^63` or more sent the round trip into infinite recursion (a crash); both fixed in v0.6061.
- **Gaussian-rational coefficients are fast** (v0.6051, issue #114). A
  multivariate polynomial with `I` in its coefficients is still factored over
  `Q(I)`, as in the reference (`Factor[Expand[(x + I)^2 (1 + y)]]` is
  `(I + x)^2 (1 + y)`), but FLINT's `fmpq_mpoly_factor` declines it, and the
  classical path was slow: a 28-term, 30-kernel irreducible input took 63 s
  and now takes 26 ms. Two changes did this. A polynomial of degree 1 in some
  variable with coprime coefficients is concluded irreducible (Gauss's lemma),
  provided every variable is a genuine indeterminate (not `Sqrt[3]`, a `Root`
  or an `AlgebraicNumber`). And the FLINT zero test now accepts Gaussian
  coefficients by mapping `I` to a fresh generator, instead of recursing
  through the classical zero test.
- **FLINT acceleration** (when built with FLINT): a univariate polynomial over Q is factored via `fmpz_poly_factor`, read in one pass (no coefficient extraction per degree), and a multivariate polynomial over Q via `fmpq_mpoly_factor`. This turns the exponential-in-variable-count classical multivariate factoriser — which hangs on inputs like `Factor[x^99 - y^99]` (> 20 s) — into a few-millisecond call. **Sign convention** (v0.6063): each irreducible factor has a positive *lexicographic* leading coefficient, the variables taken in canonical order (`a < A < aa < b < x < X < x1`, as `Sort` orders them) with the first most significant; the discarded signs fold into the content.  So `Factor[y^2 - x^2]` is `-(x - y)(x + y)` and `Factor[X^2 - a^2]` is `-((a - X)(a + X))` (it was `(-a + X)(a + X)` under a degree-lexicographic order with `strcmp`-sorted names). Inputs that are not a polynomial over Q (a denominator, fractional/symbolic exponent, or a non-polynomial head such as `Sqrt`/`Sin`) fall through to the classical path unchanged; `Extension`/`GaussianIntegers`/`Modulus` option forms travel their own branch and are untouched. See also the `` FLINT`Factor `` builtin.
- **Gaussian-integer content** (v0.6018), as in the reference. A polynomial
  with Gaussian-integer coefficients, at least one not real, takes its
  content over Z[i]:
  - `Factor[(1 + I) a + 2 b]` is `(1 + I) (a + (1 - I) b)`. The content is the
    first-quadrant associate of the gcd.
  - When what is left is purely imaginary, `I` joins the content:
    `Factor[4 I a - 4 I b]` is `(4 I) (a - b)`, and `Factor[I x^2 - I]` is
    `I (-1 + x) (1 + x)`.
  - In at most two variables the primitive part's sign moves into the
    content too (`Factor[-I a + I b]` is `-I (a - b)`). Over more variables
    the input's orientation is kept, as the reference keeps it.
  - A product with such a coefficient keeps it and factors the rest:
    `Factor[-(16 I) (56 - 30 D + 3 D^2) a b]` stays as it is. It used to come
    back as `16 a b (-56 I + (30 I) D - (3 I) D^2)`.
  - A univariate polynomial with Gaussian coefficients factors fully over
    Q(I) (v0.6069): `Factor[x^2 + 2 I x - 1]` is `(I + x)^2`, and
    `Factor[I x^2 + I]` is `I (-I + x) (I + x)` (an associate of the
    reference's `(1 + I x) (I + x)`: the unit sits in front).  A multivariate
    one still factors over Q unless `GaussianIntegers -> True` is given.
- **Number fields** (v0.6069). `Extension -> alpha | {a1, ...} | Automatic`
  and `GaussianIntegers -> True` work in any number of variables through
  `src/poly/flint_trager.c`: Trager's method, with a random shift until the
  norm is square-free, the norm from `fmpq_mpoly_resultant` and factored by
  `fmpq_mpoly_factor`, and each factor recovered with the certified modular
  field GCD.  The whole answer is checked by exact division before it is
  returned, and any decline falls back to the univariate `qafactor.c` path.
  Factors are spelled in the caller's own radicals, primitive over Z (over
  Z[i] for Q(I)), with the generator term leading and positive:
  `Factor[x^2 - 2 y^2, Extension -> Sqrt[2]]` is `-(-x + Sqrt[2] y) (x + Sqrt[2] y)`,
  `Factor[x^2 + y^2, GaussianIntegers -> True]` is `(x - I y) (x + I y)`.
  A polynomial over Q is factored over Q first and each irreducible factor
  over the field separately (v0.6071: `x^120 - 1` over Q(I) in 14 ms).
  `Extension -> Automatic` takes every generator in the input, and
  coefficients may be big rationals.  Timings against the old path
  (degree-40 products over fields of degree 2, 4, 6): 0.78 s -> 0.027 s,
  17.4 s -> 0.34 s, 87.8 s -> 0.39 s.
- When given a rational expression, first resolves dependencies over `Together` before factoring.
- **Powers with a symbolic exponent factor over their generators.** The input
  is rewritten with each `b^t` as a private symbol, factored, and put back
  together, as in the reference: `Factor[x^(2 + e) - x^e]` is
  `(-1 + x) x^e (1 + x)`, and `Factor[x^(1 + e) + x^(2 + e)]` is
  `x^(1 + e) (1 + x)`. A single power such as `(x + y)^(-3 + 3 e)` stays.
- Threads over lists, equations, inequalities and logical combinations.
- Bivariate inputs whose leading coefficient (in some variable) is the constant `-1` are handled via Wang's leading-coefficient correction, Stage 1: the input is pre-negated to make it monic, the existing monic Hensel pipeline runs on `-P`, and the overall sign is absorbed into the highest-degree factor via `Expand`.  This unlocks inputs of shape `Factor[(1 - x^k)(x - y^m)]` (and similar non-monic cases with constant `±1` LC) that previously fell back to the legacy linear-trial-division loop.
- Bivariate inputs whose leading coefficient (in some variable) is a non-unit integer constant `a` (with `|a| > 1`) are handled via Wang's leading-coefficient correction, Stage 2: the monic substitution `Q(x, y) = a^(d-1) · P(x/a, y)` makes the lift's input monic in x with integer coefficients.  After lifting `Q = G_1 · ... · G_r` via the existing pipeline, the true factors of `P` are recovered as `F_i = G_i(a·x, y) / cont_Z(G_i(a·x, y))`, where the integer content collects exactly the share of `a^(d-1)` redistributed into G_i by the substitution.  Stages 1 and 2 compose, so inputs with negative non-unit LC (e.g. `lc_x(P) = -6`) also enter the structured pipeline.  This unlocks inputs like `Factor[Expand[(2x+3y)(3x+5y)]]`, `Factor[2 a^2 - 5 a b + 3 b^2]`, and three-factor non-monic forms whose LCs are constant in y (or constant in x).
- Bivariate inputs whose leading coefficient (in some variable) is a non-constant polynomial in the other variable are handled via Wang's leading-coefficient correction, Stage 3 (predicted-LC two-factor Hensel).  When `lc_x(P)(y) = A(y)`, Mathilda factors `A` over `Z[y]`, finds `α` with `A(α) = +1` so the squarefree univariate image `P(x, α)` factors into monic Z[x] pieces `u`, `v`, then enumerates distributions of `A`'s irreducible factors between two predicted leading coefficients `q_u, q_v` (with `q_u · q_v = A` and `q_u(α) = q_v(α) = +1`).  The Hensel iteration is modified so each `Δu` correction has its leading-x coefficient PINNED to the y^k coefficient of `q_u`, keeping `lc_x(U)(y) = q_u(y)` invariant across the lift.  This unlocks inputs like `Factor[Expand[(xy+1)(xy+2)]]`, `Factor[Expand[((y²+1)x+1)(x+3)]]`, `Factor[Expand[((y+1)x+1)((y+1)x+2)]]`.  MVP scope: r = 2 (two univariate factors), both monic, `|cont(A)| = 1`, and inputs with non-trivial monomial content fall through so `heuristic_factor`'s Phase 0 path produces the canonical fully-factored form.
- **Binomial-descent soundness.** The heuristic two-term `a^k ± b^k` descent fires only when *both* terms of the `Plus` parse as genuine polynomial monomials (an integer times variables and integer powers). A term with a non-polynomial head — `Cos[r^2]`, `r^2 Sin[r^2]`, a generic `f[r^2]` — is no longer misread as the bare constant `1`. Previously the monomial parse silently dropped such factors, so `Factor[Cos[r^2] - r^2]` returned the unsound `-(-1 + r)(1 + r)`, and, because `Simplify` runs `Factor` as a candidate and keeps the shorter form, `Simplify[Laplacian[Sin[r^2], {r, t}, "Polar"]]` collapsed to `4 - 4 r^2`. Such inputs are now returned unfactored.
- **A product never comes back as a sum.** `Factor` normalises through `Together`, and `Together` expands a multivariate product: `Together[E^x (x^2 - 1)]` is `-E^x + x^2 E^x`. That expansion is not a polynomial in the single variable `collect_variables` finds (`x`) — its coefficients carry `E^x` — so there was nothing to factor and the expansion came straight back: `Factor[E^x (x^2 - 1)]` returned `-E^x + x^2 E^x`, *less* factored than its own input, and `Integrate[x^5 E^x, x] // Factor` multiplied the integrator's answer out. Factoring is multiplicative, so a **denominator-free** product is now factored factor by factor: `Factor[E^x (x^2 - 1)]` → `E^x (-1 + x) (1 + x)`, `Factor[Log[x] (x^2 - 1)]` → `(-1 + x) (1 + x) Log[x]`, and an irreducible polynomial factor is returned whole rather than expanded. The denominator case is deliberately excluded, because a sum is the right answer there: `Factor[(x^2 - 1)/(x - 1)]` is `1 + x`, and threading over that `Times` would hand back an uncancelled `(x - 1) (x + 1)/(x - 1)`. Under that guard the rewrite fires only where the polynomial path has already failed, so nothing that factors today changes route. (The root defect is `Together`'s expansion, left alone here: a great deal of the integrator, `DSolve` and `Simplify` is tuned around the shape it returns.)
- **Function kernels are opaque generators.** Every non-symbol generator (`BesselJ[1/3, x]`, `Log[x]`, `f[x]`, `E^x`) is mapped to a fresh symbol before factoring and restored afterwards. Before this, the multivariate paths specialised a *symbol* generator by structural substitution: the ≥3-generator z-split step ran `ReplaceAll[P, x -> 0]`, which also rewrote the `x` *inside* `BesselJ[-2/3, x]`. That substitution is not a specialisation of the polynomial at all. For `BesselJ[-2/3, 0]` it produced `0 ComplexInfinity`, so `Factor` (and `Simplify`, which calls it) printed `Infinity::indet`. For `f[x]` it silently specialised to `f[0]`. The answer was still correct because candidates are trial-divided, but `Check` took the failure branch. Now `Factor[x^2 BesselJ[-2/3, x]^2 + x BesselJ[1/3, x]^2 + 1]` emits no message, and the FLINT multivariate fast path, which declines on unknown heads, takes these inputs. The fresh symbols are named in the kernels' canonical order, so the sign/ordering convention of the result is unchanged (`Factor[x^2 a[x]^2 - b[x]^2]` → `(x a[x] - b[x]) (x a[x] + b[x])`).
- **Uncombined-fraction soundness.** `Factor` derives its numerator/denominator from `Together`, then factors each as a polynomial. When `Together` *declines* to combine a sum of fractions (notably a sum whose numerators share ≥2 dependent trig/hyperbolic kernels over `Q(i)`, which `Together` leaves uncombined to avoid a GCD blow-up), the "numerator" is still a `Plus` carrying `Power[linear, -1]` denominators — not a polynomial. Feeding it to the polynomial factorer treated `1/(x - I)` as an indeterminate and returned a **non-equivalent** result: `Factor[(-1/2*I) (Cosh[x] + Sinh[x])/(-I + x) + (1/2*I) (Cosh[x] + Sinh[x])/(I + x)]` collapsed to `0`, which (via the same `Simplify`-trusts-`Factor` path) made `D[Integrate[Exp[x]/(x^2 + 1), x], x] // Simplify` return `0` instead of `E^x/(1 + x^2)`. `Factor` now verifies the numerator/denominator are genuinely denominator-free before factoring and, when they are not, returns the input unchanged rather than a false factorisation. The well-combined cases (polynomials, single-trig, symbol sums, real-denominator sums, `f[x] + g[x]`, Gaussian-integer polynomials) are unaffected.

```mathematica
In[1]:= Factor[1 + 2x + x^2]
Out[1]= (1 + x)^2

In[2]:= Factor[x^10 - 1]
Out[2]= (-1 + x) (1 + x) (1 + x + x^2 + x^3 + x^4) (1 - x + x^2 - x^3 + x^4)

In[3]:= Factor[x^10 - y^10]
Out[3]= (x + y) (x - y) (x^4 + x^3 y + x^2 y^2 + x y^3 + y^4) (x^4 - x^3 y + x^2 y^2 - x y^3 + y^4)

In[4]:= Factor[2x^3 y - 2a^2 x y - 3a^2 x^2 + 3a^4]
Out[4]= (a + x) (a - x) (3 a^2 - 2 x y)

In[5]:= Factor[(x^3 + 2x^2)/(x^2 - 4y^2) - (x + 2)/(x^2 - 4y^2)]
Out[5]= ((-1 + x) (1 + x) (2 + x))/((x - 2 y) (x + 2 y))

In[6]:= Factor[x^4 + 1, Modulus -> 3]
Out[6]= (2 + x + x^2) (2 + 2 x + x^2)

In[6a]:= Factor[x^4 - y^4, Modulus -> 3]
Out[6a]= (x + y) (x + 2 y) (x^2 + y^2)

In[7]:= Factor[x^4 + 1, Extension -> Sqrt[2]]
Out[7]= (1 + Sqrt[2] x + x^2) (1 - Sqrt[2] x + x^2)

In[8]:= Factor[Integrate[f[x], x] x^2 - Integrate[f[x], x]]
Out[8]= (-1 + x) (1 + x) Integrate[f[x], x]
```

**v0.6035.** The transcendental constants Pi, E, EulerGamma and Catalan, and their fractional powers (`Sqrt[Pi]`), are indeterminates to the multivariate factoriser, so a rational expression carrying them reaches FLINT. A 13-generator Struve residual went from 61 s to 2.5 ms. `Factor[x^2 - Pi^2]` is `(-Pi + x) (Pi + x)`.

## FactorList
Gives a list of the irreducible factors of a polynomial together with their exponents.
- `FactorList[poly]`
- `FactorList[poly, GaussianIntegers -> True]`
- `FactorList[poly, Extension -> {a1, a2, ...}]`
- `FactorList[poly, Modulus -> p]`

**Features**:
- `Protected`; not `Listable` (v0.6066): `FactorList[{x^2 - 1, x}]` is `{{1, 1}, {{-1 + x^2, x}, 1}}` and an equation is likewise one factor of itself.  `Options[FactorList]` is `Factor`'s.
- A thin wrapper over `Factor`: it validates the options under its own name (`FactorList::optx`, `::modp`, ...; it used to wrap an unevaluated `Factor[...]` as a factor), factors via `Factor[poly, opts...]` and splits the product into `{factor, exponent}` pairs.
- **A power with a symbolic exponent lists its generators** as the reference
  does: `FactorList[(x + y)^(-3 + 3 e)]` is
  `{{1, 1}, {(x + y)^e, 3}, {x + y, -3}}`, and `FactorList[x^(1 + e) y^(2 - e)]`
  is `{{1, 1}, {x^e, 1}, {y^e, -1}, {x, 1}, {y, 2}}`. The generators come right
  after the constant, in order of appearance, and the ordinary factors keep
  `FactorList`'s own order.
- The first element is the overall numerical factor `{c, 1}` — `{1, 1}` when there is none.  A rational factor `p/q` is listed as `{p, 1}, {q, -1}`, the `{p, 1}` dropped when `p` is 1 (v0.6061): `FactorList[-15/2 (x^2 - 1)]` is `{{-15, 1}, {2, -1}, {-1 + x, 1}, {1 + x, 1}}`, `FactorList[(x^2 - 1)/(6 x + 6)]` is `{{6, -1}, {-1 + x, 1}}`.
- Every factor has a positive leading coefficient (lexicographic, as `CoefficientRules` orders terms), its sign going into the numerical factor: `FactorList[-x - 1]` is `{{-1, 1}, {1 + x, 1}}`, `FactorList[y - x]` is `{{-1, 1}, {x - y, 1}}`. Needed since `-(1 + x)` itself evaluates to `-1 - x` (v0.591).
- Denominator factors of a rational function appear with negative exponents.
- `Times @@ Power @@@ FactorList[poly]` reconstructs `poly` (up to `Factor`'s normal form).
- Because it delegates to `Factor`, the factorisation depth and normal form (factor ordering, `Extension`/`GaussianIntegers` handling) are exactly `Factor`'s.

```mathematica
In[1]:= FactorList[x^2 - 1]
Out[1]= {{1, 1}, {-1 + x, 1}, {1 + x, 1}}

In[2]:= FactorList[2 x^3 + 2 x^2 - 2 x - 2]
Out[2]= {{2, 1}, {-1 + x, 1}, {1 + x, 2}}

In[3]:= FactorList[(x^3 + 2 x^2)/(x^2 - 4 y^2) - (x + 2)/(x^2 - 4 y^2)]
Out[3]= {{1, 1}, {-1 + x, 1}, {1 + x, 1}, {2 + x, 1}, {x - 2 y, -1}, {x + 2 y, -1}}

In[4]:= FactorList[x^2 + 1, Extension -> I]
Out[4]= {{1, 1}, {-I + x, 1}, {I + x, 1}}
```

## FactorTerms
Pulls out an overall numerical factor in a polynomial, or factors that do
not depend on a given set of variables.
- `FactorTerms[poly]` -- pulls out any overall numerical factor in `poly`.
- `FactorTerms[poly, x]` -- pulls out any overall factor in `poly` that
  does not depend on `x` (i.e. extracts the polynomial content of `poly`
  with respect to `x`).
- `FactorTerms[poly, {x_1, x_2, ...}]` -- pulls out any overall factor
  in `poly` that does not depend on any of the `x_i`. The result is then
  recursively refined by extracting content with respect to the smaller
  subsets `{x_1, ..., x_{n-1}}`, ..., `{x_1}`.

**Features**:
- `Protected`.  `FactorTerms[poly, Modulus -> p]` reduces the coefficients mod p first (`FactorTerms[3 x + 6, Modulus -> 3]` is `0`).
- The numerical factor is signed so the rest has a positive leading term (lexicographic, variables in canonical order): `FactorTerms[-6 x^2 - 12 x y]` is `-6 (x^2 + 2 x y)` (v0.6061).
- Auto-threads over `List`, `Equal`, `Unequal`, `Less`, `LessEqual`,
  `Greater`, `GreaterEqual`, `And`, `Or`, `Not`, `Xor`.
- Together-normalises rational inputs before extracting content from the
  numerator, then divides the denominator back through, so rational
  functions round-trip.
- Numerical content over `Z` is computed via the integer GCD of monomial
  coefficients. Gaussian-integer content (e.g. `5 I` from `5 I x^2 + ...`)
  is not extracted as a Gaussian unit; the integer GCD `5` is returned
  instead, with the leading factor of `I` left inside the residue. The
  resulting factorization is mathematically equivalent.
- The content of a nonzero input is always a unit — never `0`. If an upstream
  `Together` leaves a Laurent term (e.g. `b^-1`) inside the numerator, the
  nonnegative-degree coefficient extraction sees no monomials; the content is
  then taken as `1` (nothing is factored, the residue is intact) rather than a
  spurious `0`. (A `0` content formerly multiplied the whole result to `0` and
  silently poisoned `Simplify`, which trusts a transform's literal `0`.)

```mathematica
In[1]:= FactorTerms[3 + 6x + 3x^2]
Out[1]= 3 (1 + 2 x + x^2)

In[2]:= FactorTerms[3 + 3a + 6 a x + 6 x + 12 a x^2 + 12 x^2, x]
Out[2]= 3 (1 + a) (1 + 2 x + 4 x^2)

In[3]:= FactorTerms[12 a^4 + 9 x^2 + 66 b^2]
Out[3]= 3 (4 a^4 + 22 b^2 + 3 x^2)

In[4]:= FactorTerms[7 x + (14 y + 21)/z]
Out[4]= (7 (3 + 2 y + x z))/z

In[5]:= FactorTerms[{5 x^2 - 15, 7 x^4 - 77, 8 x^8 - 24}]
Out[5]= {5 (-3 + x^2), 7 (-11 + x^4), 8 (-3 + x^8)}

In[6]:= FactorTerms[1 < 77 x^3 - 21 x + 35 < 2]
Out[6]= 1 < 7 (5 - 3 x + 11 x^3) < 2

In[7]:= f = 2 x^2 y z + 2 x^2 y + 4 x^2 z + 4 x^2 + 4 y^2 z^2 + 4 z y^2
        + 8 z^2 y + 2 z y - 6 y - 12 z - 12;
        FactorTerms[f, x]
Out[7]= 2 (2 + y + 2 z + y z) (-3 + x^2 + 2 y z)

In[8]:= FactorTerms[f, {x, y}]
Out[8]= 2 (2 + y) (1 + z) (-3 + x^2 + 2 y z)
```

## FactorTermsList
Lists the factors that `FactorTerms` would multiply together.
- `FactorTermsList[poly]` -- gives `{numerical_factor, residue}`.
- `FactorTermsList[poly, x]` -- gives
  `{numerical_factor, x_independent_factor, residue}`.
- `FactorTermsList[poly, {x_1, ..., x_n}]` -- gives a list whose first
  element is the overall numerical factor; the second is a factor that
  does not depend on any of the `x_i`; each subsequent element is a
  factor depending on progressively more of the `x_i`; and the final
  element is the residue.

**Features**:
- `Protected`.
- The product of the returned list always reproduces the input (after
  Together-normalisation), so `Apply[Times, FactorTermsList[poly]]` is a
  faithful round-trip.
- Variables in the second argument that do not actually appear in `poly`
  are filtered out, so the output never contains spurious trailing `1`s.
- The numerical factor (v0.6066) is rational and signed: `FactorTermsList[-6 x^2 - 12 x y]` is `{-6, x^2 + 2 x y}`, `FactorTermsList[3/4 x + 3/2]` is `{3/4, 2 + x}` (the integer content of a denominator belongs to it: `FactorTermsList[(4 x + 4)/(6 y)]` is `{2/3, 1/y + x/y}`), and `FactorTermsList[0]` is `{1, 0}`.  A list argument is `General::poly`.  `Modulus -> p` reduces the coefficients first.

```mathematica
In[1]:= FactorTermsList[3 + 6 x + 3 x^2]
Out[1]= {3, 1 + 2 x + x^2}

In[2]:= FactorTermsList[14 x + 21 y + 35 x y + 63]
Out[2]= {7, 9 + 2 x + 5 x y + 3 y}

In[3]:= FactorTermsList[3 + 3 a + 6 a x + 6 x + 12 a x^2 + 12 x^2, x]
Out[3]= {3, 1 + a, 1 + 2 x + 4 x^2}

In[4]:= FactorTermsList[-6 y - 6 a y + 2 x^2 y + 2 a x^2 y + 4 a y^2
                        + 4 a^2 y^2, {x, y}]
Out[4]= {2, 1 + a, y, -3 + x^2 + 2 a y}

In[5]:= Times @@ FactorTermsList[14 x + 21 y + 35 x y + 63]
Out[5]= 7 (9 + 2 x + 5 x y + 3 y)
```

## PolynomialGCD
Gives the greatest common divisor of the polynomials.
- `PolynomialGCD[poly1, poly2, ...]`
- `PolynomialGCD[poly1, poly2, ..., Extension -> alpha]`

**Features**:
- `Protected`, `Listable`.
- Handles univariate and multivariate polynomials.
- Treats algebraic numbers (like `I`) as independent variables or constants seamlessly during complex arithmetic evaluations.
- Pre-extracts common factors before falling back to a full primitive Euclidean algorithm computation.
- **Option `Extension -> alpha`** (Phase 0 of the Integrate plan): computes the GCD over `Q(alpha)` for `alpha` ∈ {`Sqrt[c]`, `c^(1/n)`, `I`} via lifting both inputs into the QAUPoly substrate (`src/poly/qaupoly.h`) and folding `qaupoly_gcd`. Extension support requires univariate inputs (after stripping the alpha-render symbol). Defaults `Extension -> None` and `Extension -> Automatic` work over the rationals and treat algebraic numbers as opaque variables. `Extension -> {alpha_1, ..., alpha_n}` (tower form) currently falls back to the no-extension path; tower-aware GCD is a Phase 0.5 follow-up.

```mathematica
In[1]:= PolynomialGCD[(1+x)^2(2+x)(4+x), (1+x)(2+x)(3+x)]
Out[1]= (1+x) (2+x)

In[2]:= PolynomialGCD[x^2+4x+4, x^2+2x+1]
Out[2]= 1

In[3]:= PolynomialGCD[x^2-1, x^3-1, x^4-1, x^5-1, x^6-1, x^7-1]
Out[3]= -1 + x

In[4]:= PolynomialGCD[x^2 - 2, x - Sqrt[2], Extension -> Sqrt[2]]
Out[4]= -Sqrt[2] + x

In[5]:= PolynomialGCD[x^3 - 2, x - 2^(1/3), Extension -> 2^(1/3)]
Out[5]= -2^(1/3) + x
```

## PolynomialExtendedGCD
Gives the extended GCD of polynomials.
- `PolynomialExtendedGCD[poly1, poly2, x]`
- `PolynomialExtendedGCD[poly1, poly2, x, Modulus -> p]`

**Features**:
- `Protected`.
- Returns `{d, {a, b}}` such that $a \cdot poly1 + b \cdot poly2 = d$.
- $d$ is the GCD, normalized to be monic.
- Efficiently handles termination when a constant remainder is reached.
- Optimized for cases where the divisor is a constant.

```mathematica
In[1]:= PolynomialExtendedGCD[2x^5-2x, (x^2-1)^2, x]
Out[1]= {-1 + x^2, {x/4, (-4 - 2 x^2)/4}}

In[2]:= PolynomialExtendedGCD[a (x+b)^2, (x+a)(x+b), x]
Out[2]= {b + x, {-(1/(a (a - b))), 1/(a - b)}}
```

## PolynomialLCM
Gives the least common multiple of the polynomials.
- `PolynomialLCM[poly1, poly2, ...]`
- `PolynomialLCM[poly1, poly2, ..., Extension -> alpha]`

**Features**:
- `Protected`, `Listable`.
- Handles univariate and multivariate polynomials.
- Treats algebraic numbers (like `I`) as independent variables or constants seamlessly during complex arithmetic evaluations.
- Preserves explicit factored forms where possible.
- **Option `Extension -> alpha`** computes the LCM over `Q(alpha)` via `lcm(a, b) = a*b / PolynomialGCD[a, b, Extension -> alpha]`, returning the monic, expanded form. Same scope and fallback as `PolynomialGCD`'s extension option.

```mathematica
In[1]:= PolynomialLCM[(1+x)^2(2+x)(4+x), (1+x)(2+x)(3+x)]
Out[1]= (1+x)^2 (2+x) (3+x) (4+x)

In[2]:= PolynomialLCM[x^4-4, x^4+4 x^2+4]
Out[2]= (-2+x^2) (4+4 x^2+x^4)

In[3]:= PolynomialLCM[x - Sqrt[2], x + Sqrt[2], Extension -> Sqrt[2]]
Out[3]= -2 + x^2
```

## PolynomialQuotient
Gives the quotient of $p$ and $q$, treated as polynomials in $x$, with any remainder dropped.
- `PolynomialQuotient[p, q, x]`
- `PolynomialQuotient[p, q, x, Extension -> alpha]` -- divide over $\mathbb{Q}(\alpha)$.

**Features**:
- `Protected`.
- Default path uses polynomial long division over the field of rational functions in the coefficients.
- Option `Extension -> alpha` (default `None`) lifts $p, q$ into $\mathbb{Q}(\alpha)[x]$ and runs the Q($\alpha$)-aware long division (`qaupoly_divrem`). Recognised forms for $\alpha$: `Sqrt[c]`, `c^(1/n)`, and `I`. `Extension -> None` and `Extension -> Automatic` are accepted and currently behave as the default. The extension path requires univariate input (a single live polynomial variable other than the alpha generator); multivariate inputs fall through to the standard path.
- Threading Extension here keeps the polynomial arithmetic in the Q($\alpha$)[x] substrate and avoids the multivariate Q[$\alpha$, x] subresultant-PRS path that is exponentially slow on Sqrt[$\alpha$]-laden coefficients.

```mathematica
In[1]:= PolynomialQuotient[x^4+2x+1, x^2+1, x]
Out[1]= -1 + x^2

In[2]:= PolynomialQuotient[x^2+2x+1, x^3, x]
Out[2]= 0

In[3]:= PolynomialQuotient[x^2+x+1, 2x+1, x]
Out[3]= 1/4 + x/2

In[4]:= PolynomialQuotient[x^2 - 2, x - Sqrt[2], x, Extension -> Sqrt[2]]
Out[4]= Sqrt[2] + x

In[5]:= PolynomialQuotient[x^3 - 2, x - 2^(1/3), x, Extension -> 2^(1/3)]
Out[5]= 2^(2/3) + 2^(1/3) x + x^2

In[6]:= PolynomialQuotient[x^2 + 1, x - I, x, Extension -> I]
Out[6]= I + x
```

## PolynomialRemainder
Gives the remainder from dividing $p$ by $q$, treated as polynomials in $x$.
- `PolynomialRemainder[p, q, x]`
- `PolynomialRemainder[p, q, x, Extension -> alpha]` -- compute the remainder over $\mathbb{Q}(\alpha)$.

**Features**:
- `Protected`.
- The degree of the result in $x$ is guaranteed to be smaller than the degree of $q$.
- If the dividend is a multiple of the divisor, then the remainder is zero.
- Option `Extension -> alpha`: see `PolynomialQuotient` for the recognised alpha forms and the fall-through rules.

```mathematica
In[1]:= PolynomialRemainder[x^4+2x+1, x^2+1, x]
Out[1]= 2 + 2 x

In[2]:= PolynomialRemainder[x^3, a x+b, x]
Out[2]= -(b^3/a^3)

In[3]:= PolynomialRemainder[x^2 - 2, x - Sqrt[2], x, Extension -> Sqrt[2]]
Out[3]= 0

In[4]:= PolynomialRemainder[x^2 - 3, x - Sqrt[2], x, Extension -> Sqrt[2]]
Out[4]= -1
```

## PolynomialQ
Tests whether an expression is a polynomial in one or more variables.
- `PolynomialQ[expr, var]`: Yields `True` if `expr` is a polynomial in `var`.
- `PolynomialQ[expr, {var1, var2, ...}]`: Tests whether `expr` is a polynomial in the set of variables.

**Features**:
- `Protected`.
- Variables can be symbols or compound expressions.
- Constants (expressions free of the specified variables) are polynomials of degree 0.
- `Power[base, exp]` is a polynomial if `exp` is a non-negative integer and `base` is a polynomial.

```mathematica
In[1]:= PolynomialQ[x^3 - 2x/y + 3x z, x]
Out[1]= True

In[2]:= PolynomialQ[x^3 - 2x/y + 3x z, y]
Out[2]= False

In[3]:= PolynomialQ[x^2 + a x y^2 - b Sin[c], {x, y}]
Out[3]= True

In[4]:= PolynomialQ[f[a] + f[a]^2, f[a]]
Out[4]= True
```

## PolynomialMod
Gives the polynomial reduced modulo `m`.
- `PolynomialMod[poly, m]`
- `PolynomialMod[poly, {m1, m2, ...}]`

**Features**:
- `Protected`, `Listable`.
- Reduces a polynomial modulo an integer, another polynomial, or a list of integers/polynomials.
- Always gives a result with minimal degree and leading coefficients.
- Handles rational division mapping perfectly scaling over exact modulo structures dynamically.

```mathematica
In[1]:= PolynomialMod[3x^2+2x+1,2]
Out[1]= 1 + x^2

In[2]:= PolynomialMod[3x^2+2x+1,x^2+1]
Out[2]= -2 + 2 x

In[3]:= PolynomialMod[35x^3+21x^2 y^2-17x y^3+55z-123,19]
Out[3]= 10 + 2 x y^3 + 16 x^3 + 2 x^2 y^2 + 17 z

In[4]:= PolynomialMod[3x^3+21x^2 y^2-7x y^3+55,{2x^2-7,x y-3, 9}]
Out[4]= 1 + 7 x + x^3 + 4 y^2
```

## Resultant
Computes the resultant of two polynomials.
- `Resultant[poly1, poly2, var]`

**Features**:
- `Protected`, `Listable`.
- Computes the resultant of polynomials `poly1` and `poly2` with respect to the variable `var`.
- The resultant is independent of common roots and vanishes exactly when the polynomials have roots in common.
- Default algorithm is Bronstein's subresultant PRS (Symbolic Integration I, p.24): a linear chain of pseudo-remainders with scalar exact divisions in the coefficient ring, avoiding the (n+m)x(n+m) Sylvester matrix construction and its O(n^3) Bareiss reduction.  For Z/Q coefficients this is materially faster than the matrix path, and on inputs with symbolic coefficients it sidesteps the O(n!) Laplace expansion that the matrix path falls back to when Bareiss exact-division certification fails.
- Inputs containing algebraic-number coefficients (e.g. `Sqrt[N]`, cube roots — any `Power[X, Rational[a,b]]` with `b > 1`) are routed to the Sylvester+Det path instead, because the subresultant chain bloats geometrically when `Power[base, k/m]` forms can't be combined with their `Times[base^q, Sqrt[base]]` equivalents by `Plus` alone.
- A size-budget guard inside the subresultant path falls back to Sylvester+Det for any pathological input where chain elements exceed ~30x the input leaf-count.
- Automatically preserves multiplicativity (e.g., $Res(A \cdot B, Q) = Res(A, Q) Res(B, Q)$ and $Res(A^k, Q) = Res(A, Q)^k$).
- **FLINT acceleration** (when built with FLINT): plain rational inputs are computed via `fmpq_mpoly_resultant`, which avoids the subresultant PRS coefficient growth on higher-degree / multivariate inputs. FLINT's convention matches the classical output exactly; inputs with algebraic-number or otherwise non-rational coefficients fall through to the paths above. The same kernel is exposed directly as `` FLINT`Resultant `` (see the FLINT` context section).

```mathematica
In[1]:= Resultant[x^2 - 2x + 7, x^3 - x + 5, x]
Out[1]= 265

In[2]:= Resultant[x^3 - 5x^2 - 7x + 14, x^3 - 8x^2 + 9x + 58, x]
Out[2]= 0
```

## Discriminant
Computes the discriminant of the polynomial with respect to the variable.
- `Discriminant[poly, var]`

**Features**:
- `Protected`, `Listable`.
- Computes the discriminant of polynomial `poly` with respect to `var`.
- The discriminant is zero if and only if the polynomial has multiple roots.
- Derived symbolically utilizing the formula $D = \frac{(-1)^{n(n-1)/2}}{a_n} Resultant(P, P', var)$.

```mathematica
In[1]:= Discriminant[a x^2 + b x + c, x]
Out[1]= b^2 - 4 a c

In[2]:= Discriminant[5 x^4 - 3 x + 9, x]
Out[2]= 23273325

In[3]:= Discriminant[(x-1)(x-2)(x-3), x]
Out[3]= 4

In[4]:= Discriminant[(x-1)(x-2)(x-1), x]
Out[4]= 0
```

## HornerForm
Puts a polynomial or rational function into Horner form.
- `HornerForm[poly]`
- `HornerForm[poly, vars]`
- `HornerForm[poly1/poly2, vars1, vars2]`

## GroebnerBasis
Computes a Gröbner basis for the ideal generated by a list of polynomials.
- `GroebnerBasis[{p1, p2, ...}, {x1, x2, ...}]`
- `GroebnerBasis[{p1, p2, ...}, {x1, x2, ...}, {y1, y2, ...}]` — eliminates
  the `yi`.

**Features**:
- `Protected`.
- Computes the *reduced* Gröbner basis under the requested monomial order
  via Buchberger's algorithm.  Pair selection uses the normal strategy
  (smallest LCM total degree first); pair pruning uses Buchberger's
  coprime-LMs criterion together with Gebauer-Möller's M / F / B chain
  criteria.  Each basis polynomial is returned in primitive-over-Z form
  with positive leading coefficient and the basis is sorted ascending by
  leading monomial — matching the reference's output convention.
- Polynomial *equations* (`Equal[a, b]`) are accepted in place of
  polynomials; each is rewritten as `a - b` and `Expand`ed before
  conversion.
- The single-variable shorthand `GroebnerBasis[polys, x]` is accepted as
  a synonym for `GroebnerBasis[polys, {x}]`.
- **Parameters**: any free symbol that appears in the polynomials but is
  not in the explicit `vars`, the elimination block, or the
  `ParameterVariables` list — and is not a known mathematical constant
  (`Pi`, `E`, `EulerGamma`, `Catalan`, `Degree`, `Infinity`, `True`,
  `False`, ...) — is automatically promoted to a parameter and appended
  to the joint variable array in last-priority lex position.  This is
  the conventional "Gröbner basis over `Z[params][main_vars]`" form;
  consistency constraints on the parameter values themselves (purely-
  parameter polynomials) are kept in the output.
- **Cooperative abort**: the Buchberger loop calls `tc_check_deadline()`
  once per S-pair iteration so `TimeConstrained[GroebnerBasis[...], t]`
  aborts cleanly when the budget expires (returns `$Aborted`).
- **Performance note**: pure-lex Gröbner-basis computation is subject to
  exponential coefficient blowup on many systems with three or more
  variables; if `GroebnerBasis[polys, vars]` is slow, switching to
  `MonomialOrder -> DegreeReverseLexicographic` is often dramatically
  faster.  `Method -> "GroebnerWalk"` computes a cheap grevlex basis and
  then converts it to the requested order along the Gröbner fan; for hard
  lex (or weight-matrix) targets this is usually much faster than direct
  Buchberger while returning the identical reduced basis.

**Options**:

| Option               | Default          | Supported                                           |
| -------------------- | ---------------- | --------------------------------------------------- |
| `MonomialOrder`      | `Lexicographic`  | `Lexicographic`, `DegreeLexicographic`, `DegreeReverseLexicographic`, `EliminationOrder`, or an integer **weight matrix** `{{...}, ...}` (one column per variable) defining a custom term order. The `Negative*` orders are not well-founded (global) term orders, so they fall back to `Lexicographic` with a note. |
| `CoefficientDomain`  | `Rationals`      | `Rationals`, `Automatic`, `Integers` (strong Gröbner basis over Z), `RationalFunctions` (basis over the field Q(params)), `Polynomials[x, ...]` (named symbols become coefficient-ring parameters), `InexactNumbers[p]` (monic basis with `p`-digit real coefficients) |
| `Method`             | `Automatic`      | `"Buchberger"`, `"GroebnerWalk"`, `Automatic`        |
| `Sort`               | `False`          | `True` reverses the main-variable list before the basis is computed (matching the reference's empirical behaviour). |
| `ParameterVariables` | (unset)          | Symbol or list of symbols treated as coefficient-ring parameters; with no positional `vars` argument the main-variable list is auto-derived from the polynomials. |
| `Modulus`            | `0`              | `Modulus -> p` with **prime** `p` in `[2, 2^31)` computes the basis over the finite field `GF(p)` (`src/poly/gbmod.c`, `Lexicographic` / `DegreeReverseLexicographic`). A non-prime `p`, an out-of-range `p`, or a coefficient with no image in `GF(p)` falls back to the rational basis with a one-shot note. Not combined with elimination or `RationalFunctions` parameters in this version. |

A **weight-matrix** `MonomialOrder` must be a rectangular list of integer
rows with exactly one column per main variable, and must define a valid
term order: full column rank (so distinct monomials never tie) and a
positive first non-zero entry in every column (so each variable exceeds
1, making the order well-founded).  Matrices that fail either condition,
or that are used together with elimination/parameter variables, emit
`GroebnerBasis::nimpl` and fall back to `Lexicographic`.

Settings outside this table emit a `GroebnerBasis::nimpl` diagnostic to
`stderr` and fall back to the default.  `MonomialOrder ->
DegreeLexicographic` and the other named orders are built from a shared
named-order weight matrix (`gb_build_order_matrix` in `src/poly/groebner.c`,
also used by `PolynomialReduce` and the `CoefficientRules` / `MonomialList`
family), so a named order and its explicit weight-matrix spelling return
the identical basis.

`CoefficientDomain -> Integers` computes a strong Gröbner basis over the
ring Z (any parameter symbol is treated as an ordinary variable of the
polynomial ring Z[vars]).  Because the leading coefficients are not units,
the strong basis generally carries more generators than the rational
basis — completion forms, per pair, both the S-polynomial (coefficient
lcm) and a Bézout G-polynomial (coefficient gcd), and reduction cancels a
term only when the divisor's leading coefficient divides it.  Polynomial
content is preserved (it is ideal data over Z), so a constant generator
`c` denotes the ideal `c·Z[x]` rather than `{1}`.

`CoefficientDomain -> Polynomials[x, ...]` makes the named symbols
coefficient-ring variables (equivalent to listing them in
`ParameterVariables`).  `CoefficientDomain -> InexactNumbers[p]` computes
the exact basis, makes each polynomial monic, and renders coefficients as
`p`-digit machine/MPFR reals; inexact (`Real`) input coefficients are
rationalised before the computation.  This is exact-compute-then-round,
not a genuine floating-point Buchberger run.

`CoefficientDomain -> RationalFunctions` computes the basis over the field
Q(params), where the parameters (the non-main symbols) are units.  The
ring basis over Q[params][main vars] is autoreduced — generators whose
main leading monomial is a multiple of another's are dropped — yielding
the (typically smaller) field basis; coefficients are kept as primitive
integer polynomials in the parameters.  Honoured in the plain form only
(no elimination block) and equivalent to `Rationals` when there are no
parameters.

```mathematica
In[1]:= GroebnerBasis[{x^2 - 2 y^2, x y - 3}, {x, y}]
Out[1]= {-9 + 2 y^4, 3 x - 2 y^3}

In[2]:= GroebnerBasis[{x + y, x^2 - 1, y^2 - 2 x}, {x, y}]
Out[2]= {1}

In[3]:= GroebnerBasis[{(x-1)(x-2), (x-2)(x-3)}, x]
Out[3]= {-2 + x}            (* univariate -> PolynomialGCD *)

In[4]:= GroebnerBasis[{x + 2 y + 3 z, 4 x + 5 y + 6 z, 7 x + 8 y + 9 z},
                      {x, y, z}]
Out[4]= {y + 2 z, x - z}    (* linear -> Gaussian elimination *)

In[5]:= GroebnerBasis[{x^2 - 2 y^2 == 1, x y == 3}, {x, y}]
Out[5]= {-9 + y^2 + 2 y^4, 3 x - y - 2 y^3}

In[6]:= GroebnerBasis[{x^2 + y^2 + z^2 - 1, x y - z + 2, z^2 - 2 x + 3 y},
                      {x, y}, {z}]
Out[6]= {28 + 8 y + 41 y^2 - 26 y^3 + 55 y^4 - 8 y^5 + 7 y^6 - 6 y^7 + y^8,
         16 + 20 x - 10 y + 28 y^2 - 57 y^3 + 4 y^4 - 6 y^5 + 6 y^6 - y^7}

In[7]:= GroebnerBasis[{x^2 + y^2 + z^2 - 1, x - z + 2, z^2 - x y},
                      {x, y, z},
                      MonomialOrder -> DegreeReverseLexicographic]
Out[7]= {2 + x - z, -2 y + y z - z^2, 3 + y^2 - 4 z + 2 z^2,
         -6 + 4 y + 11 z - 6 z^2 + 3 z^3}

(* Implicit parameter: `a` is not in vars -> auto-promoted. *)
In[8]:= GroebnerBasis[{a x^2 + 5 x - 1, 2 x + 3 x y + y^2}, {x, y}]
Out[8]= {-4 - 12 y - 19 y^2 - 15 y^3 + a y^4,
         18 + 4 a x + 27 y + 45 y^2 + 2 a y^2 - 3 a y^3,
         2 x + 3 x y + y^2}

(* Subset-of-variables: `x` is auto-promoted to a parameter. *)
In[9]:= GroebnerBasis[
          {-5 x^2 + y z - x - 1, 2 x + 3 x y + y^2, x - 3 y + x z - 2 z^2},
          {y, z}]
Out[9]= {2 + 8 x + 113 x^2 + 212 x^3 + 938 x^4 + 364 x^5 + 1980 x^6
           - 1000 x^7 + 1625 x^8,
         -289274 - 1958323 x - 7458871 x^2 - 14140053 x^3 - 20245223 x^4
           - 19904560 x^5 + 4016500 x^6 - 24451375 x^7 + 202126 z,
         -46518 + 1236672 x + ... + 202126 y}

(* ParameterVariables -> ...: main vars auto-derived from polys. *)
In[10]:= GroebnerBasis[polys, ParameterVariables -> x]
Out[10]= (* same basis as In[9] *)

(* CoefficientDomain -> RationalFunctions: basis over the field Q(a). *)
In[11]:= GroebnerBasis[{a x^2 + 5 x - 1, 2 x + 3 x y + y^2}, {x, y},
                       CoefficientDomain -> RationalFunctions]
Out[11]= {-4 - 12 y - 19 y^2 - 15 y^3 + a y^4,
          18 + 4 a x + 27 y + 45 y^2 + 2 a y^2 - 3 a y^3}

(* CoefficientDomain -> InexactNumbers[p]: monic, approximate output. *)
In[12]:= GroebnerBasis[{a x^2 + 5 x - 1, 2 x + 3 x y + y^2}, {x, y},
                       CoefficientDomain -> InexactNumbers[20]]
Out[12]= {-4.0 - 12.0 y - 19.0 y^2 - 15.0 y^3 + 1.0 a y^4,
          4.5 + 1.0 a x + 6.75 y + 11.25 y^2 + 0.5 a y^2 - 0.75 a y^3,
          0.666... x + 1.0 x y + 0.333... y^2}

(* CoefficientDomain -> Integers: strong Gröbner basis over Z. *)
In[13]:= GroebnerBasis[{a x^2 + 5 x - 1, 2 x + 3 x y + y^2}, {x, y},
                       CoefficientDomain -> Integers]
Out[13]= {-4 - 12 y - 19 y^2 - 15 y^3 + a y^4,
          18 + 4 a x + 27 y + 45 y^2 + 2 a y^2 - 3 a y^3,
          2 x + 3 x y + y^2,
          6 + 2 a x + 9 y + a x y + 15 y^2 + a y^2 - a y^3,
          -1 + 5 x + a x^2}

(* Sort -> True reverses the variable list. *)
In[11]:= GroebnerBasis[polys, {x, y, z}, Sort -> True] ===
         GroebnerBasis[polys, {z, y, x}]
Out[11]= True

(* TimeConstrained correctly aborts the Buchberger run. *)
In[12]:= TimeConstrained[GroebnerBasis[hard_polys, {x, y, z}], 3]
Out[12]= $Aborted

(* Custom term order from a weight matrix (one column per variable). *)
In[13]:= GroebnerBasis[{-5x^2 + y z - x - 1, 2x + 3 x y + y^2,
                        x - 3y + x z - 2 z^2}, {x, y, z},
                       MonomialOrder -> {{1, 3, 1}, {-1, 2, 0}, {4, -3, 2}}]
Out[13]= (* Gröbner basis under the weight-matrix order *)

(* The Gröbner walk returns the identical reduced basis to Buchberger. *)
In[14]:= GroebnerBasis[polys, {x, y, z}, Method -> "GroebnerWalk"] ===
         GroebnerBasis[polys, {x, y, z}]
Out[14]= True
```

## PolynomialReduce
Multivariate polynomial division with cofactors.
- `PolynomialReduce[poly, {p1, p2, ...}, {x1, x2, ...}]`
- `PolynomialReduce[poly, {p1, p2, ...}]` — variables from `Variables`.
- `PolynomialReduce[poly, {p1, ...}, x]` — single-variable shorthand.

Returns `{{a1, a2, ...}, b}` such that `a1 p1 + a2 p2 + ... + b == poly`
exactly, with `b` *minimal*: no term of `b` is divisible by any leading term
of the `pi` under the chosen monomial order.  When the `pi` are a Gröbner
basis this uniquely determines the remainder `b`, so `b == 0` is an exact
ideal-membership test.

**Features**:
- `Protected`.
- The engine is `gb_divmod` (`src/poly/groebner.c`), the same normal-form
  division Buchberger uses; over Q with `Lexicographic` order it dispatches
  to FLINT's `fmpq_mpoly_divrem_ideal` when FLINT is present.
- Shares `GroebnerBasis`'s monomial orders (including `DegreeLexicographic`
  via the shared named-order weight matrix) and its parameter discovery: a
  free symbol outside the variable list is a coefficient-field parameter.
- Exact constants outside Q, such as `Pi`, `E` and `Sqrt[2]`, are field
  coefficients even when no free parameter is present (v0.537), so
  `PolynomialReduce[b c, {Pi - b c}, {b, c}]` → `{{-1}, Pi}`. Before this, the
  rational engine declined such input and the call stayed unevaluated.
- The quotients and the remainder depend on the monomial order and on the
  order of the `pi` (as in `GroebnerBasis` / the multivariate division
  algorithm).
- Inputs are `Expand`ed and (inexact leaves) rationalised before division.

**Options** (as for `GroebnerBasis`):

| Option               | Default              | Supported                                          |
| -------------------- | -------------------- | -------------------------------------------------- |
| `MonomialOrder`      | `Lexicographic`      | `Lexicographic`, `DegreeLexicographic`, `DegreeReverseLexicographic`, or an integer weight matrix. Non-well-founded (`Negative*`) orders fall back to `Lexicographic`. |
| `CoefficientDomain`  | `RationalFunctions`  | `RationalFunctions` (the default; over the field `Q(params)`, so cofactors may carry parameter denominators), `Rationals` (over `Q`; requires every symbol to be a variable). |
| `Modulus`            | `0`                  | `Modulus -> p` (prime `p` in `[2, 2^31)`) reduces over `GF(p)` (`gfp_divmod`, `Lexicographic` / `DegreeReverseLexicographic`; not combined with parameters). |
| `ParameterVariables` | `{}`                 | Symbol or list of symbols treated as coefficient-field parameters. |
| `Tolerance`          | `0`                  | Only `0` is supported. |

**Deferred** (accepted, then the head declines with a
`PolynomialReduce::nimpl` note): `CoefficientDomain -> Integers` and
`InexactNumbers`, and a nonzero `Tolerance`.

```mathematica
In[1]:= PolynomialReduce[x^3 + y^3, {x^2 - y^2 - 1, x + 2 y - 7}, {x, y}]
Out[1]= {{x, 1 + y^2}, 7 - 2 y + 7 y^2 - y^3}

(* f is a linear combination of the divisors plus the remainder. *)
In[2]:= f == First[%].{x^2 - y^2 - 1, x + 2 y - 7} + Last[%] // Expand
Out[2]= True

(* Not a Gröbner basis -> nonzero remainder even though f is in the ideal. *)
In[3]:= PolynomialReduce[2 x^3 + y^3 + 3 y, {x^2 + y^2 - 1, x y - 2}, {x, y}]
Out[3]= {{2 x, -2 y}, 2 x - y + y^3}

(* Ideal membership: remainder modulo a Gröbner basis is 0. *)
In[4]:= gb = GroebnerBasis[{x^2 - y^3 - 5, y^2 - x^3 + 7}, {x, y}];
        PolynomialReduce[x^6 - 14 x^3 - y^4 + 49, gb, {x, y}][[2]]
Out[4]= 0

(* Variable order changes the result. *)
In[5]:= PolynomialReduce[x^5 + (x + y)^2, {x^2 - y^3 - 5, y^3 - x^3 + 7}, {y, x}]
Out[5]= {{0, 0}, x^2 + x^5 + 2 x y + y^2}

(* RationalFunctions (default): cofactors over Q(a) carry a in the denominator. *)
In[6]:= gb1 = GroebnerBasis[{a x^2 + 5 y - 1, 2 x + x y - y^2}, {x, y},
                            CoefficientDomain -> RationalFunctions];
        PolynomialReduce[a^2 x - x y + y^2 - 3, gb1, {x, y}][[1]]
Out[6]= {1/(4 a), a/4 - y/(4 a)}    (* signs relative to gb1's normalisation *)

(* Reduce over GF(7). *)
In[7]:= PolynomialReduce[x^3, {x^2 + 1}, {x}, Modulus -> 7]
Out[7]= {{x}, 6 x}
```

## Eliminate
Eliminates variables between a set of simultaneous equations.
- `Eliminate[eqns, vars]`
- `Eliminate[eqns, var]` — single-symbol shorthand for `Eliminate[eqns, {var}]`.

`eqns` is a list `{eq1, eq2, ...}` or a conjunction `eq1 && eq2 && ...` of
`lhs == rhs` equations (a single `Equal[lhs, rhs]` is also accepted).
Surviving constraints are returned as balanced `Equal[posPart, negPart]`
equations, joined with `&&` when more than one remains; the empty
elimination ideal collapses to `True` and an inconsistent system to
`False`.

**Features**:
- `Protected`.
- Drives the lex-order Buchberger engine (`GroebnerBasis`) with an
  elimination block: `vars` are placed in the high-priority block so
  that any polynomial mentioning them is filtered out before output.
- A single-layer principal-branch inverse-function pre-pass handles
  equations of the form `f[u] == v` (or `v == f[u]`) for `f` in
  `{Sin, Cos, Tan, Sinh, Cosh, Tanh, Exp, Log}` and their inverses,
  rewriting to `u == InverseF[v]`.  Whenever a rewrite fires Eliminate
  emits the `Eliminate::ifun` diagnostic, signalling that the principal
  branch was used and some solutions may be missed — use `Reduce` for
  complete solution information.
- An algebraisation pre-pass handles inputs containing
  `Power[base, p/q]` with `q > 1` (i.e. `Sqrt[...]`, `(...)^(1/3)`,
  `(...)^(3/2)`, `1/Sqrt[...]`, ...) whose base mentions an elim
  variable.  Each unique base gets a fresh auxiliary symbol `$elN$`
  with the constraint `Power[$elN$, L] == base` (L the LCM of all
  rational-exponent denominators for that base), the power is rewritten
  as `Power[$elN$, p*L/q]`, and `$elN$` is appended to the elim list.
  Nested radicals (e.g. `Sqrt[x + Sqrt[x]]`) work in one shot because
  the constraint substitutes inner aux vars as well.  When the pass
  fires, `Eliminate::alg` is emitted to flag that the returned
  polynomial relation is the cross-multiplied generic consequence —
  sign / branch information may be lost.
- A transcendental algebraisation pre-pass (the general sibling of the
  radical pass) handles elim variables appearing inside circular /
  hyperbolic trig, exponential (`b^x`), or logarithmic functions — even
  in multiple places, where the single-layer inverse pre-pass cannot
  help.  Each equation is first expanded (`TrigExpand` for trig;
  `b^(p+q) -> b^p b^q`, `b^(k m) -> (b^m)^k` for exp; `Log[a b] ->
  Log[a]+Log[b]`, `Log[a^n] -> n Log[a]` for log) so every kernel lands
  on an *atomic* argument; thus `Sin[x]` and `Sin[3x]`, or `Exp[x]` and
  `Exp[2x]`, collapse onto one shared aux.  Exponentials whose exponents
  differ only by a *rational* factor of a shared monomial are commensurate
  and also collapse: for each `(base, monomial)` group the atomic kernel is
  `base^((1/L) m)` with `L` the LCM of the coefficient denominators, so
  `E^(x/3)` and `E^(x/2)` both become powers of `E^(x/6)`
  (`E^(x/3) = (E^(x/6))^2`, `E^(x/2) = (E^(x/6))^3`) — without this they
  register as independent kernels and Gate A rejects the system as `nlin`.  Each `Sin[θ]`/`Cos[θ]`
  (resp. `Sinh`/`Cosh`) becomes a fresh aux pair `$tsN$`,`$tcN$` with
  the Pythagorean constraint `s^2 + c^2 == 1` (circular) /
  `c^2 - s^2 == 1` (hyperbolic); each `b^θ` / `Log[θ]` becomes a single
  algebraically-free aux `$teN$` / `$tlN$`.  The auxes join the elim
  list and `Eliminate::ifun` is emitted (branch / sign information may be
  lost).  Two conservative soundness gates bail to `Eliminate::nlin`
  rather than emit a wrong relation: when an elim variable is shared
  across two distinct kernels (e.g. `Sin[x]` and `Sin[x*y]`, or `Sin[x]`
  and `Exp[x]`), or when it still appears as a bare polynomial atom after
  substitution (e.g. `x` alongside `Sin[x]`).  This is what lets the
  integration-by-substitution pattern
  `Eliminate[{Dt[y] == Cos[x] Sqrt[1-Sin[x]] Dt[x], u == Sin[x],
  Dt[u] == Cos[x] Dt[x]}, {x, Dt[x]}]` return `Dt[y]^2 == (1-u) Dt[u]^2`.
- An inverse-function substitution pre-pass handles the mirror shape,
  where an *inverse* trig / hyperbolic function of an elim variable sits
  buried inside a product (so the single-layer inverse pre-pass cannot
  reach it) but a *defining* equation `M == ArcF[x]` is present (`M`
  elim-free, `x` a single elim symbol, `ArcF` in `{ArcSin, ArcCos,
  ArcTan, ArcSinh, ArcCosh, ArcTanh}`).  The pass propagates
  `ArcF[x] -> M` through the whole system, rewrites the companion base
  `c + s·x^2 = co[M]^(2·co_sign)` to a power of the co-function of the
  *main* angle for **any** rational exponent — the half-integer (radical)
  case `Sqrt[1-x^2] -> Cos[M]` (ArcSin; `-> Sin[M]` ArcCos, `Sqrt[1+x^2]
  -> Cosh[M]` ArcSinh, `Sqrt[x^2-1] -> Sinh[M]` ArcCosh, `-> Sec[M]`/
  `Sech[M]` ArcTan/ArcTanh) **and** the integer (rational-denominator) case
  `1/(1+x^2) -> Cos[M]^2` (ArcTan), `1/(1-x^2) -> Cosh[M]^2` (ArcTanh) — then
  pins `x == F[M]` (F the forward function).  Because the co-function is a
  trig call of the retained variable `M`, any factor of it left after
  elimination is a main-variable monomial that is divided out cleanly (or
  cancels directly in the Groebner step) — so the u-substitution shapes
  `Eliminate[{Dt[y] == x ArcSin[x]/Sqrt[1-x^2] Dt[x], u == ArcSin[x],
  Dt[u] == Dt[x]/Sqrt[1-x^2]}, {x, Dt[x]}]` and
  `Eliminate[{Dt[y] == ArcTan[x]/(1+x^2) Dt[x], u == ArcTan[x],
  Dt[u] == Dt[x]/(1+x^2)}, {x, Dt[x]}]` return `u Sin[u] Dt[u] == Dt[y]`
  and the fully-reduced `u Dt[u] == Dt[y]` respectively (with
  `Eliminate::ifun`).  A `Log` defining
  equation `M == Log[x]` is handled by the same pass (forward function
  `Exp`, no companion radical) — but only as a *fallback*, when `x` also
  occurs as a genuine polynomial atom (e.g. `1/x`), which is exactly the
  case the forward log pass cannot capture (Gate B).  Rather than
  substituting `x -> E^M`, it pins `x == E^M` so `E^M` stays a single
  main-variable atom (the Groebner term reader treats a `Power[base, exp]`
  with symbolic exponent as one opaque indeterminate); when `x` appears
  solely inside logs the forward pass gives a tidier answer and is used
  instead.  So `Eliminate[{Dt[y] == Log[x]/x Dt[x], u == Log[x],
  Dt[u] == -Dt[x]/x^2}, {x, Dt[x]}]` returns `u Dt[u] E^u + Dt[y] == 0`
  (i.e. `-Dt[y] == E^u u Dt[u]`).  This is the pass that unlocks
  `Integrate[x ArcSin[x]/Sqrt[1-x^2], x]` and its inverse-trig /
  inverse-hyperbolic kin — and now `Integrate[Log[x]/x, x]`-style log
  substitutions — through `DerivativeDivides`.
- Equations are normalised through `Numerator[Together[lhs - rhs]]`
  before Buchberger, clearing any `Power[t, -k]` denominators
  introduced by the algebraisation pre-passes.  Surviving basis
  polynomials are stripped of any common monomial factor so the
  reported equation matches the primitive form (no spurious factor
  of a main variable from a `u == Sqrt[...]` equation).
- Equations that remain non-polynomial after both pre-passes (e.g.
  `Sin[x*y] == ...` with `x` in the elim set) return the expression
  unevaluated with `Eliminate::nlin`.
- Main-block variables (the implicit free parameters) are discovered
  automatically by walking the equations; mathematical constants such
  as `Pi`, `E`, and `EulerGamma` flow through as parameter symbols.
  Function-shaped subexpressions whose head is not one of
  `{Plus, Times, Power, List, And, Equal, Or, Not}` are treated as
  single polynomial atoms (so `Dt[y]`, `f[a]`, `Sin[k]` etc. are
  whole-expression variables rather than products of their heads and
  arguments).

```mathematica
In[1]:= Eliminate[{x == 2 + y, y == z}, y]
Out[1]= 2 + z == x

In[2]:= Eliminate[{f == x^5 + y^5, a == x + y, b == x y}, {x, y}]
Out[2]= a^5 + 5 a b^2 == 5 a^3 b + f                  (* equiv to f == a^5 - 5 a^3 b + 5 a b^2 *)

In[3]:= Eliminate[(2 x + 3 y + 4 z == 1) && (9 x + 8 y + 7 z == 2), z]
Out[3]= 22 x + 11 y == 1

In[4]:= Eliminate[{x^2 + y^2 + z^2 == 1, x - y + z == 2,
                   x^3 - y^2 == z + 1}, z]
Out[4]= 27 + 4 x^2 + 8 x^4 + 4 x^5 + 4 x^6 == 18 x + 28 x^3
        && 12 + 2 x + 5 x^2 + y == 8 x^3 + 4 x^4 + 2 x^5

In[5]:= Eliminate[{x^2 + y^2 + z^2 == 1, x - y + z == 2,
                   x^3 - y^2 == z + 1}, {y, z}]
Out[5]= 27 + 4 x^2 + 8 x^4 + 4 x^5 + 4 x^6 == 18 x + 28 x^3

In[6]:= Eliminate[{x - a == 0, x - b == 0}, x]  (* common-root condition *)
Out[6]= b == a

In[7]:= Eliminate[1 == 2, x]                    (* inconsistent -> False *)
Out[7]= False

In[8]:= Eliminate[x + y == 0, y]                (* solvable -> True *)
Out[8]= True

In[9]:= Eliminate[{u == Sqrt[x^2 + 1], v == 1/Sqrt[x^2 + 1]}, x]
Out[9]= u v == 1                                (* shared base, single aux *)

In[10]:= Eliminate[{Dt[y] == x^3/Sqrt[x^2 + 1] Dt[x], u == x^2 + 1,
                    Dt[u] == 2 x Dt[x]}, {Dt[x], x}]
Out[10]= Dt[u]^2 + u^2 Dt[u]^2 == 2 u Dt[u]^2 + 4 u Dt[y]^2
         (* algebraisation pre-pass; equivalent to the reference's
            u^2 Dt[u]^2 + u(-2 Dt[u]^2 - 4 Dt[y]^2) == -Dt[u]^2 *)

In[11]:= Eliminate[{u == x^(1/3), v == x^(2/3)}, x]
Out[11]= v == u^2                               (* exponent LCM = 3 *)
```

**Features**:
- `Protected`.
- Nests multiplications instead of using powers (e.g., $a + x(b + c x)$ instead of $a + bx + cx^2$).
- Identifies variables using `Variables` if not explicitly specified.
- Issues an error and returns unevaluated if the expression is not a polynomial or rational function in the target variables.

```mathematica
In[1]:= HornerForm[11 x^3 - 4 x^2 + 7 x + 2]
Out[1]= 2 + x (7 + x (-4 + 11 x))

In[2]:= HornerForm[a + b x + c x^2, x]
Out[2]= a + x (b + c x)

In[3]:= HornerForm[(11 x^3 - 4 x^2 + 7 x + 2)/(x^2 - 3 x + 1)]
Out[3]= (2 + x (7 + x (-4 + 11 x))) / (1 + x (-3 + x))

In[4]:= HornerForm[1 + x^a, x]
HornerForm::poly: 1+x^a is not a polynomial.
Out[4]= HornerForm[1 + x^a, x]
```


## FLINT` context (direct FLINT kernels)

The `FLINT`` context exposes FLINT's kernels directly at the REPL, in parallel
with the transparent acceleration of the corresponding `System`` builtins. Each
is `Protected`, carries a docstring (`?FLINT`Factor`), and its symbol is interned
in `sym_names.c`. All require a build with FLINT (`USE_FLINT`); without it — or
when the argument is out of FLINT's scope (non-polynomial / non-rational
coefficients, a symbolic matrix, a symbolic function argument or a pole) — the
call is returned unevaluated rather than approximated.

Polynomial kernels (`fmpq_mpoly`):

- `` FLINT`PolynomialGCD[a, b] `` — the **monic** GCD of `a` and `b` over the
  rationals (`fmpq_mpoly_gcd`), multivariate. Note this is FLINT's monic
  representative (e.g. `` FLINT`PolynomialGCD[2 x - 1, 4 x^2 - 1] `` is
  `x - 1/2`), which differs from `PolynomialGCD`'s primitive-over-Z form
  (`2 x - 1`).
- `` FLINT`Resultant[a, b, x] `` — the resultant eliminating `x`
  (`fmpq_mpoly_resultant`); other free symbols are treated as coefficients.
- `` FLINT`Factor[p] `` — irreducible factorisation over Q
  (`fmpq_mpoly_factor`), returned as `Times[const, factor^exp, …]`,
  multivariate. (`System`Factor` on a **univariate** polynomial over Z is
  transparently routed through FLINT's `fmpz_poly_factor`, giving primitive
  positive-leading factors + integer content — the reference's convention. This
  also fixes a bug where a univariate polynomial with a bignum coefficient came
  back unfactored. Transparent *multivariate* Factor is not rewired: FLINT's
  factor sign-form / ordering differs from the classical output.)
- `` FLINT`FactorSquareFree[p] `` — squarefree factorisation
  (`fmpq_mpoly_factor_squarefree`).

Exact linear algebra (integer / rational matrices, `fmpq_mat`); each also
transparently accelerates the corresponding `System`` builtin on all-rational
matrices (the results — determinant, inverse, RREF, unique solution, rank — are
unique / basis-independent, so they match the classical output):

- `` FLINT`Det[m] `` — exact determinant via `fmpq_mat_det` (see *Linear
  Algebra → Det*).
- `` FLINT`Inverse[m] `` — exact inverse via `fmpq_mat_inv`; unevaluated for a
  singular or non-rational matrix.
- `` FLINT`RowReduce[m] `` — reduced row echelon form via `fmpq_mat_rref`.
- `` FLINT`LinearSolve[m, b] `` — square nonsingular solve via `fmpq_mat_solve`
  (vector or matrix RHS); unevaluated otherwise.
- `` FLINT`MatrixRank[m] `` — rank via `fmpq_mat_rref`.

Arbitrary-precision special functions (arb/acb ball arithmetic; result
precision follows the arguments, machine precision for exact input):

- `` FLINT`Zeta[s] `` — Riemann ζ (`acb_dirichlet_zeta`).
- `` FLINT`HurwitzZeta[s, a] `` — Hurwitz ζ (`acb_dirichlet_hurwitz`).
- `` FLINT`PolyGamma[n, z] `` — ψ⁽ⁿ⁾, `n = 0` is digamma (`acb_polygamma`).
- `` FLINT`StieltjesGamma[n] `` / `` [n, a] `` — Stieltjes constants
  (`acb_dirichlet_stieltjes`); non-negative integer `n`.

```mathematica
In[1]:= FLINT`Factor[x^2 - 2 x y + y^2]
Out[1]= (x - y)^2

In[2]:= FLINT`Resultant[x^2 - a, x^2 - b, x]
Out[2]= a^2 - 2 a b + b^2

In[3]:= FLINT`Det[{{1/2, 1/3}, {1/5, 1/7}}]
Out[3]= 1/210

In[4]:= FLINT`Zeta[2] - N[Pi^2/6]
Out[4]= 0.
```
