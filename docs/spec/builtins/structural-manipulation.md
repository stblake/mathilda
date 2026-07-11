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
  The spec must be a permutation of `{1, ..., r}` where `r` is the depth of
  `list`. A repeated index (e.g. `{1, 1}`) extracts the corresponding
  diagonal.

**Features**:
- `Protected`.
- Works only on rectangular arrays.
- `Transpose[m, {1, 1}]` extracts the diagonal of a square matrix.

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
  does not change the shape of `vec` (matches the Mathematica convention).
- For symbolic entries, `Conjugate[x]` is left wrapped around `x`.
- Works on higher-rank tensors with the same `spec` semantics as
  `Transpose`.

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
- `Flatten[list, n]`: Flattens up to level `n`.
- `Flatten[list, n, h]`: Flattens subexpressions with head `h`.

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

```mathematica
In[1]:= Split[{a, a, a, b, b, a, a, c}]
Out[1]= {{a, a, a}, {b, b}, {a, a}, {c}}

In[2]:= Split[{1, 2, 3, 4, 3, 2, 1}, Less]
Out[2]= {{1, 2, 3, 4}, {3}, {2}, {1}}
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
- Canonical order:
    - Real numbers by numerical value.
    - Complex numbers by real part, then imaginary part magnitude.
    - Strings in dictionary order (lowercase before uppercase).
    - Symbols by name.
    - Expressions by length, then head, then parts depth-first.
- Polynomial order: `x^n` sorts relative to `x`.
- Numeric coefficient stripping: when a `Times` term has a leading numeric factor (Integer, Real, BigInt, MPFR, `Rational[n,d]`, `Complex[re,im]`, or a radical such as `Sqrt[2] = Power[2, 1/2]`), that factor is ignored when computing the term's main factor. As a result, `1 + x^2 + Sqrt[2] x` is canonicalised to `1 + Sqrt[2] x + x^2`, matching Mathematica's main-factor-first ordering.
- Custom ordering function `p` can return `1`, `0`, `-1`, `True`, or `False`.

```mathematica
In[1]:= Sort[{d, b, c, a}]
Out[1]= {a, b, c, d}

In[2]:= Sort[{Pi, E, 2, 3, 1, Sqrt[2]}, Less]
Out[2]= {1, Sqrt[2], 2, E, 3, Pi}
```

## Union
Gives a sorted list of all distinct elements from one or more expressions.
- `Union[list]`
- `Union[list1, list2, ...]`
- `Union[..., SameTest -> test]`

**Features**:
- `Flat`, `OneIdentity`, `Protected`, `ReadProtected`.
- All expressions must have the same head.
- Result has the same head as the inputs.

```mathematica
In[1]:= Union[{1, 2, 1, 3, 6, 2, 2}]
Out[1]= {1, 2, 3, 6}

In[2]:= Union[{a, b, a, c}, {d, a, e, b}, {c, a}]
Out[2]= {a, b, c, d, e}
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

```mathematica
In[1]:= Min[9, 2]
Out[1]= 2

In[2]:= Min[{4, 1, 7, 2}]
Out[2]= 1

In[3]:= Max[Infinity, 5]
Out[3]= Infinity
```

## MinMax
`MinMax[list]` gives `{Min[list], Max[list]}` in one shot — the range of the
data, handy for plot bounds. Over an association it uses the values. Delegates
to `Min` / `Max`, so it inherits their numeric handling.

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
- Compound expressions have depth `1 + Max(depths of arguments)`.
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
- Lists subexpressions in post-order (depth-first), resulting in lexicographic ordering of indices.

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
- Implements an efficient binary-splitting algorithm for distributing products and repeated squaring for powers.
- `Expand[expr, patt]` leaves unexpanded any parts of `expr` that are free of the pattern `patt`.

```mathematica
In[1]:= Expand[(x+3)(x+2)]
Out[1]= 6 + 5 x + x^2

In[2]:= Expand[(x+y)^2 (x-y)^2]
Out[2]= x^4 - 2 x^2 y^2 + y^4

In[3]:= Expand[(x+1)^2 + (y+1)^2, x]
Out[3]= 1 + 2x + x^2 + (1+y)^2
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

## Coefficient
Gives the coefficient of a specific form in a polynomial.
- `Coefficient[expr, form]`
- `Coefficient[expr, form, n]`

**Features**:
- `Protected`, `Listable`.
- `Coefficient[expr, form, 0]` picks out terms that do NOT contain `form`.
- Works whether or not `expr` is explicitly given in expanded form (it automatically expands internally).
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

```mathematica
In[1]:= Collect[a x + b y + c x + d y, x]
Out[1]= b y + d y + (a + c) x

In[2]:= Collect[a x^4 + b x^4 + 2 a^2 x - 3 b x + x - 7, x]
Out[2]= -7 + (1 + 2 a^2 - 3 b) x + (a + b) x^4

In[3]:= Collect[a Sqrt[x] + Sqrt[x] + x^(2/3) - c x + 3x - 2b x^(2/3) + 5, x]
Out[3]= 5 + (1 + a) Sqrt[x] + (1 - 2b) x^(2/3) + (3 - c) x
```

## FactorSquareFree
Pulls out any multiple factors in a polynomial using Yun's algorithm.
- `FactorSquareFree[poly]`

**Features**:
- `Listable`, `Protected`.
- Automatically threads over lists, as well as equations, inequalities and logic functions.
- Works on both univariate and multivariate polynomials.
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
Out[4]= {(-1 + x)^2 (1 + x), (-1 + x^2)^2 (1 + x^2)}
```

## Factor
Factors a polynomial over the integers.
- `Factor[poly]`

**Features**:
- `Listable`, `Protected`.
- **FLINT acceleration** (when built with FLINT): a univariate polynomial over Z is factored via `fmpz_poly_factor`, and a genuine multivariate polynomial over Q via `fmpq_mpoly_factor`, both from a guarded fast path at the top of `builtin_factor` (plain single-argument form only). This turns the exponential-in-variable-count classical multivariate factoriser — which hangs on inputs like `Factor[x^99 - y^99]` (> 20 s) — into a few-millisecond call. Each irreducible factor is normalised to a positive leading coefficient (highest total degree, deglex tie-break — Mathematica's convention), with the discarded sign folded into a separated rational content, so e.g. `Factor[y^2 - x^2]` → `-(x - y)(x + y)`. Inputs that are not a polynomial over Q (a denominator, fractional/symbolic exponent, or a non-polynomial head such as `Sqrt`/`Sin`) fall through to the classical path unchanged; `Extension`/`GaussianIntegers`/`Modulus` option forms travel their own branch and are untouched. See also the `` FLINT`Factor `` builtin.
- When given a rational expression, first resolves dependencies over `Together` before factoring.
- Uses exact root isolation (Rational Root Theorem limits) and binomial descents structured identically to Zassenhaus recombination, evaluating combinations exact and memory safe.
- Threads natively across lists, logic structures, and numeric groupings perfectly.
- Bivariate inputs whose leading coefficient (in some variable) is the constant `-1` are handled via Wang's leading-coefficient correction, Stage 1: the input is pre-negated to make it monic, the existing monic Hensel pipeline runs on `-P`, and the overall sign is absorbed into the highest-degree factor via `Expand`.  This unlocks inputs of shape `Factor[(1 - x^k)(x - y^m)]` (and similar non-monic cases with constant `±1` LC) that previously fell back to the legacy linear-trial-division loop.
- Bivariate inputs whose leading coefficient (in some variable) is a non-unit integer constant `a` (with `|a| > 1`) are handled via Wang's leading-coefficient correction, Stage 2: the monic substitution `Q(x, y) = a^(d-1) · P(x/a, y)` makes the lift's input monic in x with integer coefficients.  After lifting `Q = G_1 · ... · G_r` via the existing pipeline, the true factors of `P` are recovered as `F_i = G_i(a·x, y) / cont_Z(G_i(a·x, y))`, where the integer content collects exactly the share of `a^(d-1)` redistributed into G_i by the substitution.  Stages 1 and 2 compose, so inputs with negative non-unit LC (e.g. `lc_x(P) = -6`) also enter the structured pipeline.  This unlocks inputs like `Factor[Expand[(2x+3y)(3x+5y)]]`, `Factor[2 a^2 - 5 a b + 3 b^2]`, and three-factor non-monic forms whose LCs are constant in y (or constant in x).
- Bivariate inputs whose leading coefficient (in some variable) is a non-constant polynomial in the other variable are handled via Wang's leading-coefficient correction, Stage 3 (predicted-LC two-factor Hensel).  When `lc_x(P)(y) = A(y)`, Mathilda factors `A` over `Z[y]`, finds `α` with `A(α) = +1` so the squarefree univariate image `P(x, α)` factors into monic Z[x] pieces `u`, `v`, then enumerates distributions of `A`'s irreducible factors between two predicted leading coefficients `q_u, q_v` (with `q_u · q_v = A` and `q_u(α) = q_v(α) = +1`).  The Hensel iteration is modified so each `Δu` correction has its leading-x coefficient PINNED to the y^k coefficient of `q_u`, keeping `lc_x(U)(y) = q_u(y)` invariant across the lift.  This unlocks inputs like `Factor[Expand[(xy+1)(xy+2)]]`, `Factor[Expand[((y²+1)x+1)(x+3)]]`, `Factor[Expand[((y+1)x+1)((y+1)x+2)]]`.  MVP scope: r = 2 (two univariate factors), both monic, `|cont(A)| = 1`, and inputs with non-trivial monomial content fall through so `heuristic_factor`'s Phase 0 path produces the canonical fully-factored form.

```mathematica
In[1]:= Factor[1 + 2x + x^2]
Out[1]= (1 + x)^2

In[2]:= Factor[x^10 - 1]
Out[2]= (-1 + x) (1 + x) (1 - x + x^2 - x^3 + x^4) (1 + x + x^2 + x^3 + x^4)

In[3]:= Factor[x^10 - y^10]
Out[3]= (x - y) (x + y) (x^4 - x^3 y + x^2 y^2 - x y^3 + y^4) (x^4 + x^3 y + x^2 y^2 + x y^3 + y^4)

In[4]:= Factor[2x^3 y - 2a^2 x y - 3a^2 x^2 + 3a^4]
Out[4]= (a - x) (a + x) (3 a^2 - 2 x y)

In[5]:= Factor[(x^3 + 2x^2)/(x^2 - 4y^2) - (x + 2)/(x^2 - 4y^2)]
Out[5]= ((-1 + x) (1 + x) (2 + x)) / ((x - 2 y) (x + 2 y))
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
- `Protected`.
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
Out[7]= 2 (-3 + x^2 + 2 y z) (2 + y + 2 z + y z)

In[8]:= FactorTerms[f, {x, y}]
Out[8]= 2 (1 + z) (2 + y) (-3 + x^2 + 2 y z)
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
  leading monomial — matching Mathematica's output convention.
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
| `MonomialOrder`      | `Lexicographic`  | `Lexicographic`, `DegreeReverseLexicographic`, `EliminationOrder`, or an integer **weight matrix** `{{...}, ...}` (one column per variable) defining a custom term order |
| `CoefficientDomain`  | `Rationals`      | `Rationals`, `Automatic`, `Integers` (strong Gröbner basis over Z), `RationalFunctions` (basis over the field Q(params)), `Polynomials[x, ...]` (named symbols become coefficient-ring parameters), `InexactNumbers[p]` (monic basis with `p`-digit real coefficients) |
| `Method`             | `Automatic`      | `"Buchberger"`, `"GroebnerWalk"`, `Automatic`        |
| `Sort`               | `False`          | `True` reverses the main-variable list before the basis is computed (matching Mathematica's empirical behaviour). |
| `ParameterVariables` | (unset)          | Symbol or list of symbols treated as coefficient-ring parameters; with no positional `vars` argument the main-variable list is auto-derived from the polynomials. |

A **weight-matrix** `MonomialOrder` must be a rectangular list of integer
rows with exactly one column per main variable, and must define a valid
term order: full column rank (so distinct monomials never tie) and a
positive first non-zero entry in every column (so each variable exceeds
1, making the order well-founded).  Matrices that fail either condition,
or that are used together with elimination/parameter variables, emit
`GroebnerBasis::nimpl` and fall back to `Lexicographic`.

Settings outside this table emit a `GroebnerBasis::nimpl` diagnostic to
`stderr` and fall back to the default.  Deferred (not implemented in this
release):

- `MonomialOrder -> DegreeLexicographic`.
- `Modulus -> n` (modular Gröbner bases): emits a
  `GroebnerBasis::modnotimpl` diagnostic and falls back to the rational
  basis.

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
  `ArcF[x] -> M` through the whole system, rewrites the companion radical
  to the co-function of the *main* angle (`Sqrt[1-x^2] -> Cos[M]` for
  ArcSin, `-> Sin[M]` for ArcCos, `Sqrt[1+x^2] -> Cosh[M]` for ArcSinh,
  `Sqrt[x^2-1] -> Sinh[M]` for ArcCosh, `-> Sec[M]`/`Sech[M]` for
  ArcTan/ArcTanh), and pins `x == F[M]` (F the forward function).  Because
  the co-function is a trig call of the retained variable `M`, any factor
  of it left after elimination is a main-variable monomial that is divided
  out cleanly — so the u-substitution shape
  `Eliminate[{Dt[y] == x ArcSin[x]/Sqrt[1-x^2] Dt[x], u == ArcSin[x],
  Dt[u] == Dt[x]/Sqrt[1-x^2]}, {x, Dt[x]}]` returns
  `u Sin[u] Dt[u] == Dt[y]` (with `Eliminate::ifun`).  A `Log` defining
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
         (* algebraisation pre-pass; equivalent to Mathematica's
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
  positive-leading factors + integer content — Mathematica's convention. This
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
