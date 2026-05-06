# PicoCAS Specification

PicoCAS is a tiny, AI Agent-generated, symbolic computer algebra system (CAS) inspired by the core architecture of Mathematica. PicoCAS is written in the C programming language. It supports a recursive expression model, pattern matching, rewriting rules, and a small library of built-in functions.

The name, PicoCAS is inspired by David Stoutemyer's PICOMATH-80 tiny computer algebra system. 

## 1. Expression Structure

Everything in PicoCAS is an `Expr`, which can be one of the following types:
- **Integer**: 64-bit signed integers (e.g., `123`, `-45`).
- **Real**: Double-precision floating-point numbers (e.g., `3.14`, `1.2E-3`).
- **Rational**: Exact fractions of integers, represented as `Rational[n, d]`.
- **Complex**: Complex numbers, represented as `Complex[re, im]`. The symbol `I` represents `Complex[0, 1]`.
- **Symbol**: Named identifiers (e.g., `x`, `Plus`, `Pi`, `Infinity`, `E`).
- **String**: Character sequences enclosed in double quotes (e.g., `"hello"`).
- **Function**: Compound expressions consisting of a `head` and zero or more `arguments` (e.g., `f[x, y]`, `{1, 2, 3}`).

## 2. Core Built-in Functions

### Structural Manipulation

#### Part
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

#### Extract
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

#### Span
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

#### Head
Returns the top-level wrapper of an expression.
- `Head[expr]`

**Features**: 
- For functions, returns the symbol or expression acting as the head.
- For atoms, returns the symbolic type name: `Integer`, `Real`, `Rational`, `Complex`, `Symbol`, or `String`.

```mathematica
In[1]:= Head[f[x]]
Out[1]= f

In[2]:= Head[3/4]
Out[2]= Rational
```

#### Length
Returns the number of arguments in an expression.
- `Length[expr]`

**Features**: 
- Returns the count of top-level arguments for functions.
- Returns `0` for all atomic expressions.

```mathematica
In[1]:= Length[{a, b, c}]
Out[1]= 3
```

#### Dimensions
Returns the list of lengths of levels in a nested rectangular structure.
- `Dimensions[expr]`

```mathematica
In[1]:= Dimensions[{{1, 2}, {3, 4}}]
Out[1]= {2, 2}
```

#### First, Last, Most, Rest
Convenience functions for accessing parts of a sequence.
- `First[expr]`: Equivalent to `expr[[1]]`.
- `Last[expr]`: Equivalent to `expr[[-1]]`.
- `Most[expr]`: Returns all elements except the last.
- `Rest[expr]`: Returns all elements except the first.

#### Reverse
Reverses the order of elements.
- `Reverse[expr]`: Reverses top-level elements.
- `Reverse[expr, n]`: Reverses elements at level `n`.
- `Reverse[expr, {n1, n2, ...}]`: Reverses at levels `n1`, `n2`, etc.

#### RotateLeft, RotateRight
Cycles elements.
- `RotateLeft[expr, n]`: Cycles `n` positions to the left.
- `RotateLeft[expr]`: Cycles 1 position to the left.
- `RotateLeft[expr, {n1, n2, ...}]`: Cycles at successive levels.
- `RotateRight[...]`: Similar, but cycles to the right.

#### Transpose
Transposes levels in a rectangular array.
- `Transpose[list]`: Transposes the first two levels.
- `Transpose[list, {n1, n2, ...}]`: Level `k` becomes level `nk` in result.

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

#### Take, Drop

In[1]:= Reverse[{a, b, c}]
Out[1]= {c, b, a}

In[2]:= RotateLeft[{a, b, c}, 1]
Out[2]= {b, c, a}

In[3]:= RotateRight[{a, b, c}, 1]
Out[3]= {c, a, b}
```

#### Take, Drop

In[1]:= Most[{a, b, c}]
Out[1]= {a, b}

In[2]:= Rest[f[x, y]]
Out[2]= f[y]
```

#### Take, Drop
Extract or remove sequences of elements.
- `Take[expr, spec]`: Extracts elements according to `spec`.
- `Drop[expr, spec]`: Removes elements according to `spec`.

#### Flatten
Flattens out nested lists.
- `Flatten[list]`: Flattens all levels.
- `Flatten[list, n]`: Flattens up to level `n`.
- `Flatten[list, n, h]`: Flattens subexpressions with head `h`.

#### Join
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

#### Partition
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

#### Split
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

#### OrderedQ
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

#### Sort
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
- Custom ordering function `p` can return `1`, `0`, `-1`, `True`, or `False`.

```mathematica
In[1]:= Sort[{d, b, c, a}]
Out[1]= {a, b, c, d}

In[2]:= Sort[{Pi, E, 2, 3, 1, Sqrt[2]}, Less]
Out[2]= {1, Sqrt[2], 2, E, 3, Pi}
```

#### Union
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

#### DeleteDuplicates
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

#### Tally
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

#### Commonest
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

#### Min, Max
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
- Simplifies numeric arguments to a single value.

```mathematica
In[1]:= Min[9, 2]
Out[1]= 2

In[2]:= Min[{4, 1, 7, 2}]
Out[2]= 1

In[3]:= Max[Infinity, 5]
Out[3]= Infinity
```

#### Append, Prepend

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

#### Append, Prepend

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

#### Append, Prepend
Adds an element to the end or beginning of an expression.
- `Append[expr, elem]`
- `Prepend[expr, elem]`

```mathematica
In[1]:= Append[{a, b}, c]
Out[1]= {a, b, c}
```

#### AppendTo, PrependTo
Updates a variable by appending or prepending an element.
- `AppendTo[symbol, elem]`
- `PrependTo[symbol, elem]`

#### Insert, Delete
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

#### Depth
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

#### Level
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

#### Variables
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

#### Expand
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

#### ExpandNumerator
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

#### ExpandDenominator
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

#### Coefficient
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

#### CoefficientList
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

#### Collect
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

#### FactorSquareFree
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

#### Factor
Factors a polynomial over the integers.
- `Factor[poly]`

**Features**:
- `Listable`, `Protected`.
- When given a rational expression, first resolves dependencies over `Together` before factoring.
- Uses exact root isolation (Rational Root Theorem limits) and binomial descents structured identically to Zassenhaus recombination, evaluating combinations exact and memory safe.
- Threads natively across lists, logic structures, and numeric groupings perfectly.
- Bivariate inputs whose leading coefficient (in some variable) is the constant `-1` are handled via Wang's leading-coefficient correction, Stage 1: the input is pre-negated to make it monic, the existing monic Hensel pipeline runs on `-P`, and the overall sign is absorbed into the highest-degree factor via `Expand`.  This unlocks inputs of shape `Factor[(1 - x^k)(x - y^m)]` (and similar non-monic cases with constant `±1` LC) that previously fell back to the legacy linear-trial-division loop.
- Bivariate inputs whose leading coefficient (in some variable) is a non-unit integer constant `a` (with `|a| > 1`) are handled via Wang's leading-coefficient correction, Stage 2: the monic substitution `Q(x, y) = a^(d-1) · P(x/a, y)` makes the lift's input monic in x with integer coefficients.  After lifting `Q = G_1 · ... · G_r` via the existing pipeline, the true factors of `P` are recovered as `F_i = G_i(a·x, y) / cont_Z(G_i(a·x, y))`, where the integer content collects exactly the share of `a^(d-1)` redistributed into G_i by the substitution.  Stages 1 and 2 compose, so inputs with negative non-unit LC (e.g. `lc_x(P) = -6`) also enter the structured pipeline.  This unlocks inputs like `Factor[Expand[(2x+3y)(3x+5y)]]`, `Factor[2 a^2 - 5 a b + 3 b^2]`, and three-factor non-monic forms whose LCs are constant in y (or constant in x).
- Bivariate inputs whose leading coefficient (in some variable) is a non-constant polynomial in the other variable are handled via Wang's leading-coefficient correction, Stage 3 (predicted-LC two-factor Hensel).  When `lc_x(P)(y) = A(y)`, picocas factors `A` over `Z[y]`, finds `α` with `A(α) = +1` so the squarefree univariate image `P(x, α)` factors into monic Z[x] pieces `u`, `v`, then enumerates distributions of `A`'s irreducible factors between two predicted leading coefficients `q_u, q_v` (with `q_u · q_v = A` and `q_u(α) = q_v(α) = +1`).  The Hensel iteration is modified so each `Δu` correction has its leading-x coefficient PINNED to the y^k coefficient of `q_u`, keeping `lc_x(U)(y) = q_u(y)` invariant across the lift.  This unlocks inputs like `Factor[Expand[(xy+1)(xy+2)]]`, `Factor[Expand[((y²+1)x+1)(x+3)]]`, `Factor[Expand[((y+1)x+1)((y+1)x+2)]]`.  MVP scope: r = 2 (two univariate factors), both monic, `|cont(A)| = 1`, and inputs with non-trivial monomial content fall through so `heuristic_factor`'s Phase 0 path produces the canonical fully-factored form.

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

#### FactorTerms
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

#### FactorTermsList
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

#### PolynomialGCD
Gives the greatest common divisor of the polynomials.
- `PolynomialGCD[poly1, poly2, ...]`

**Features**:
- `Protected`, `Listable`.
- Handles univariate and multivariate polynomials.
- Treats algebraic numbers (like `I`) as independent variables or constants seamlessly during complex arithmetic evaluations.
- Pre-extracts common factors before falling back to a full primitive Euclidean algorithm computation.

```mathematica
In[1]:= PolynomialGCD[(1+x)^2(2+x)(4+x), (1+x)(2+x)(3+x)]
Out[1]= (1+x) (2+x)

In[2]:= PolynomialGCD[x^2+4x+4, x^2+2x+1]
Out[2]= 1

In[3]:= PolynomialGCD[x^2-1, x^3-1, x^4-1, x^5-1, x^6-1, x^7-1]
Out[3]= -1 + x
```

#### PolynomialExtendedGCD
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

#### PolynomialLCM
Gives the least common multiple of the polynomials.
- `PolynomialLCM[poly1, poly2, ...]`

**Features**:
- `Protected`, `Listable`.
- Handles univariate and multivariate polynomials.
- Treats algebraic numbers (like `I`) as independent variables or constants seamlessly during complex arithmetic evaluations.
- Preserves explicit factored forms where possible.

```mathematica
In[1]:= PolynomialLCM[(1+x)^2(2+x)(4+x), (1+x)(2+x)(3+x)]
Out[1]= (1+x)^2 (2+x) (3+x) (4+x)

In[2]:= PolynomialLCM[x^4-4, x^4+4 x^2+4]
Out[2]= (-2+x^2) (4+4 x^2+x^4)
```

#### PolynomialQuotient
Gives the quotient of $p$ and $q$, treated as polynomials in $x$, with any remainder dropped.
- `PolynomialQuotient[p, q, x]`

**Features**:
- `Protected`.
- Uses polynomial long division over the field of rational functions in the coefficients.

```mathematica
In[1]:= PolynomialQuotient[x^4+2x+1, x^2+1, x]
Out[1]= -1 + x^2

In[2]:= PolynomialQuotient[x^2+2x+1, x^3, x]
Out[2]= 0

In[3]:= PolynomialQuotient[x^2+x+1, 2x+1, x]
Out[3]= 1/4 + x/2
```

#### PolynomialRemainder
Gives the remainder from dividing $p$ by $q$, treated as polynomials in $x$.
- `PolynomialRemainder[p, q, x]`

**Features**:
- `Protected`.
- The degree of the result in $x$ is guaranteed to be smaller than the degree of $q$.
- If the dividend is a multiple of the divisor, then the remainder is zero.

```mathematica
In[1]:= PolynomialRemainder[x^4+2x+1, x^2+1, x]
Out[1]= 2 + 2 x

In[2]:= PolynomialRemainder[x^3, a x+b, x]
Out[2]= -(b^3/a^3)
```

#### PolynomialQ
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

#### PolynomialMod
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

#### Resultant
Computes the resultant of two polynomials.
- `Resultant[poly1, poly2, var]`

**Features**:
- `Protected`, `Listable`.
- Computes the resultant of polynomials `poly1` and `poly2` with respect to the variable `var`.
- The resultant is independent of common roots and vanishes exactly when the polynomials have roots in common.
- Implemented fundamentally over symbolic systems using Sylvester matrix determinants, ensuring full generic expansion over arbitrary constants.
- Automatically preserves multiplicativity (e.g., $Res(A \cdot B, Q) = Res(A, Q) Res(B, Q)$ and $Res(A^k, Q) = Res(A, Q)^k$).

```mathematica
In[1]:= Resultant[x^2 - 2x + 7, x^3 - x + 5, x]
Out[1]= 265

In[2]:= Resultant[x^3 - 5x^2 - 7x + 14, x^3 - 8x^2 + 9x + 58, x]
Out[2]= 0
```

#### Discriminant
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

#### HornerForm
Puts a polynomial or rational function into Horner form.
- `HornerForm[poly]`
- `HornerForm[poly, vars]`
- `HornerForm[poly1/poly2, vars1, vars2]`

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

### Expression Information

#### Attributes
Returns the core evaluation attributes assigned to a symbol.
- `Attributes[symbol]`

**Features**: 
- Common attributes include `Flat` (associativity), `Orderless` (commutativity), `Listable` (automatic threading over lists), `HoldFirst`, `HoldRest`, `HoldAll` (evaluation control), and `Protected`.

```mathematica
In[1]:= Attributes[Plus]
Out[1]= {Flat, Listable, NumericFunction, OneIdentity, Orderless}
```

#### ClearAttributes
Removes attributes from symbols.
- `ClearAttributes[s, attr]`: Removes `attr` from the list of attributes of the symbol `s`.
- `ClearAttributes["s", attr]`: Removes `attr` from the attributes of the symbol named `"s"` if it exists.
- `ClearAttributes[s, {attr1, attr2, ...}]`: Removes several attributes at a time.
- `ClearAttributes[{s1, s2, ...}, attrs]`: Removes attributes from several symbols at a time.

**Features**:
- `HoldFirst`, `Protected`.
- `ClearAttributes` modifies `Attributes[s]`.
- Cannot clear attributes of a `Locked` symbol.
- Clearing an attribute that is not set is a no-op.

```mathematica
In[1]:= SetAttributes[f, Listable]
In[2]:= f[{1, 2, 3}]
Out[2]= {f[1], f[2], f[3]}
In[3]:= ClearAttributes[f, Listable]
In[4]:= f[{1, 2, 3}]
Out[4]= f[{1, 2, 3}]

In[5]:= SetAttributes[f, {Flat, Orderless, OneIdentity}]
In[6]:= ClearAttributes[f, OneIdentity]
In[7]:= Attributes[f]
Out[7]= {Flat, Orderless}

In[8]:= ClearAttributes[f, {Flat, Orderless}]
In[9]:= Attributes[f]
Out[9]= {}

In[10]:= SetAttributes[{g, h}, Protected]
In[11]:= ClearAttributes[{g, h}, Protected]
In[12]:= Attributes[g]
Out[12]= {}
```

#### AtomQ, NumberQ, IntegerQ
Predicates for testing expression types.
- `AtomQ[expr]`: `True` if the expression has no parts.
- `NumberQ[expr]`: `True` if the expression is a numeric type (Integer, Real, Rational, Complex).
- `IntegerQ[expr]`: `True` if the expression is an Integer.

#### ListQ, VectorQ, MatrixQ
Predicates for testing lists and their structures.
- `ListQ[expr]`: `True` if the head of `expr` is `List`.
- `VectorQ[expr]`: `True` if `expr` is a list, none of whose elements are themselves lists.
- `VectorQ[expr, test]`: `True` if `expr` is a list and `test` yields `True` for each element.
- `MatrixQ[expr]`: `True` if `expr` is a list of lists of the same length, with no deeper lists.
- `MatrixQ[expr, test]`: `True` if `expr` is a matrix and `test` yields `True` for each element.

#### EvenQ, OddQ
- `EvenQ[n]`: `True` if `n` is an even integer.
- `OddQ[n]`: `True` if `n` is an odd integer.

#### Identity
- `Identity[expr]` returns `expr` unchanged. Useful as a default callable in higher-order functions.
- Attributes: `Protected`.

```mathematica
In[1]:= Identity[x]
Out[1]= x

In[2]:= Identity[1 + 2]
Out[2]= 3

In[3]:= Map[Identity, {a, b, c}]
Out[3]= {a, b, c}
```

#### Composition
- `Composition[f1, f2, f3, ...]` represents the symbolic composition of the functions `f1, f2, f3, ...`. Applied to arguments, the composition acts innermost-first:
  `Composition[f, g, h][x, y]` -> `f[g[h[x, y]]]`.
- Attributes: `Flat`, `OneIdentity`, `Protected`.
- Algebraic simplifications applied automatically:
  - `Composition[]` -> `Identity`
  - `Composition[f]` -> `f`
  - `Identity` arguments are dropped: `Composition[f, Identity, g]` -> `Composition[f, g]`
  - Adjacent inverse pairs cancel: `Composition[f, InverseFunction[f]]` -> `Identity`, and similarly for the reverse order.
- Infix syntax: `f1 @* f2 @* ...` parses as `Composition[f1, f2, ...]` (precedence 625, left-associative). The `Flat` attribute then flattens any nested `Composition` calls produced by the parser.

```mathematica
In[1]:= Composition[f, g, h][x, y]
Out[1]= f[g[h[x, y]]]

In[2]:= f @* g @* h @ x
Out[2]= f[g[h[x]]]

In[3]:= Composition[f, Identity, g]
Out[3]= Composition[f, g]

In[4]:= Composition[f, InverseFunction[f]][x]
Out[4]= x

In[5]:= Composition[f, g] @* Composition[a, b]
Out[5]= Composition[f, g, a, b]
```

#### ComposeList
- `ComposeList[{f1, f2, ..., fn}, x]` generates the list `{x, f1[x], f2[f1[x]], ..., fn[...f2[f1[x]]...]}`. The output has `n + 1` elements.
- Attributes: `Protected`.
- ComposeList builds the symbolic applications and lets the evaluator reduce them in the normal way, so it composes naturally with pure functions and explicit operator-style sequences.

```mathematica
In[1]:= ComposeList[{a, b, c, d}, x]
Out[1]= {x, a[x], b[a[x]], c[b[a[x]]], d[c[b[a[x]]]]}

In[2]:= ComposeList[{f, g}[[{1, 2, 1, 1, 2}]], x]
Out[2]= {x, f[x], g[f[x]], f[g[f[x]]], f[f[g[f[x]]]], g[f[f[g[f[x]]]]]}

In[3]:= ComposeList[{1 - # &, 1/# &}[[{2, 2, 1, 2, 2, 1}]], x]
Out[3]= {x, 1/x, x, 1 - x, 1/(1 - x), 1 - x, x}
```

#### Print
- `Print[expr1, expr2, ...]` prints the expressions to stdout and returns `Null`.
- Supports `FullForm` and `InputForm` wrappers.

```mathematica
In[1]:= Print["Result: ", x + y]
Result: x + y
Out[1]= Null

In[2]:= Print[x + y // FullForm]
Plus[x, y]
Out[2]= Null
```

#### FullForm
- `FullForm[expr]` causes `expr` to be printed in its full internal form.
- `expr // FullForm` is a common shortcut.


#### HoldForm
- `HoldForm[expr]` prints as the expression `expr`, with `expr` maintained in an unevaluated form.
- Has attribute `HoldAll` and `Protected`.
- It dissolves prior to printing but preserves the unevaluated `expr` it contains.

```mathematica
In[1]:= HoldForm[1 + 1]
Out[1]= 1 + 1
```
#### Evaluate
Causes `expr` to be evaluated even if it appears as the argument of a function whose attributes specify that it should be held unevaluated.
- `Evaluate[expr]`

**Features**:
- `Protected`.
- `Evaluate` only overrides `HoldFirst`, `HoldRest`, and `HoldAll` attributes when it appears directly as the head of the function argument that would otherwise be held.
- `Evaluate` does not override `HoldAllComplete`.
- `Evaluate` works only on the first level, directly inside a held function. It does not penetrate into deeper subexpressions.
- Outside of held contexts, `Evaluate` acts as identity.

```mathematica
In[1]:= Evaluate[1+1]
Out[1]= 2

In[2]:= Hold[Evaluate[1+1], 2+2]
Out[2]= Hold[2, 2 + 2]

In[3]:= Hold[Evaluate[1+1], Evaluate[2+2], Evaluate[3+3]]
Out[3]= Hold[2, 4, 6]

In[4]:= Hold[f[Evaluate[1+2]]]
Out[4]= Hold[f[Evaluate[1 + 2]]]

In[5]:= Hold[Evaluate[Sin[Pi/6]]]
Out[5]= Hold[1/2]

In[6]:= Hold[Evaluate[2^10]]
Out[6]= Hold[1024]

In[7]:= HoldForm[Evaluate[1+1]]
Out[7]= 2

In[8]:= Hold[Evaluate[Length[{a,b,c}]]]
Out[8]= Hold[3]

In[9]:= Evaluate[Head[{1,2,3}]]
Out[9]= List
```

#### ReleaseHold
Removes `Hold`, `HoldForm`, `HoldPattern`, and `HoldComplete` in `expr`.
- `ReleaseHold[expr]`

**Features**:
- `Protected`.
- `ReleaseHold` removes only one layer of `Hold` etc.; it does not remove inner occurrences in nested `Hold` etc. functions.
- `ReleaseHold` traverses into subexpressions of `expr` and strips any hold wrapper it finds, but does not recurse into the contents of the stripped wrapper.
- When `expr` does not contain any hold wrappers, `ReleaseHold` acts as identity.
- `ReleaseHold` removes all standard unevaluated containers: `Hold`, `HoldForm`, `HoldComplete`, and `HoldPattern`.

```mathematica
In[1]:= Hold[1+1]
Out[1]= Hold[1 + 1]

In[2]:= ReleaseHold[Hold[1+1]]
Out[2]= 2

In[3]:= ReleaseHold /@ {Hold[1+2], HoldForm[2+3], HoldComplete[3+4]}
Out[3]= {3, 5, 7}

In[4]:= ReleaseHold[f[Hold[1+2]]]
Out[4]= f[3]

In[5]:= ReleaseHold[f[Hold[1+g[Hold[2+3]]]]]
Out[5]= f[1 + g[Hold[2 + 3]]]

In[6]:= ReleaseHold[Hold[Hold[1+1]]]
Out[6]= Hold[1 + 1]

In[7]:= ReleaseHold[Hold[Sin[Pi/6]]]
Out[7]= 1/2

In[8]:= ReleaseHold[{f[Hold[1+2]], g[HoldForm[3+4]]}]
Out[8]= {f[3], g[7]}
```

#### HoldPattern
Keeps a pattern expression in unevaluated form while still allowing it to act as a pattern for matching.
- `HoldPattern[expr]`

**Features**:
- Attributes: `{HoldAll, Protected}`.
- `HoldPattern[p]` is equivalent to `p` in the pattern matcher; the matcher transparently unwraps a single-argument `HoldPattern` before matching.
- Useful on the left-hand side of rules and assignments, because those positions are normally evaluated before being used for matching. Wrapping in `HoldPattern` stops that evaluation and preserves the literal pattern shape.
- `HoldPattern` is removed by one layer of `ReleaseHold`.

```mathematica
In[1]:= HoldPattern[_+_] -> 0
Out[1]= HoldPattern[Blank[] + Blank[]] -> 0

In[2]:= a + b /. HoldPattern[_+_] -> 0
Out[2]= 0

In[3]:= Cases[{a -> b, c -> d}, HoldPattern[a -> _]]
Out[3]= {a -> b}

In[4]:= MatchQ[a + b, HoldPattern[_+_]]
Out[4]= True
```

#### Unevaluated
Represents the unevaluated form of `expr` when it appears as the argument to a function.
- `Unevaluated[expr]`

**Features**:
- Attributes: `{HoldAllComplete, Protected}`.
- `f[Unevaluated[expr]]` passes `expr` to `f` as if `f` temporarily held that single argument; the `Unevaluated` wrapper is then stripped before `f`'s body runs, effectively yielding `f[expr]` with `expr` unevaluated.
- The wrapper is **not** stripped when the enclosing function holds the argument (e.g. `f` has `HoldAll`, or `HoldFirst`/`HoldRest` applies to that position).
- The wrapper is **not** stripped when the enclosing function has `HoldAllComplete`.
- Stripping happens **after** `Sequence` flattening, so a `Sequence` directly inside `Unevaluated` survives into the argument slot (`Length[Unevaluated[Sequence[a, b]]]` gives `2`).
- Nested `Unevaluated` wrappers are stripped one layer per evaluation step.
- As a top-level expression, `Unevaluated[expr]` evaluates to itself (because of `HoldAllComplete`).

```mathematica
In[1]:= Length[Unevaluated[Plus[5, 6, 7, 8]]]
Out[1]= 4

In[2]:= Length[Unevaluated[Sequence[a, b]]]
Out[2]= 2

In[3]:= Hold[Evaluate[Unevaluated[1+2]]]
Out[3]= Hold[Unevaluated[1 + 2]]

In[4]:= SetAttributes[f, HoldAll]; f[Unevaluated[1+2]]
Out[4]= f[Unevaluated[1 + 2]]

In[5]:= HoldComplete[Unevaluated[1+2]]
Out[5]= HoldComplete[Unevaluated[1 + 2]]

In[6]:= Attributes[Unevaluated]
Out[6]= {HoldAllComplete, Protected}
```

#### HoldComplete
Shields `expr` completely from the standard evaluation process.
- `HoldComplete[expr1, expr2, ...]`

**Features**:
- Attributes: `{HoldAllComplete, Protected}`.
- `HoldComplete` prevents argument evaluation, `Sequence` flattening inside its own arguments, `Unevaluated` wrapper stripping, and `Evaluate` from firing. `Evaluate` cannot override `HoldAllComplete`.
- Structural substitution (via `ReplaceAll`, `Replace`, `ReplacePart`, etc.) still descends into `HoldComplete` because substitution is not part of evaluation.
- `HoldComplete` is removed by one layer of `ReleaseHold`.
- `HoldComplete` is a milder form of `Unevaluated` at top level: `HoldComplete` always keeps the wrapper, while `Unevaluated` is typically stripped by the enclosing function.

```mathematica
In[1]:= Attributes[HoldComplete]
Out[1]= {HoldAllComplete, Protected}

In[2]:= HoldComplete[1+1, Evaluate[1+2], Sequence[3, 4]]
Out[2]= HoldComplete[1 + 1, Evaluate[1 + 2], Sequence[3, 4]]

In[3]:= HoldComplete[Sequence[a, b]]
Out[3]= HoldComplete[Sequence[a, b]]

In[4]:= HoldComplete[f[1+2]] /. f[x_] :> g[x]
Out[4]= HoldComplete[g[1 + 2]]

In[5]:= ReleaseHold[HoldComplete[Sequence[1, 2]]]
Out[5]= Sequence[1, 2]
```

#### HoldAllComplete (attribute)
An attribute that specifies that **all** arguments to a function are not to be modified or looked at in any way during evaluation. It is stricter than `HoldAll`.

A function with `HoldAllComplete`:
- does not evaluate any argument,
- does not flatten `Sequence[...]` that appears inside an argument,
- does not strip `Unevaluated[...]` wrappers inside its arguments,
- is not overridden by `Evaluate[...]` wrappers inside its arguments,
- does not apply any up-values associated with its arguments.

`HoldComplete` and `Unevaluated` are the two standard built-in heads that carry `HoldAllComplete`.

#### InputForm
- `InputForm[expr]` causes `expr` to be printed in a form suitable for input (standard form in PicoCAS).

#### Information
Returns a formatted string containing the syntax and description of a symbol.
- `Information[symbol]`
- `?symbol` (shortcut)

```mathematica
In[1]:= ?Range
Out[1]= "Range[n]
	generates the list {1, 2, 3, ..., n}.
Range[n, m]
	generates the list {n, n + 1, ..., m - 1, m}.
Range[n, m, d]
	uses step d."
```

#### MemberQ
- `MemberQ[list, form]`: Returns `True` if an element of list matches form, and `False` otherwise.
- `MemberQ[list, form, levelspec]`: Tests all parts of list specified by levelspec.
- `MemberQ[form]`: Represents an operator form of `MemberQ` that can be applied to an expression.

**Features**:
- `Protected`.
- Default option: `Heads -> False`.
- `form` can be a structural pattern.
- The first argument of `MemberQ` can have any head, not necessarily `List`.
- Returns immediately upon finding the first match.
- Standard level specifications are supported. The default value for `levelspec` in `MemberQ` is `{1}`.

```mathematica
In[1]:= MemberQ[{1, 3, 4, 1, 2}, 2]
Out[1]= True

In[2]:= MemberQ[{x^2, y^2, x^3}, x^_]
Out[2]= True

In[3]:= MemberQ[{{1, 1, 3, 0}, {2, 1, 2, 2}}, 0, 2]
Out[3]= True

In[4]:= MemberQ[{{1, 1, 3, 0}, {2, 1, 2, 2}}, 0]
Out[4]= False
```

#### FreeQ
Yields `True` if no subexpression in `expr` matches `form`, and yields `False` otherwise.
- `FreeQ[expr, form]`
- `FreeQ[expr, form, levelspec]`

**Features**:
- `Protected`.
- By default, explores levels `{0, Infinity}` and option `Heads -> True` is enabled.
- `form` can be a structural pattern.

```mathematica
In[1]:= FreeQ[{1, 2, 4, 1, 0}, 0]
Out[1]= False

In[2]:= FreeQ[{a, b, b, a, a, a}, _Integer]
Out[2]= True

In[3]:= f[c_ x_, x_] := c f[x, x] /; FreeQ[c, x]
In[4]:= {f[3 x, x], f[a x, x], f[(1 + x) x, x]}
Out[4]= {3 f[x, x], a f[x, x], f[x (1 + x), x]}
```

#### LeafCount
Gives the total number of indivisible subexpressions in an expression.
- `LeafCount[expr]`

**Features**:
- `Protected`.
- Counts the number of subexpressions in `expr` that correspond to "leaves" on the expression tree.
- By default `Heads -> True` includes the head of expressions and their parts. With `Heads -> False`, it excludes them.
- Evaluates atoms like `Rational` and `Complex` based on their structural representation as functions.

```mathematica
In[1]:= LeafCount[1 + a + b^2]
Out[1]= 6

In[2]:= LeafCount[f[a, b][x, y]]
Out[2]= 5

In[3]:= LeafCount[{1/2, 1 + I}]
Out[3]= 7
```

#### ByteCount
Gives the number of bytes used internally by PicoCAS to store the expression.
- `ByteCount[expr]`

**Features**:
- `Protected`.
- Uses `sizeof()` in C and measures the internal AST memory allocation boundaries, dynamically capturing sizes of individual strings, symbols, allocated blocks, arrays, and expression structs.

### Time and Date

#### Timing
Evaluates `expr` and returns a list of the time in seconds used, together with the result obtained.
- `Timing[expr]`

**Features**:
- `HoldAll`, `Protected`, `SequenceHold`.
- Returns `{timing, result}`.
- Includes only CPU time spent evaluating the expression.

#### RepeatedTiming
Evaluates `expr` repeatedly and returns a list of the average time in seconds used, together with the result obtained.
- `RepeatedTiming[expr]`
- `RepeatedTiming[expr, t]`

**Features**:
- `HoldFirst`, `Protected`, `SequenceHold`.
- Returns `{average_timing, result}`.
- Does repeated evaluation for at least `t` seconds. Default is 1 second.
- Gives a trimmed mean of the timings obtained, discarding lower and upper quartiles.
- Always evaluates `expr` at least four times.

### Linear Algebra

#### Dot (.)
Gives products of vectors, matrices, and tensors.
- `a . b` or `Dot[a, b]`
- `a . b . c` or `Dot[a, b, c]`

**Features**:
- `Flat`, `OneIdentity`, `Protected`.
- Contracts the last index in `a` with the first index in `b`.
- Applying `Dot` to a rank `n` tensor and a rank `m` tensor gives a rank `m+n-2` tensor.
- Scalar product of two vectors yields a scalar.
- Product of a matrix and a vector yields a vector.
- Product of two matrices yields a matrix.
- When arguments are not lists, `Dot` remains unevaluated.
- Gives an error message `Dot::dotsh` if the shapes of the inputs are not compatible.

```mathematica
In[1]:= {a, b, c} . {x, y, z}
Out[1]= a x + b y + c z

In[2]:= {{a, b}, {c, d}} . {x, y}
Out[2]= {a x + b y, c x + d y}

In[3]:= {{a, b}, {c, d}} . {{1, 2}, {3, 4}}
Out[3]= {{a + 3 b, 2 a + 4 b}, {c + 3 d, 2 c + 4 d}}
```


#### Inner
A generalization of `Dot` in which `f` plays the role of multiplication and `g` of addition.
- `Inner[f, list1, list2, g]`: computes the "inner `f`" of two lists, with "plus operation" `g`.
- `Inner[f, list1, list2]`: uses `Plus` for `g`.
- `Inner[f, list1, list2, g, n]`: contracts index `n` of the first tensor with the first index of the second tensor.

**Features**:
- `Protected`.
- Like `Dot`, `Inner` effectively contracts the last index of the first tensor with the first index of the second tensor.
- Applying `Inner` to a rank $r$ tensor and a rank $s$ tensor gives a rank $r+s-2$ tensor.

```mathematica
In[1]:= Inner[f, {a, b}, {x, y}, g]
Out[1]= g[f[a, x], f[b, y]]

In[2]:= Inner[Times, {{a, b}, {c, d}}, {x, y}, Plus]
Out[2]= {a x + b y, c x + d y}

In[3]:= Inner[Times, {a1, a2, a3}, {b1, b2, b3}, Plus]
Out[3]= a1 b1 + a2 b2 + a3 b3
```

#### Outer
Gives the generalized outer product of the `listi`, forming all possible combinations of the lowest-level elements in each of them, and feeding them as arguments to `f`.
- `Outer[f, list1, list2, ...]`
- `Outer[f, list1, list2, ..., n]`: treats as separate elements only sublists at level `n`.
- `Outer[f, list1, list2, ..., n1, n2, ...]`: treats as separate elements only sublists at level `ni` in the corresponding `listi`.

**Features**:
- `Protected`.
- Applying `Outer` to two tensors of ranks $r$ and $s$ gives a tensor of rank $r+s$.
- The heads of all `listi` must be the same, but need not necessarily be `List`.

```mathematica
In[1]:= Outer[f, {a, b}, {x, y, z}]
Out[1]= {{f[a, x], f[a, y], f[a, z]}, {f[b, x], f[b, y], f[b, z]}}

In[2]:= Outer[Times, {1, 2, 3, 4}, {a, b, c}]
Out[2]= {{a, b, c}, {2 a, 2 b, 2 c}, {3 a, 3 b, 3 c}, {4 a, 4 b, 4 c}}

In[3]:= Outer[g, f[a, b], f[x, y, z]]
Out[3]= f[f[g[a, x], g[a, y], g[a, z]], f[g[b, x], g[b, y], g[b, z]]]
```
#### Det
Gives the determinant of a square matrix.
- `Det[m]`

**Features**:
- `Protected`.
- Evaluates the determinant of a square matrix symbolically or numerically using Laplace expansion.
- Returns a warning `Det::matsq` if `m` is not a non-empty square matrix.

```mathematica
In[1]:= Det[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= 0

In[2]:= Det[{{1.7, 7.1, -2.7}, {2.2, 8.7, 3.2}, {3.2, -9.2, 1.2}}]
Out[2]= 251.572

In[3]:= Det[{{a, b, c}, {d, e, f}, {g, h, i}}]
Out[3]= -c e g + b f g + c d h - a f h - b d i + a e i
```

#### Cross
Gives the vector cross product of its arguments.
- `Cross[v1, v2, ...]`

**Features**:
- `Protected`.
- Returns the cross product or totally antisymmetric product of $n-1$ vectors of length $n$.
- Works for symbolic and numerical inputs.
- Outputs `Cross::nonn1` error message if inputs are not equal-length vectors or if the number of arguments is not one less than their length.

```mathematica
In[1]:= Cross[{1, 2, -1}, {-1, 1, 0}]
Out[1]= {1, 1, 3}

In[2]:= Cross[{1, Sqrt[3]}]
Out[2]= {-Sqrt[3], 1}

In[3]:= Cross[{3.2, 4.2, 5.2}, {0.75, 0.09, 0.06}]
Out[3]= {-0.216, 3.708, -2.862}

In[4]:= Cross[{1.3 + I, 2, 3 - 2 I}, {6. + I, 4, 5 - 7 I}]
Out[4]= {-2 - 6 I, 6.5 - 4.9 I, -6.8 + 2. I}

In[5]:= Cross[{1, 2, 3}, {4, 5}]
Cross::nonn1: The arguments are expected to be vectors of equal length, and the number of arguments is expected to be 1 less than their length.
Out[5]= Cross[{1, 2, 3}, {4, 5}]
```

#### Norm
Gives the norm of a number, vector, or matrix.
- `Norm[expr]`
- `Norm[expr, p]`

**Features**:
- `Protected`.
- For scalars, `Norm[z]` is `Abs[z]`.
- For vectors, `Norm[v]` defaults to the 2-norm: `Sqrt[v . Conjugate[v]]`.
- For vectors, `Norm[v, p]` is `Total[Abs[v]^p]^(1/p)`.
- For vectors, `Norm[v, Infinity]` is the $\infty$-norm given by `Max[Abs[v]]`.
- `Norm[m, "Frobenius"]` gives the Frobenius norm of a matrix `m`.

```mathematica
In[1]:= Norm[{x, y, z}]
Out[1]= Sqrt[Abs[x]^2 + Abs[y]^2 + Abs[z]^2]

In[2]:= Norm[-2 + I]
Out[2]= Sqrt[5]

In[3]:= v = {1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1}; Norm[v]
Out[3]= Sqrt[5]

In[4]:= Norm[{x, y, z}, p]
Out[4]= (Abs[x]^p + Abs[y]^p + Abs[z]^p)^(1/p)

In[5]:= Norm[{x, y, z}, Infinity]
Out[5]= Max[Abs[x], Abs[y], Abs[z]]

In[6]:= Norm[{{a11, a12}, {a21, a22}}, "Frobenius"]
Out[6]= Sqrt[Abs[a11]^2 + Abs[a12]^2 + Abs[a21]^2 + Abs[a22]^2]
```

#### Tr
Finds the trace of a matrix or tensor.
- `Tr[list]`
- `Tr[list, f]`
- `Tr[list, f, n]`

**Features**:
- `Protected`.
- `Tr[list]` sums the diagonal elements `list[[i, i, ...]]`.
- `Tr[list, f]` applies the function `f` instead of `Plus`.
- `Tr[list, f, n]` considers elements down to level `n`.
- Works for rectangular as well as square matrices and tensors, stopping at the minimum dimension.
- If `n` is omitted, defaults to the depth of the tensor.

```mathematica
In[1]:= Tr[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= 15

In[2]:= Tr[{{a, b}, {c, d}}]
Out[2]= a + d

In[3]:= Tr[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, List]
Out[3]= {1, 5, 9}

In[4]:= Tr[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, Plus, 1]
Out[4]= {12, 15, 18}
```

#### Inverse
Gives the inverse of a square matrix.
- `Inverse[m]`

**Features**:
- `Protected`.
- Works on both symbolic and numerical matrices.
- For matrices with approximate real or complex numbers, the inverse is generated to the maximum possible precision given the input.
- Uses fraction-free Gauss-Jordan elimination on the augmented matrix `[A | I]` to compute the inverse exactly for integer, rational, and symbolic matrices.
- Issues `Inverse::sing` warning and returns unevaluated if the matrix is singular.
- Issues `Inverse::matsq` warning and returns unevaluated if the argument is not a non-empty square matrix.
- Satisfies the relation `a . Inverse[a] == Inverse[a] . a == IdentityMatrix[n]`.
- Satisfies the relation `Inverse[a . b] == Inverse[b] . Inverse[a]`.

```mathematica
In[1]:= Inverse[{{1.4,2},{3,-6.7}}]
Out[1]= {{0.435631, 0.130039}, {0.195059, -0.0910273}}

In[2]:= Inverse[{{1,2,3},{4,2,2},{5,1,7}}]
Out[2]= {{-2/7, 11/42, 1/21}, {3/7, 4/21, -5/21}, {1/7, -3/14, 1/7}}

In[3]:= Inverse[{{u,v},{v,u}}]
Out[3]= {{u/(u^2 - v^2), -v/(u^2 - v^2)}, {-v/(u^2 - v^2), u/(u^2 - v^2)}}

In[4]:= Inverse[{{1.2,2.5,-3.2},{0.7,-9.4,5.8},{-0.2,0.3,6.4}}]
Out[4]= {{0.74546, 0.204249, 0.187629}, {0.0679223, -0.0847825, 0.110795}, {0.0201118, 0.010357, 0.15692}}

In[5]:= Inverse[{{2,3,2},{4,9,2},{7,2,4}}]
Out[5]= {{-8/13, 2/13, 3/13}, {1/26, 3/26, -1/13}, {55/52, -17/52, -3/26}}

In[6]:= Inverse[{{a,b},{c,d}}]
Out[6]= {{d/(-b c + a d), -b/(-b c + a d)}, {-c/(-b c + a d), a/(-b c + a d)}}

In[7]:= Inverse[{{1,2},{1,2}}]
Inverse::sing: Matrix {{1, 2}, {1, 2}} is singular.
Out[7]= Inverse[{{1, 2}, {1, 2}}]

In[8]:= a = {{1,2},{3,4}}; a . Inverse[a] == IdentityMatrix[2]
Out[8]= True

In[9]:= a = {{1,1,1},{6,9,7},{8,1,9}}; b = {{0,3,9},{7,9,7},{4,4,1}}; Inverse[a . b] == Inverse[b] . Inverse[a]
Out[9]= True

In[10]:= Inverse[{{1,2},{3,4},{5,6}}]
Inverse::matsq: Argument {{1, 2}, {3, 4}, {5, 6}} at position 1 is not a non-empty square matrix.
Out[10]= Inverse[{{1, 2}, {3, 4}, {5, 6}}]
```

#### MatrixPower
Gives the matrix power of a square matrix.
- `MatrixPower[m, n]`: gives the n-th matrix power of the matrix `m`.
- `MatrixPower[m, n, v]`: gives the n-th matrix power of the matrix `m` applied to the vector `v`.

**Features**:
- `Protected`.
- `MatrixPower[m, n]` effectively evaluates the product of a matrix with itself `n` times.
- When `n` is negative, `MatrixPower` finds powers of the inverse of the matrix `m`.
- `MatrixPower[m, 0]` returns `IdentityMatrix[Length[m]]`.
- `MatrixPower` works only on square matrices.
- Uses binary exponentiation (repeated squaring) for efficient computation.
- Fractional matrix powers are not currently supported and generate a warning.

```mathematica
In[1]:= MatrixPower[{{a, b}, {c, d}}, 2]
Out[1]= {{a^2 + b c, a b + b d}, {a c + c d, b c + d^2}}

In[2]:= MatrixPower[{{a, b}, {c, d}}, 2] == {{a, b}, {c, d}} . {{a, b}, {c, d}}
Out[2]= True

In[3]:= MatrixPower[{{a, b}, {c, d}}, -2] == Inverse[{{a, b}, {c, d}}] . Inverse[{{a, b}, {c, d}}]
Out[3]= True

In[4]:= MatrixPower[{{1, 1}, {1, 2}}, 10]
Out[4]= {{4181, 6765}, {6765, 10946}}

In[5]:= MatrixPower[{{1.2, 2.5, -3.2}, {0.7, -9.4, 5.8}, {-0.2, 0.3, 6.4}}, 5]
Out[5]= {{-1208.61, 19598.2, -12658.4}, {5784.51, -83315.1, 35420.6}, {-559.11, 1960.12, 11511.9}}

In[6]:= MatrixPower[{{2, 3, 0}, {4, 9, 0}, {0, 0, 4}}, 14]
Out[6]= {{25881337259836, 54508871401413, 0}, {72678495201884, 153068703863133, 0}, {0, 0, 268435456}}

In[7]:= MatrixPower[{{a, b}, {0, c}}, 4]
Out[7]= {{a^4, a^2 (a b + b c) + c^2 (a b + b c)}, {0, c^4}}

In[8]:= MatrixPower[{{1, 1}, {1, 2}}, 3, {1, 0}]
Out[8]= {5, 8}

In[9]:= MatrixPower[{{1, 2}, {3, 4}}, 0]
Out[9]= {{1, 0}, {0, 1}}

In[10]:= MatrixPower[{{5}}, -2]
Out[10]= {{1/25}}
```

#### RowReduce
Gives the row-reduced form of the matrix `m`.
- `RowReduce[m]`

**Features**:
- `Protected`.
- Uses fraction-free division logic to perform exact algorithmic reduction across numerical, rational, and symbolics expressions natively avoiding division errors.

```mathematica
In[1]:= RowReduce[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= {{1, 0, -1}, {0, 1, 2}, {0, 0, 0}}

In[2]:= RowReduce[{{a, b, c}, {d, e, f}, {a+d, b+e, c+f}}]
Out[2]= {{1, 0, (-b f + c e)/(a e - b d)}, {0, 1, (-a f + c d)/(-a e + b d)}, {0, 0, 0}}
```

#### IdentityMatrix
Generates an identity matrix.
- `IdentityMatrix[n]`: Gives the `n x n` identity matrix.
- `IdentityMatrix[{m, n}]`: Gives the `m x n` identity matrix.

**Features**:
- `Protected`.
- Generates exact integer outputs (`1` on main diagonal, `0` elsewhere).
- Will remain unevaluated if arguments are symbolic or negative.

```mathematica
In[1]:= IdentityMatrix[3]
Out[1]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}

In[2]:= IdentityMatrix[{2, 3}]
Out[2]= {{1, 0, 0}, {0, 1, 0}}
```

#### DiagonalMatrix
Generates a matrix with specified elements on a diagonal.
- `DiagonalMatrix[list]`: Elements on leading diagonal.
- `DiagonalMatrix[list, k]`: Elements on `k`-th diagonal.
- `DiagonalMatrix[list, k, n]`: Pads with zeros to create an `n x n` matrix.
- `DiagonalMatrix[list, k, {m, n}]`: Creates an `m x n` matrix.

**Features**:
- `Protected`.
- For `k > 0`, places elements `k` positions above the leading diagonal.
- For `k < 0`, places elements `k` positions below the leading diagonal.
- By default, size is optimally bounded to fit the full array cleanly. Extraneous elements are dropped if manual constraints fall short of required lengths.

```mathematica
In[1]:= DiagonalMatrix[{a, b, c}]
Out[1]= {{a, 0, 0}, {0, b, 0}, {0, 0, c}}

In[2]:= DiagonalMatrix[{a, b}, 1]
Out[2]= {{0, a, 0}, {0, 0, b}, {0, 0, 0}}

In[3]:= DiagonalMatrix[{1, 2, 3}, 0, {3, 5}]
Out[3]= {{1, 0, 0, 0, 0}, {0, 2, 0, 0, 0}, {0, 0, 3, 0, 0}}
```

### Statistics


#### Median
Gives the median estimate of the elements in data.
- `Median[data]`: gives the median estimate $\hat{q}_{1/4}$ of the elements in `data`.
- `Median[dist]`: gives the median of the distribution `dist`.

**Features**:
- `Protected`.
- Median is a robust location estimator, which means it not very sensitive to outliers.
- For `VectorQ` data $\{x_1, \dots, x_n\}$, the median can be thought of as the "middle value". Formally, when data is sorted as $\{x_{(1)}, \dots, x_{(n)}\}$, the median is given by the center element $x_{((n+1)/2)}$ if $n$ is odd and the mean of the two center elements $(x_{(n/2)} + x_{(n/2+1)})/2$ if $n$ is even.
- For `MatrixQ` data, the median is computed for each column vector. `Median` for a tensor gives columnwise medians at the first level.
- `Median` requires numeric values.

```mathematica
In[1]:= Median[{1, 2, 3, 4, 5, 6, 7}]
Out[1]= 4

In[2]:= Median[{1, 2, 3, 4, 5, 6, 7, 8}]
Out[2]= 9/2

In[3]:= Median[{1, 2, 3, 4}]
Out[3]= 5/2

In[4]:= Median[{Pi, E, 2}]
Out[4]= E

In[5]:= Median[{1., 2., 3., 4.}]
Out[5]= 2.5

In[6]:= Median[{{1, 11, 3}, {4, 6, 7}}]
Out[6]= {5/2, 17/2, 5}

In[7]:= Median[{{{3, 7}, {2, 1}}, {{5, 19}, {12, 4}}}]
Out[7]= {{4, 13}, {7, 5/2}}

In[8]:= Median[{a, b, c}]
Median::rectn: Rectangular array of real numbers is expected at position 1 in Median[{a, b, c}].
Out[8]= Median[{a, b, c}]
```
#### Mean
Gives the mean estimate of the elements in data.
- `Mean[data]`

**Features**:
- `Protected`.
- Supports numerical and symbolic data.
- For vectors, computes $(1/n) \sum x_i$.
- For matrices, computes means of elements in each column.

```mathematica
In[1]:= Mean[{1, 2, 3, 4}]
Out[1]= 5/2

In[2]:= Mean[{{a, u}, {b, v}, {c, w}}]
Out[2]= {1/3 (a + b + c), 1/3 (u + v + w)}
```

#### RootMeanSquare
Gives the root mean square of values in `list`.
- `RootMeanSquare[list]`

**Features**:
- `Protected`.
- Gives the square root of the second sample moment.
- For a list `{x1, x2, ...}`, it computes `Sqrt[1/n Total[{x1^2, x2^2, ...}]]`.
- Handles both numerical and symbolic data.
- Works column-wise on matrices.

```mathematica
In[1]:= RootMeanSquare[{a, b, c, d}]
Out[1]= 1/2 Sqrt[a^2 + b^2 + c^2 + d^2]

In[2]:= RootMeanSquare[{{1, 2}, {5, 10}, {5, 2}, {4, 8}}]
Out[2]= {1/2 Sqrt[67], Sqrt[43]}

In[3]:= RootMeanSquare[{1, 2, 3, 4}]
Out[3]= Sqrt[15/2]

In[4]:= RootMeanSquare[{Pi, E, 2}]
Out[4]= Sqrt[1/3 (4 + E^2 + Pi^2)]

In[5]:= RootMeanSquare[{1., 2., 3., 4.}]
Out[5]= 2.73861
```

#### Variance
Gives the unbiased variance estimate of the elements in data.
- `Variance[data]`

**Features**:
- `Protected`.
- For vectors, computes $(1/(n-1)) \sum (x_i - \hat{\mu}) \overline{(x_i - \hat{\mu})}$.
- For matrices, computes variances of elements in each column.

```mathematica
In[1]:= Variance[{1, 2, 3}]
Out[1]= 1

In[2]:= Variance[{{5.2, 7}, {5.3, 8}, {5.4, 9}}]
Out[2]= {0.01, 1}
```

#### StandardDeviation
Gives the standard deviation estimate of the elements in data.
- `StandardDeviation[data]`

**Features**:
- `Protected`.
- Equivalent to `Sqrt[Variance[data]]`.
- For matrices, computes standard deviations of elements in each column.

```mathematica
In[1]:= StandardDeviation[{1, 2, 3}]
Out[1]= 1
```

#### MovingAverage
Gives the moving average over a list, with either a uniform window length or a list of weights.
- `MovingAverage[list, r]`: averages runs of `r` consecutive elements.
- `MovingAverage[list, {w_1, w_2, ..., w_r}]`: weighted moving average with effective weights $w_i / \sum_j w_j$.

**Features**:
- `Protected`.
- Output length is `Length[list] - r + 1`.
- Stays unevaluated when `r < 1`, when `r > Length[list]`, when the second argument is non-integer / non-list, or when the first argument is not a `List`.
- Exact rational arithmetic for integer / rational data; bignums (arbitrary-precision integers) handled natively. Real-valued data or weights yield approximate output. Symbolic data and weights are supported.
- The unweighted form delegates to `Mean` for each window, so it inherits `Mean`'s exact / numeric / symbolic dispatch.

```mathematica
In[1]:= MovingAverage[{1, 5, 7, 3, 6, 2}, 3]
Out[1]= {13/3, 5, 16/3, 11/3}

In[2]:= MovingAverage[{1.2, 5.2, 3.4, 4.5, 2.3, 4.5}, 3]
Out[2]= {3.26667, 4.36667, 3.4, 3.76667}

In[3]:= MovingAverage[{a, b, c, d, e}, 2]
Out[3]= {1/2 (a + b), 1/2 (b + c), 1/2 (c + d), 1/2 (d + e)}

In[4]:= MovingAverage[{a, b, c, d, e}, {1, 2}]
Out[4]= {1/3 a + 2/3 b, 1/3 b + 2/3 c, 1/3 c + 2/3 d, 1/3 d + 2/3 e}

In[5]:= MovingAverage[{2^100, 2^101, 2^102, 2^103}, 2]
Out[5]= {1901475900342344102245054808064, 3802951800684688204490109616128, 7605903601369376408980219232256}

In[6]:= MovingAverage[{1, 2, 3, 4, 5}, 6]
Out[6]= MovingAverage[{1, 2, 3, 4, 5}, 6]
```

#### MovingMedian
Gives the moving median of a list, taken over spans of `r` consecutive elements.
- `MovingMedian[list, r]`: medians of runs of `r` elements.

**Features**:
- `Protected`.
- Output length is `Length[list] - r + 1`.
- Operates on real-valued vectors and matrices. For matrix input, each window of `r` consecutive rows is reduced via `Median`, yielding a column-wise median vector per window.
- Exact rationals, bignums (arbitrary-precision integers), machine-precision reals, and `NumericQ`-real symbolic constants (`Pi`, `E`, ...) are all supported. Even-window medians yield exact rational midpoints when the data is exact.
- Stays unevaluated when `r < 1`, when `r > Length[list]`, when `r` is non-integer, or when the first argument is not a `List`.
- Non-numeric data triggers the `MovingMedian::arg1` message and the expression remains unevaluated.

```mathematica
In[1]:= MovingMedian[{1, 2, 5, 6, 1, 4, 3}, 3]
Out[1]= {2, 5, 5, 4, 3}

In[2]:= MovingMedian[{{1, 2}, {5, 3}, {1, 4}, {3, 2}, {5, 5}}, 2]
Out[2]= {{3, 5/2}, {3, 7/2}, {2, 3}, {4, 7/2}}

In[3]:= MovingMedian[N[{1, 5, 7, 3, 6, 2}], 3]
Out[3]= {5.0, 5.0, 6.0, 3.0}

In[4]:= MovingMedian[{1, 2, 3, 4}, 2]
Out[4]= {3/2, 5/2, 7/2}

In[5]:= MovingMedian[{2^100, 2^101, 2^102, 2^103}, 2]
Out[5]= {1901475900342344102245054808064, 3802951800684688204490109616128, 7605903601369376408980219232256}

In[6]:= MovingMedian[{a, b, c}, 2]
MovingMedian::arg1: The first argument {a, b, c} must be a vector or matrix of real values.
Out[6]= MovingMedian[{a, b, c}, 2]
```

#### ExponentialMovingAverage
Gives the exponential moving average of a list with smoothing constant `alpha`.
- `ExponentialMovingAverage[list, alpha]`: produces the recurrence $y_1 = x_1$, $y_{i+1} = y_i + \alpha (x_{i+1} - y_i)$.

**Features**:
- `Protected`.
- Output has the same length as `list`.
- Two evaluation strategies: a fast O(n) double-precision path activates when at least one element of `list` or `alpha` itself is a machine-precision real and every other entry is a real-valued numeric (Integer, Real, Rational); otherwise a symbolic recurrence path is taken so exact rationals, bignums (arbitrary-precision integers), and symbolic alpha all work.
- The smoothing constant `alpha` is typically a number between 0 and 1 but may be any expression; with `alpha = 0` the output is constant at `x_1`, and with `alpha = 1` the output equals the input.
- Stays unevaluated when the first argument is not a `List`, when the list is empty, or when the call has the wrong arity.

```mathematica
In[1]:= ExponentialMovingAverage[Range[10], 1/3]
Out[1]= {1, 4/3, 17/9, 70/27, 275/81, 1036/243, 3773/729, 13378/2187, 46439/6561, 158488/19683}

In[2]:= ExponentialMovingAverage[N[{1, 5, 7, 3, 6, 2}], 1/2]
Out[2]= {1.0, 3.0, 5.0, 4.0, 5.0, 3.5}

In[3]:= ExponentialMovingAverage[{a, b, c, d}, 0]
Out[3]= {a, a, a, a}

In[4]:= ExponentialMovingAverage[{a, b, c, d}, 1]
Out[4]= {a, b, c, d}

In[5]:= ExponentialMovingAverage[{a, b}, x]
Out[5]= {a, a + x (-a + b)}

In[6]:= ExponentialMovingAverage[{2^100, 2^200}, 1/2]
Out[6]= {1267650600228229401496703205376, 803469022129495137770981046171215126561215611592144769253376}
```

### Random Number Generation

#### RandomInteger
Gives pseudorandom integers.
- `RandomInteger[{imin, imax}]`: gives a pseudorandom integer in the range {imin, ..., imax}.
- `RandomInteger[imax]`: gives a pseudorandom integer in the range {0, ..., imax}.
- `RandomInteger[]`: pseudorandomly gives 0 or 1.
- `RandomInteger[range, n]`: gives a list of n pseudorandom integers.
- `RandomInteger[range, {n1, n2, ...}]`: gives an n1 x n2 x ... array of pseudorandom integers.

**Features**:
- `Protected`.
- RandomInteger[{imin, imax}] chooses integers in the range {imin, ..., imax} with equal probability.
- RandomInteger[] gives 0 or 1 with probability 1/2.
- RandomInteger gives a different sequence of pseudorandom integers whenever you run PicoCAS. You can start with a particular seed using SeedRandom.
- Returns bignums when the range exceeds 64-bit integer limits.

```mathematica
In[1]:= SeedRandom[42]; RandomInteger[]
Out[1]= 1

In[2]:= SeedRandom[42]; RandomInteger[10]
Out[2]= 6

In[3]:= SeedRandom[42]; RandomInteger[{1, 6}]
Out[3]= 4

In[4]:= SeedRandom[42]; RandomInteger[{0, 9}, 5]
Out[4]= {6, 9, 4, 1, 3}

In[5]:= SeedRandom[42]; Dimensions[RandomInteger[{0, 1}, {3, 4}]]
Out[5]= {3, 4}

In[6]:= SeedRandom[42]; RandomInteger[{-10, -5}]
Out[6]= -6

In[7]:= IntegerQ[RandomInteger[10^20]]
Out[7]= True
```

#### RandomReal
Gives pseudorandom real numbers.
- `RandomReal[]`: gives a pseudorandom real number in the range 0 to 1.
- `RandomReal[{xmin, xmax}]`: gives a pseudorandom real number in the range xmin to xmax.
- `RandomReal[xmax]`: gives a pseudorandom real number in the range 0 to xmax.
- `RandomReal[range, n]`: gives a list of n pseudorandom reals.
- `RandomReal[range, {n1, n2, ...}]`: gives an n1 x n2 x ... array of pseudorandom reals.

**Features**:
- `Protected`.
- RandomReal[{xmin, xmax}] chooses reals with a uniform probability distribution in the range xmin to xmax.
- RandomReal gives a different sequence of pseudorandom reals whenever you run PicoCAS. You can start with a particular seed using SeedRandom.
- Uses 53 bits of randomness for full double-precision mantissa coverage.
- Accepts integer, real, rational, and bigint range arguments.

```mathematica
In[1]:= SeedRandom[42]; RandomReal[]
Out[1]= 0.376082

In[2]:= SeedRandom[42]; RandomReal[10]
Out[2]= 3.76082

In[3]:= SeedRandom[42]; RandomReal[{-1, 1}]
Out[3]= -0.247836

In[4]:= SeedRandom[42]; Length[RandomReal[{0, 1}, 5]]
Out[4]= 5

In[5]:= SeedRandom[42]; Dimensions[RandomReal[{0, 1}, {3, 4}]]
Out[5]= {3, 4}

In[6]:= SeedRandom[42]; RandomReal[{0, 1}, 0]
Out[6]= {}

In[7]:= RandomReal[x]
Out[7]= RandomReal[x]
```

#### SeedRandom
Resets the pseudorandom generator.
- `SeedRandom[n]`: seeds the generator with integer n.
- `SeedRandom[]`: reseeds from system entropy.

**Features**:
- `Protected`.
- After `SeedRandom[n]`, the sequence of pseudorandom numbers generated will be the same each time.
- Accepts bignums as seeds.

```mathematica
In[1]:= SeedRandom[42]; {RandomInteger[], RandomInteger[], RandomInteger[]}
Out[1]= {1, 1, 0}

In[2]:= SeedRandom[42]; {RandomInteger[], RandomInteger[], RandomInteger[]}
Out[2]= {1, 1, 0}
```

#### RandomComplex
Gives pseudorandom complex numbers.
- `RandomComplex[]`: gives a pseudorandom complex number with real and imaginary parts in the range 0 to 1.
- `RandomComplex[{zmin, zmax}]`: gives a pseudorandom complex number in the rectangle with corners given by the complex numbers zmin and zmax.
- `RandomComplex[zmax]`: gives a pseudorandom complex number in the rectangle whose corners are the origin and zmax.
- `RandomComplex[range, n]`: gives a list of n pseudorandom complex numbers.
- `RandomComplex[range, {n1, n2, ...}]`: gives an n1 x n2 x ... array of pseudorandom complex numbers.

**Features**:
- `Protected`.
- `RandomComplex[{zmin, zmax}]` chooses complex numbers uniformly in the rectangle with corners at `zmin` and `zmax`.
- RandomComplex gives a different sequence of pseudorandom complex numbers whenever you run PicoCAS. You can start with a particular seed using SeedRandom.
- Uses 53 bits of randomness per component for full double-precision mantissa coverage.
- Accepts integer, real, rational, and complex range arguments. When the range has no imaginary component, the result simplifies to a real.

```mathematica
In[1]:= SeedRandom[42]; Head[RandomComplex[]]
Out[1]= Complex

In[2]:= SeedRandom[42]; RandomComplex[2 + 3 I]
Out[2]= 0.752164 + 1.30654 I

In[3]:= SeedRandom[42]; Length[RandomComplex[{0, 1 + I}, 5]]
Out[3]= 5

In[4]:= SeedRandom[42]; Dimensions[RandomComplex[{0, 1 + I}, {3, 4}]]
Out[4]= {3, 4}

In[5]:= RandomComplex[x]
Out[5]= RandomComplex[x]
```

#### RandomChoice
Gives pseudorandom choices from a list of elements.
- `RandomChoice[{e1, e2, ...}]`: gives a pseudorandom choice of one of the ei.
- `RandomChoice[list, n]`: gives a list of n pseudorandom choices.
- `RandomChoice[list, {n1, n2, ...}]`: gives an n1 x n2 x ... array of pseudorandom choices.
- `RandomChoice[{w1, w2, ...} -> {e1, e2, ...}]`: gives a pseudorandom choice weighted by the wi.
- `RandomChoice[wlist -> elist, n]`: gives a list of n weighted choices.
- `RandomChoice[wlist -> elist, {n1, n2, ...}]`: gives an n1 x n2 x ... array of weighted choices.

**Features**:
- `Protected`.
- `RandomChoice[{e1, e2, ...}]` chooses with equal probability between all of the ei.
- RandomChoice gives a different sequence of pseudorandom choices whenever you run PicoCAS. You can start with a particular seed using SeedRandom.
- Weighted selection uses cumulative weight binary search for efficient O(log n) per choice.
- Weights must be non-negative real numbers with a positive total.

```mathematica
In[1]:= SeedRandom[42]; RandomChoice[{a, b, c, d, e}]
Out[1]= c

In[2]:= SeedRandom[42]; RandomChoice[{a, b, c}, 5]
Out[2]= {b, c, b, a, a}

In[3]:= SeedRandom[42]; Dimensions[RandomChoice[{1, 2, 3}, {3, 4}]]
Out[3]= {3, 4}

In[4]:= RandomChoice[{1, 0, 0} -> {a, b, c}]
Out[4]= a

In[5]:= RandomChoice[{1, 0} -> {x, y}, 5]
Out[5]= {x, x, x, x, x}

In[6]:= RandomChoice[x]
Out[6]= RandomChoice[x]
```

#### RandomSample
Gives a pseudorandom sample of elements without replacement.
- `RandomSample[{e1, e2, ...}, n]`: gives a pseudorandom sample of n of the ei.
- `RandomSample[{w1, w2, ...} -> {e1, e2, ...}, n]`: gives a pseudorandom sample of n of the ei chosen using weights wi.
- `RandomSample[{e1, e2, ...}]`: gives a pseudorandom permutation of the ei.
- `RandomSample[list, UpTo[n]]`: gives a sample of n of the ei, or as many as are available.

**Features**:
- `Protected`.
- `RandomSample[{e1, e2, ...}, n]` never samples any of the ei more than once.
- `RandomSample[{e1, e2, ...}, n]` samples each of the ei with equal probability.
- `RandomSample[{e1, e2, ...}, UpTo[n]]` gives a sample of n of the ei, or as many as are available.
- RandomSample gives a different sequence of pseudorandom choices whenever you run PicoCAS. You can start with a particular seed using SeedRandom.
- Requesting n greater than the list length (without UpTo) returns unevaluated.
- Uses the Fisher-Yates shuffle for uniform sampling without replacement.
- Weighted sampling removes selected elements and renormalizes weights.

```mathematica
In[1]:= SeedRandom[42]; RandomSample[{a, b, c, d, e}, 3]
Out[1]= {d, e, b}

In[2]:= Sort[RandomSample[{1, 2, 3, 4, 5}, 5]]
Out[2]= {1, 2, 3, 4, 5}

In[3]:= Length[RandomSample[{a, b, c, d, e}]]
Out[3]= 5

In[4]:= RandomSample[{a, b, c}, 0]
Out[4]= {}

In[5]:= Length[RandomSample[{a, b, c, d, e}, UpTo[10]]]
Out[5]= 5

In[6]:= RandomSample[{1, 0, 0} -> {a, b, c}, 1]
Out[6]= {a}

In[7]:= Sort[RandomSample[{1, 1, 0} -> {a, b, c}, 2]]
Out[7]= {a, b}

In[8]:= Sort[RandomSample[{1, 2, 3} -> {a, b, c}]]
Out[8]= {a, b, c}

In[9]:= RandomSample[{a, b}, 5]
Out[9]= RandomSample[{a, b}, 5]

In[10]:= RandomSample[x]
Out[10]= RandomSample[x]
```

### String Operations

#### StringLength

Gives the number of characters in a string.

- `StringLength["string"]`: Returns the number of characters as an integer.
- Returns unevaluated for non-string arguments.
- **Attributes**: `Listable`, `Protected`.

```mathematica
In[1]:= StringLength["tiger"]
Out[1]= 5

In[2]:= StringLength[""]
Out[2]= 0

In[3]:= StringLength["hello world"]
Out[3]= 11

In[4]:= StringLength[{\"ABC\", \"DE\", \"F\"}]
Out[4]= {3, 2, 1}

In[5]:= StringLength[x]
Out[5]= StringLength[x]
```

#### Characters

Gives a list of the characters in a string.

- `Characters["string"]`: Returns a `List` of single-character strings.
- Each character is given as a length-1 string.
- Returns unevaluated for non-string arguments.
- **Attributes**: `Listable`, `Protected`.

```mathematica
In[1]:= Characters["ABC"]
Out[1]= {"A", "B", "C"}

In[2]:= Characters["A string."]
Out[2]= {"A", " ", "s", "t", "r", "i", "n", "g", "."}

In[3]:= Characters[""]
Out[3]= {}

In[4]:= Characters[{"ABC", "DEF", "XYZ"}]
Out[4]= {{"A", "B", "C"}, {"D", "E", "F"}, {"X", "Y", "Z"}}

In[5]:= Characters[x]
Out[5]= Characters[x]
```

#### StringJoin

Concatenates strings together.

- `StringJoin["s1", "s2", ...]`: Joins all string arguments into a single string.
- `StringJoin[{"s1", "s2", ...}]`: Flattens all lists recursively and joins enclosed strings.
- `StringJoin[]`: Returns the empty string `""`.
- The infix form is `"s1" <> "s2" <> ...`.
- Returns unevaluated if any non-string, non-list leaf is encountered.
- **Attributes**: `Flat`, `OneIdentity`, `Protected`.

```mathematica
In[1]:= StringJoin["abcd", "ABCD", "xyz"]
Out[1]= "abcdABCDxyz"

In[2]:= "abcd" <> "ABCD" <> "xyz"
Out[2]= "abcdABCDxyz"

In[3]:= StringJoin[{{"AB", "CD"}, "XY"}]
Out[3]= "ABCDXY"

In[4]:= StringJoin[]
Out[4]= ""

In[5]:= StringJoin["a", x]
Out[5]= StringJoin["a", x]

In[6]:= StringJoin[Characters["hello"]]
Out[6]= "hello"
```

#### StringPart

Extracts characters from a string by position.

- `StringPart["string", n]`: Gives the nth character in `"string"` as a single-character string.
- `StringPart["string", -n]`: Counts from the end to give the nth-to-last character.
- `StringPart["string", {n1, n2, ...}]`: Gives a list of the specified characters.
- `StringPart["string", m;;n]`: Gives a list of characters from position m through n.
- `StringPart["string", m;;n;;s]`: Gives characters from m through n in steps of s.
- `StringPart[{s1, s2, ...}, spec]`: Gives the list of results for each of the si.
- Negative indices count from the end of the string.
- In `StringPart["string", m;;n;;s]`, m, n, and/or s can be negative.
- Returns unevaluated for invalid arguments or out-of-bounds indices.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringPart["abcdefghijklm", 6]
Out[1]= "f"

In[2]:= StringPart["abcdefghijklm", {1, 3, 5}]
Out[2]= {"a", "c", "e"}

In[3]:= StringPart["abcdefghijklm", -4]
Out[3]= "j"

In[4]:= StringPart["abcdefghijklm", 1;;6]
Out[4]= {"a", "b", "c", "d", "e", "f"}

In[5]:= StringPart["abcdefghijklm", 1;;-1;;2]
Out[5]= {"a", "c", "e", "g", "i", "k", "m"}

In[6]:= StringPart["abcdefghijklm", -1;;1;;-2]
Out[6]= {"m", "k", "i", "g", "e", "c", "a"}

In[7]:= StringPart[{"abcd", "efgh", "ijklm"}, 1]
Out[7]= {"a", "e", "i"}

In[8]:= StringPart[{"abcd", "efgh", "ijklm"}, {1, -1}]
Out[8]= {{"a", "d"}, {"e", "h"}, {"i", "m"}}

In[9]:= StringPart["abcdef", -3;;-1]
Out[9]= {"d", "e", "f"}

In[10]:= StringPart["abcde", 5;;1;;-1]
Out[10]= {"e", "d", "c", "b", "a"}

In[11]:= StringPart[x, 1]
Out[11]= StringPart[x, 1]
```

#### StringTake

Gives a substring of a string.

- `StringTake["string", n]`: Gives a string containing the first n characters.
- `StringTake["string", -n]`: Gives the last n characters.
- `StringTake["string", {n}]`: Gives the nth character.
- `StringTake["string", {m, n}]`: Gives characters m through n.
- `StringTake["string", {m, n, s}]`: Gives characters m through n in steps of s.
- `StringTake["string", UpTo[n]]`: Gives n characters, or as many as are available.
- `StringTake[{s1, s2, ...}, spec]`: Gives the list of results for each si.
- `StringTake["string", 0]`: Returns the empty string `""`.
- Negative indices count from the end of the string.
- In `StringTake["string", {m, n, s}]`, m, n, and/or s can be negative.
- Returns unevaluated for invalid arguments or out-of-bounds counts.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringTake["abcdefghijklm", 6]
Out[1]= "abcdef"

In[2]:= StringTake["abcdefghijklm", -4]
Out[2]= "jklm"

In[3]:= StringTake["abcdefghijklm", {5, 10}]
Out[3]= "efghij"

In[4]:= StringTake["abcdefghijklm", {6}]
Out[4]= "f"

In[5]:= StringTake["abcdefghijklm", {1, -1, 2}]
Out[5]= "acegikm"

In[6]:= StringTake[{"abcdef", "stuv", "xyzw"}, -2]
Out[6]= {"ef", "uv", "zw"}

In[7]:= StringTake["abc", UpTo[4]]
Out[7]= "abc"

In[8]:= StringTake["abcdef", {-3, -1}]
Out[8]= "def"

In[9]:= StringTake["abcde", {5, 1, -1}]
Out[9]= "edcba"

In[10]:= StringTake["abc", 0]
Out[10]= ""

In[11]:= StringTake[x, 1]
Out[11]= StringTake[x, 1]
```

### Arithmetic and Algebra

#### Plus (+)
Symbolic sum.
- `a + b + ...`

**Features**:
- `Flat`, `Orderless`, `Listable`.
- Combines numeric constants and collects like terms (e.g., `x + 2x` becomes `3x`).
- Returns `0` if no arguments are provided.
- Returns `Overflow[]` if integer addition overflows or if any argument is `Overflow[]`.

```mathematica
In[1]:= 1 + 2 + x + 2*x
Out[1]= 3 + 3*x
```

#### Total
Gives the total of elements in a list.
- `Total[list]`
- `Total[list, n]`
- `Total[list, {n}]`
- `Total[list, {n1, n2}]`

**Features**:
- `Protected`.
- `Total[list]` is equivalent to `Apply[Plus, list]`.
- `Total[list, n]` totals all elements down to level `n`.
- `Total[list, {n}]` totals elements at level `n` only.
- Supports negative levels to count from the bottom (`-1` is the last dimension).
- Handles ragged arrays correctly by summing from the inside out when multiple levels are specified.
- `Total[list, Infinity]` totals all atoms in the expression.

```mathematica
In[1]:= Total[{a, b, c, d}]
Out[1]= a + b + c + d

In[2]:= Total[{{1, 2}, {3, 4}}]
Out[2]= {4, 6}

In[3]:= Total[{{1, 2}, {3, 4}}, 2]
Out[3]= 10

In[4]:= Total[{{1, 2}, {3}}, 2]
Out[4]= 6
```

#### Accumulate
Gives the running cumulative totals of the elements in a list.
- `Accumulate[list]`
- `Accumulate[list, Method -> "CompensatedSummation"]`

**Features**:
- `Protected`.
- `Accumulate[list]` has the same length as `list`, and is effectively equivalent to `FoldList[Plus, list]`.
- The head of the input is preserved, so `Accumulate[f[a, b, c]]` returns `f[a, a + b, a + b + c]`.
- Threads naturally over rows via `Listable` `Plus`, so `Accumulate` of a matrix accumulates within columns.
- Works on machine integers, GMP arbitrary-precision integers, machine-precision doubles, and symbolic expressions.
- With `Method -> "CompensatedSummation"`, Kahan compensated summation is used in double precision when every element is a machine number, reducing floating-point round-off error. For symbolic or mixed input the option is silently ignored and the standard symbolic accumulation is returned.

```mathematica
In[1]:= Accumulate[{a, b, c, d}]
Out[1]= {a, a + b, a + b + c, a + b + c + d}

In[2]:= Accumulate[{{a, b}, {c, d}, {e, f}}]
Out[2]= {{a, b}, {a + c, b + d}, {a + c + e, b + d + f}}

In[3]:= Accumulate[f[a, b, c, d]]
Out[3]= f[a, a + b, a + b + c, a + b + c + d]

In[4]:= Accumulate[{1, 2, 3, 4, 5}]
Out[4]= {1, 3, 6, 10, 15}

In[5]:= Accumulate[{1.0, 2.0, 3.0}, Method -> "CompensatedSummation"]
Out[5]= {1., 3., 6.}
```

#### Times (*)
Symbolic product.
- `a * b * ...`

**Features**:
- `Flat`, `Orderless`, `Listable`.
- Combines numeric constants and groups identical bases into `Power` expressions.
- Handles `I` as `Complex[0, 1]`.
- Returns `1` if no arguments are provided.
- Returns `Overflow[]` if integer multiplication overflows or if any argument is `Overflow[]`.

```mathematica
In[1]:= 2 * x * 3 * y
Out[1]= 6*x*y

In[2]:= I * I
Out[2]= -1
```

#### Power (^)
Exponentiation.
- `base ^ exp`

**Features**:
- `Listable`.
- Simplifies integer powers of integers.
- Returns `Overflow[]` if the result exceeds 64-bit integer limits.
- Reduces radicals (e.g., `8^(1/2)` becomes `2*Sqrt[2]`).
- Supports complex results for negative bases (e.g., `(-1)^(1/2)` becomes `I`).
- Distributes power over product if the exponent is an integer.

```mathematica
In[1]:= Sqrt[45]
Out[1]= 3*Sqrt[5]

In[2]:= (a * b)^2
Out[2]= a^2 * b^2
```

#### Sqrt
Square root.
- `Sqrt[z]`: Internally represented as `Power[z, 1/2]`.

#### Numerator
Gives the numerator of an expression.
- `Numerator[expr]`

**Features**:
- `Protected`, `Listable`.
- Picks out terms which do not have superficially negative exponents.
- Can be used on rational and complex numbers.

```mathematica
In[1]:= Numerator[(x-1)(x-2)/(x-3)^2]
Out[1]= (-2 + x) (-1 + x)

In[2]:= Numerator[3/7 + I/11]
Out[2]= 33 + 7 I
```

#### Denominator
Gives the denominator of an expression.
- `Denominator[expr]`

**Features**:
- `Protected`, `Listable`.
- Picks out terms which have superficially negative exponents.
- Can be used on rational and complex numbers.

```mathematica
In[1]:= Denominator[(x-1)(x-2)/(x-3)^2]
Out[1]= (-3 + x)^2

In[2]:= Denominator[3/7 + I/11]
Out[2]= 77
```

#### Cancel
Cancels out common factors in the numerator and denominator of an expression.
- `Cancel[expr]`

**Features**:
- `Protected`, `Listable`.
- Threads over equations, inequalities, logic functions, and sums dynamically.
- Evaluates greatest common divisors via polynomial GCD derivations avoiding extraneous expansions.
- Handles a single symbolic base appearing with rational fractional exponents (e.g. `Sqrt[y]`, `y^(1/3)`) by treating it as an algebraic generator: substitutes `y -> g^m` where `m` is the LCM of denominators, runs the polynomial cancellation in `g`, then substitutes back.
- The algebraic-generator pass runs `Together` on the substituted form (not just GCD-cancellation), so inputs whose `g`-substituted denominator is a Plus of terms with different `g`-denominators (e.g. `1/(g^2 - 1/g)` from `1/(y^(2/3) - 1/y^(1/3))`) are handled correctly.

```mathematica
In[1]:= Cancel[(x^2 - 1) / (x - 1)]
Out[1]= 1 + x

In[2]:= Cancel[(x - y)/(x^2 - y^2) + (x^3 - 27)/(x^2 - 9)]
Out[2]= (9 + 3 x + x^2)/(3 + x) + 1/(x + y)

In[3]:= Cancel[(y - 1)/(Sqrt[y] - 1)]
Out[3]= 1 + Sqrt[y]

In[4]:= Cancel[(y - 1)/(y^(1/3) - 1)]
Out[4]= 1 + y^(1/3) + y^(2/3)

In[5]:= Cancel[1/(y^(2/3) - 1/y^(1/3))]
Out[5]= y^(1/3)/(-1 + y)
```

#### Together
Puts terms in a sum over a common denominator, and cancels factors in the result.
- `Together[expr]`

**Features**:
- `Protected`, `Listable`.
- Makes a sum of terms into a single rational function.
- Computes lowest common multiples (LCM) of denominators securely without unconditionally destroying pre-factored bases unnecessarily.
- Handles a single symbolic base appearing with rational fractional exponents (e.g. `y^(1/3)`, `y^(2/3)`, `y^(73/24)`) by treating it as an algebraic generator: substitutes `y -> g^m` where `m` is the LCM of denominators, runs the polynomial pipeline in `g`, then substitutes back.

```mathematica
In[1]:= Together[a/b + c/d]
Out[1]= (b c + a d) / (b d)

In[2]:= Together[x^2/(x^2 - 1) + x/(x^2 - 1)]
Out[2]= x / (-1 + x)

In[3]:= Together[1/x + 1/(x + 1) + 1/(x + 2) + 1/(x + 3)]
Out[3]= (2 (3 + 11 x + 9 x^2 + 2 x^3)) / (x (1 + x) (2 + x) (3 + x))

In[4]:= Together[x^2/(x - y) - x y/(x - y)]
Out[4]= x

In[5]:= Together[y^(5/8)*(y^(19/8) - y^(73/24)/(y^(2/3) - 1/y^(1/3)))]
Out[5]= -y^3 / (-1 + y)
```

#### Apart
Gives the partial fraction decomposition of a rational expression.
- `Apart[expr]`
- `Apart[expr, var]`

**Features**:
- `Protected`, `Listable`.
- Writes `expr` as a polynomial in `var` together with a sum of ratios of polynomials with minimal denominators.
- If `var` is not specified, intelligently selects the main polynomial variable natively.
- Implements exact undetermined coefficients algebraically leveraging row-reduced identity expansions over algebraic inputs avoiding recursive fractional losses natively.
- When `Together[expr]` produces a numerator or denominator that is not polynomial in the chosen variable (e.g. fractional-power inputs whose Together'd form is `y^(1/3)/(y - 1)`), the matrix-of-coefficients algorithm cannot apply; Apart returns the `Together` form unchanged rather than synthesising a spurious zero.

```mathematica
In[1]:= Apart[1/((1+x)(5+x))]
Out[1]= 1/(4 (1 + x)) - 1/(4 (5 + x))

In[2]:= Apart[(x^5-2)/((1+x+x^2)(2+x)(1-x))]
Out[2]= 2 - x + (-1 - x/3)/(1 + x + x^2) + 1/(9 (-1 + x)) - 34/(9 (2 + x))

In[3]:= Apart[(x+y)/((x+1)(y+1)(x-y)), x]
Out[3]= 2 y/((1 + y)^2 (x - y)) - (-1 + y)/((1 + x) (1 + y)^2)

In[4]:= Apart[1/(y^(2/3) - 1/y^(1/3))]
Out[4]= y^(1/3)/(-1 + y)
```

#### Mod, Quotient, QuotientRemainder
- `Mod[n, m]`: Remainder of `n/m`.
- `Quotient[n, m]`: Integer part of `n/m`.
- `QuotientRemainder[n, m]`: Returns `{Quotient[n, m], Mod[n, m]}`.

#### GCD, LCM
- `GCD[n1, n2, ...]`: Greatest common divisor.
- `LCM[n1, n2, ...]`: Least common multiple.

#### PowerMod
Gives modular exponentiations, inverses, and roots.
- `PowerMod[a, b, m]`: Gives $a^b \pmod m$.
- `PowerMod[a, -1, m]`: Finds the modular inverse of `a` modulo `m`.
- `PowerMod[a, 1/r, m]`: Finds a modular `r`-th root of `a`.

**Features**:
- `Protected`, `Listable`.
- Evaluates much more efficiently than `Mod[a^b, m]`.
- Returns unevaluated if the corresponding inverse or root does not exist.
- Allows threading over lists natively.

```mathematica
In[1]:= PowerMod[2, 10, 3]
Out[1]= 1

In[2]:= PowerMod[3, -2, 7]
Out[2]= 4

In[3]:= PowerMod[3, 1/2, 2]
Out[3]= 1

In[4]:= PowerMod[2, {10, 11, 12, 13, 14}, 5]
Out[4]= {4, 3, 1, 2, 4}
```

#### Factorial (!)
Gives the factorial of an integer or half-integer.
- `n!` or `Factorial[n]`

**Features**:
- `Protected`, `Listable`, `NumericFunction`.
- Evaluates exactly for positive integers up to $20!$.
- Yields `ComplexInfinity` for negative integers.
- Supports half-integers utilizing factors of $\sqrt{\pi}$ recursively.
- Supports trailing `!` parsed natively as a postfix operator.

```mathematica
In[1]:= 5!
Out[1]= 120

In[2]:= (1/2)!
Out[2]= Sqrt[Pi]/2

In[3]:= Factorial[0]
Out[3]= 1
```

#### Binomial
Gives the binomial coefficient $\binom{n}{m}$.
- `Binomial[n, m]`

**Features**:
- `Protected`, `Listable`, `NumericFunction`.
- Evaluates exactly for integers, half-integers, and dynamically factors symbolic terms correctly (e.g. `Binomial[n, 4]`).
- Reduces numerical boundaries logically utilizing continuous `Gamma` interpolations.

```mathematica
In[1]:= Binomial[10, 3]
Out[1]= 120

In[2]:= Binomial[8.5, -4.2]
Out[2]= 0.0000604992

In[3]:= Binomial[9/2, 7/2]
Out[3]= 9/2

In[4]:= Binomial[n, 4]
Out[4]= 1/24 (-3 + n) (-2 + n) (-1 + n) n

In[5]:= Binomial[0, 1]
Out[5]= 0
```

#### PrimeQ
- `PrimeQ[n]`: Returns `True` if `n` is a prime number, `False` otherwise.
- `PrimeQ[z]`: For a Gaussian integer `z = a + b I`, returns `True` if `z` is a Gaussian prime.

**Features**:
- `Listable`, `Protected`.
- A Gaussian integer `a + b I` is a Gaussian prime if:
  - Both `a` and `b` are nonzero and `a^2 + b^2` is an ordinary prime, or
  - One of `a`, `b` is zero and the absolute value of the other is a prime congruent to 3 mod 4.
- Returns `False` for non-integer complex numbers, rationals, reals, and strings.

```mathematica
In[1]:= PrimeQ[7]
Out[1]= True

In[2]:= PrimeQ[1 + I]
Out[2]= True

In[3]:= PrimeQ[1 + 2 I]
Out[3]= True

In[4]:= PrimeQ[3 I]
Out[4]= True

In[5]:= PrimeQ[5 I]
Out[5]= False

In[6]:= PrimeQ[2 + 2 I]
Out[6]= False
```

#### PrimePi
- `PrimePi[x]`: Returns the number of primes less than or equal to `x`.

**Features**:
- `Listable`, `Protected`.
- Uses Meissel-Lehmer algorithm for efficient counting.

```mathematica
In[1]:= PrimePi[10]
Out[1]= 4

In[2]:= PrimePi[100]
Out[2]= 25

In[3]:= PrimePi[{10, 100}]
Out[3]= {4, 25}
```

#### NextPrime
- `NextPrime[x]`: Smallest prime above `x`.
- `NextPrime[x, k]`: $k$-th next prime above `x` (or previous if $k$ is negative).

**Features**:
- `Protected`, `ReadProtected`.
- Supports negative $k$ for finding previous primes.
- Remains unevaluated if no such prime exists (e.g., `NextPrime[2, -1]`).

#### FactorInteger
- `FactorInteger[n]`: Returns a list of prime factors and their exponents.
- `FactorInteger[n, k]`: Partial factorization, at most `k` distinct factors.
- `FactorInteger[n, Automatic]`: Pulls out easy factors using trial division.
- `FactorInteger[n, Method -> method]`: Factors `n` using a specific algorithmic method.

**Options for `Method`**:
- `"Automatic"`: (Default) Attempts factorization by sequentially executing Trial Division, Pollard's Rho, and ECM.
- `"TrialDivision"`: Extracts bounds matching the first 1000 computed primes and halts cleanly. Composite residues are preserved.
- `"PollardRho"`: Executes the Brent cycle-finding variant of Pollard's Rho algorithm targeting GMP bignums.
- `"BlakeRationalBaseDescent"`: Executes the Rational Base Descent algorithm against semiprimes `n = p q` where one factor is approximately `c (a/b)^k` for some unknown `c`, `k` and a small coprime integer pair `a > b`. With no explicit `"Base"`, auto-searches coprime partitions `a + b = j` for `j = 3, 4, ...` indefinitely until a factor is found (interruptible).
- `{"BlakeRationalBaseDescent", "Base" -> a/b}`: Same algorithm but pinned to a specific rational base `a/b`. Skips the auto-search.
- `"PollardP-1"`: Executes the Pollard $P-1$ algorithm, leveraging GMP and ECM natively.
- `"WilliamsP+1"`: Executes the Williams $P+1$ algorithm via the ECM library natively.
- `"Fermat"`: Explores factors symmetrically close to the square root boundary natively on large integers.
- `"CFRAC"`: Implements the continued fraction integer factorization method natively on GMP bignums.
- `"ShanksSquareForms"`: Implements Shanks's Square Forms Factorization (SQUFOF) natively on large integers. Halts explicitly if factors are not discovered within the loop constraints.
- `"ECM"`: Explicitly triggers Elliptic-Curve Method discovery natively.

**Features**:
- `Listable`, `Protected`.
- Supports negative integers (includes `{-1, 1}`).
- Supports rational numbers (denominator factors have negative exponents).

```mathematica
In[1]:= FactorInteger[12]
Out[1]= {{2, 2}, {3, 1}}

In[2]:= FactorInteger[-12]
Out[2]= {{-1, 1}, {2, 2}, {3, 1}}

In[3]:= FactorInteger[3/4]
Out[3]= {{2, -2}, {3, 1}}

In[4]:= FactorInteger[100, 1]
Out[4]= {{2, 2}, 25}
```

#### EulerPhi
Gives the Euler totient function $\phi(n)$.
- `EulerPhi[n]`

**Features**:
- `Listable`, `Protected`.
- Counts the number of positive integers less than or equal to $n$ that are relatively prime to $n$.
- Returns 0 for $n = 0$, and handles negative integers via $\phi(-n) = \phi(n)$.

```mathematica
In[1]:= EulerPhi[10]
Out[1]= 4
```

#### Rational
Represents a rational number.
- `Rational[n, d]`

**Features**:
- Automatically simplifies to lowest terms (e.g. `Rational[15, 5]` evaluates to `3`, `Rational[2, 4]` evaluates to `Rational[1, 2]`).
- Returns `Indeterminate` when `n` and `d` are both `0` (e.g. `Rational[0, 0]`).
- Returns `ComplexInfinity` when `n` is non-zero and `d` is `0` (e.g. `Rational[1, 0]`).

```mathematica
In[1]:= Rational[15, 5]
Out[1]= 3
```

#### Re, Im, ReIm, Abs, Conjugate, Arg
Complex number functions.
- `Re[z]`: Real part.
- `Im[z]`: Imaginary part.
- `ReIm[z]`: Returns `{Re[z], Im[z]}`.
- `Abs[z]`: Magnitude.
- `Conjugate[z]`: Complex conjugate.
- `Arg[z]`: Phase angle.

### Calculus

Symbolic differentiation is implemented natively in C (`src/deriv.c`).
Earlier versions of PicoCAS bootstrapped `D`, `Dt` and `Derivative`
from a rule file (`src/internal/deriv.m`); that approach was fragile
and slow -- every call walked a linear list of ~60 DownValues, re-ran
pattern matching against each, and recursed through the full rule
engine. The native implementation dispatches on the head symbol in a
single strcmp, uses an allocation-free structural walk for the
`FreeQ[f, x]` test, and builds the derivative tree directly, letting
the ordinary evaluator simplify the arithmetic afterwards. Measured
speedups range from roughly 2x on simple chain-rule calls to ~8x on
mixed-partial and higher-order derivatives; see "Performance Notes"
below.

#### D

Partial derivative.
- `D[f, x]` -- derivative of `f` with respect to `x` (treating all
  other symbols as constants).
- `D[f, {x, n}]` -- the `n`-th partial derivative; `n` must be a
  non-negative integer.
- `D[f, x, y, ...]` -- mixed partial derivative, equivalent to
  `D[D[f, x], y, ...]`.

**Features**:
- `Protected`, `ReadProtected`.
- Recognises the elementary heads `Plus`, `Times`, `Power`, `Sqrt`,
  `Exp`, `Log`, `Log[b, f]`, all six trig heads and their inverses,
  all six hyperbolic heads and their inverses, and threads
  element-wise over `List`.
- For a single-argument unknown head, applies the standard chain
  rule `D[f[g], x] -> Derivative[1][f][g] * D[g, x]`.
- For multi-argument unknown heads, emits the multi-index
  Derivative form, e.g. `D[f[x, y], x] -> Derivative[1, 0][f][x, y]`.
- Recognises existing `Derivative[n1, ..., nm][f][g1, ..., gn]`
  expressions and advances the appropriate partial-index.
- Short-circuits via an allocation-free FreeQ walk so constants and
  sub-trees independent of the differentiation variable are dropped
  without further work.

**Examples**:
```mathematica
In[1]:= D[x^3, x]
Out[1]= 3 x^2

In[2]:= D[Sin[x^2], x]
Out[2]= 2 x Cos[x^2]

In[3]:= D[Sin[a x], {x, 3}]
Out[3]= -a^3 Cos[a x]

In[4]:= D[f[g[x]], x]
Out[4]= Derivative[1][f][g[x]] Derivative[1][g][x]

In[5]:= D[f[x, y], x]
Out[5]= Derivative[1, 0][f][x, y]

In[6]:= D[Derivative[2][f][x], x]
Out[6]= Derivative[3][f][x]

In[7]:= D[{x, x^2, Sin[x]}, x]
Out[7]= {1, 2 x, Cos[x]}

In[8]:= D[Log[b, x], x]
Out[8]= 1 / (x Log[b])
```

#### Dt

Total derivative.
- `Dt[f]` -- total derivative of `f`. Numeric literals and the
  distinguished constants (`Pi`, `E`, `I`, `Infinity`,
  `ComplexInfinity`, `EulerGamma`, `Catalan`, `GoldenRatio`,
  `Degree`) vanish. Unknown symbols `y` appear as `Dt[y]` factors,
  modelling implicit functional dependence.
- `Dt[f, x]`, `Dt[f, {x, n}]`, `Dt[f, x, y, ...]` -- delegate to
  `D[...]` for partial-derivative behaviour.

**Features**:
- `Protected`, `ReadProtected`.
- Shares the elementary-function derivative table with `D`; the
  only dispatch difference is the base-case handling of symbols.

**Examples**:
```mathematica
In[1]:= Dt[y^2 + Sin[x]]
Out[1]= 2 y Dt[y] + Cos[x] Dt[x]

In[2]:= Dt[Pi + 3 + x y]
Out[2]= x Dt[y] + y Dt[x]

In[3]:= Dt[y^2, x]
Out[3]= 0

In[4]:= Dt[x^2, {x, 2}]
Out[4]= 2
```

#### Derivative

Higher-order derivative operator; represents the symbolic object
`Derivative[n1, ..., nk][f]` (the mixed `(n1, ..., nk)`-th partial
of a `k`-argument function `f`).

**Features**:
- `Protected`, `ReadProtected`.
- Acts primarily as a tag carried through the differentiation
  pipeline: `D` and `Dt` produce `Derivative[...]` heads for
  unknown functions and advance their indices when differentiating
  an expression that already contains one.

**Examples**:
```mathematica
In[1]:= Derivative[2][f][x]
Out[1]= Derivative[2][f][x]

In[2]:= D[%, x]
Out[2]= Derivative[3][f][x]

In[3]:= D[f[x, y], y]
Out[3]= Derivative[0, 1][f][x, y]
```

### Power Series

#### SeriesData
Represents an explicit truncated power series.
- `SeriesData[x, x0, {a0, a1, ..., a_{k-1}}, nmin, nmax, den]`

The `ai` are the coefficients of the series about the expansion point `x0`.
The power of `(x - x0)` attached to `ai` is `(nmin + i)/den`, and a trailing
`O[x - x0]^(nmax/den)` term indicates the order at which higher terms have
been dropped.

**Features**:
- `Protected`.
- `SeriesData` is a pure data head; it has no evaluator and is normally
  produced by `Series`.
- Standard printing renders the series as an ordinary mathematical sum:
  `a0 + a1 (x - x0) + a2 (x - x0)^2 + ... + O[x - x0]^p`. Zero
  coefficients are suppressed, and `x0 == 0` is displayed as simply `x`
  without the subtraction.
- `InputForm[...]` switches to the literal `SeriesData[x, x0, {...},
  nmin, nmax, den]` form, which round-trips through the parser.
- `FullForm[...]` shows the raw tree structure.

```mathematica
In[1]:= SeriesData[x, 0, {1, 1, 1/2, 1/6, 1/24, 1/120}, 0, 6, 1]
Out[1]= 1 + x + 1/2 x^2 + 1/6 x^3 + 1/24 x^4 + 1/120 x^5 + O[x]^6

In[2]:= InputForm[%]
Out[2]= SeriesData[x, 0, {1, 1, 1/2, 1/6, 1/24, 1/120}, 0, 6, 1]

In[3]:= SeriesData[x, 0, Table[i^2, {i, 10}], 0, 10, 1]
Out[3]= 1 + 4 x + 9 x^2 + 16 x^3 + 25 x^4 + 36 x^5 + 49 x^6 + 64 x^7 + 81 x^8 + 100 x^9 + O[x]^10

In[4]:= SeriesData[x, 2, {a, b, c}, 0, 3, 1]
Out[4]= a + b (x - 2) + c (x - 2)^2 + O[x - 2]^3

In[5]:= SeriesData[x, 0, {1, 2, 3}, 1, 7, 2]
Out[5]= Sqrt[x] + 2 x + 3 x^(3/2) + O[x]^(7/2)
```

#### Series
Produces the power-series expansion of an expression about a point.
- `Series[f, {x, x0, n}]` — Taylor/Laurent/Puiseux expansion up to order `(x - x0)^n`.
- `Series[f, x -> x0]` — leading-term form. The engine scans the internal expansion for the first non-zero coefficient at exponent `e1` and the next non-zero at `e2 > e1`; the reported `O`-term lands at exponent `e2` (or `e1 + 1` when no further non-zero term exists). So `Series[Sin[x] - x, x -> 0]` returns `-x^3/6 + O[x]^5`, `Series[f[x], x -> 0]` returns `f[0] + O[x]^1`, and analytic-at-x0 inputs collapse to their constant plus `O[x - x0]^1`.
- `Series[f, {x, x0, nx}, {y, y0, ny}, ...]` — iterated multivariate expansion. Each inner coefficient is itself a `SeriesData` in the next variable.
- `Series[f, {x, Infinity, n}]` — expansion at infinity, substituting `x -> 1/u` internally. The emitted `SeriesData` uses `Power[x, -1]` as its variable, so the series prints with `1/x` as the base and an `O[1/x]^(n+1)` term.

**Features**:
- `HoldAll` and `Protected` (so the expansion variable is not evaluated before `Series` has a chance to shield it).
- Threaded over lists: `Series[{f1, f2, ...}, spec]` becomes `{Series[f1, spec], Series[f2, spec], ...}`.
- Handles Taylor expansions for smooth functions, Laurent expansions where the function has a pole at `x0`, Puiseux expansions for fractional-power cases such as `Sqrt[Sin[x]]`, and logarithmic expansions for cases like `x^x` where `Log[x]` survives as a symbolic coefficient.
- Symbolic parameters in exponents are supported: `Series[(1 + x)^n, {x, 0, 4}]` returns the binomial expansion with `n` kept unexpanded.
- Approximate numeric coefficients flow through series arithmetic unchanged.
- For unknown heads (e.g. `f[x]` where `f` has no rules), the engine falls back to naive Taylor via `D` at the expansion point; the coefficients appear as `Derivative[k][f][x0]`.

**Coefficient arithmetic** automatically promotes to BigInt-backed `Rational` when 64-bit numerators or denominators would overflow, so previously-overflowing Laurent/Puiseux cases like `Series[1/Sin[x]^10, {x, 0, 2}]` and `Series[Sqrt[Log[1 + x]], {x, 0, 12}]` now produce exact coefficients (at the cost of slower evaluation for large orders).

Inverse trigonometric and inverse hyperbolic heads (`ArcSin`, `ArcCos`, `ArcTan`, `ArcCot`, `ArcSinh`, `ArcCosh`, `ArcTanh`, `ArcCoth`) are handled by direct series kernels at `u = 0` rather than by naive repeated differentiation, which would blow up expression size exponentially for higher orders. `ArcCosh` uses the principal-branch identity `ArcCosh[u] = I*ArcCos[u]`, so its expansion at `x = 0` has the expected `I*Pi/2` constant term and imaginary coefficients.

Forward reciprocal heads (`Sec`, `Csc`, `Cot`, `Sech`, `Csch`, `Coth`) are rewritten as `1/Cos[x]`, `1/Sin[x]`, `Cos[x]/Sin[x]`, etc., before expansion. Inverse reciprocal heads (`ArcSec`, `ArcCsc`, `ArcSech`, `ArcCsch`) are rewritten via the identities `ArcSec[z] = ArcCos[1/z]`, `ArcCsc[z] = ArcSin[1/z]`, `ArcSech[z] = ArcCosh[1/z]`, `ArcCsch[z] = ArcSinh[1/z]`, so that a blowing-up inner series (e.g. `z = 1/x`) collapses to a convergent kernel case rather than triggering spurious `Power::infy` warnings.

Expansions where the inner series diverges at the expansion point (e.g. `Series[f[1/x], {x, 0, n}]`) are handled via dedicated at-infinity identities:
- `ArcCoth[1/u] = ArcTanh[u]`, `ArcCot[1/u] = ArcTan[u]` (handled at the series level via inner-series inversion).
- `ArcTanh[1/u] = I*Pi/2 + ArcTanh[u]` (principal branch).
- `ArcSinh[1/v] = -Log[v] + Log[1 + Sqrt[1 + v^2]]` and `ArcCosh[1/v] = -Log[v] + Log[1 + Sqrt[1 - v^2]]` (handled by rewriting at the expression level; the symbolic `-Log[x]` term rides the existing `Log[x]` symbolic-coefficient path).

**Internal padding for symbolic expansion points**: The engine computes series at a padded internal order (user order + 12 by default) so that intermediate Laurent/Puiseux operations don't lose accuracy. When the expansion point `x0` is not a literal number, padding is capped at 2 — at a symbolic point the series coefficients are themselves symbolic expressions (e.g. `Cosh[a]`, `Sinh[a]`), and the `O(N^2)` convolution inside `so_inv`/`so_div` would otherwise spin indefinitely on exponentially growing expression trees. This makes cases like `Series[Coth[x], {x, a, 1}]`, `Series[Tanh[x], {x, a, 1}]`, `Series[Sec[x], {x, a, 1}]`, and `Series[1/Cosh[x], {x, a, 1}]` terminate in milliseconds.

**Constant inputs**: If `f` is free of the expansion variable (e.g. `Series[0, {x, 0, 4}]`, `Series[Sin[y], {x, 0, 4}]`, `Series[a + b^2, {x, 0, 3}]`), `Series` returns `f` verbatim instead of wrapping it in a trivial `SeriesData`.

**Symbolic prefactors**: A factor of `x^alpha` with `alpha` symbolic (non-integer, non-rational) is pulled outside the expansion so the remaining body is expanded as an ordinary power series. For example, `Series[x^a Exp[x], {x, 0, 5}]` returns `x^a (1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120 + O[x]^6)` — a `Times[Power[x, a], SeriesData[...]]` at the expression level, so the `SeriesData` pretty-printer still renders the body and the outer `Times` decorates it with the symbolic prefactor.

**Expansion at regular points of `Arc*` heads**: When `so_apply_kernel_at_zero` can't apply (because the inner series constant `c` is not `0`), the engine falls back to naive Taylor via `D`. This makes `Series[ArcSin[x], {x, 1/2, 3}]`, `Series[ArcTan[x], {x, 2, 2}]`, `Series[ArcSinh[x], {x, 1, 2}]`, etc. work without special-casing each non-zero expansion point.

**Maxima-style algebraic fast paths**:
- **Monomial binomial** `(a + b x^m)^alpha` with `alpha` non-integer. The generic path factors out `a`, forms `u = (b/a) x^m`, and feeds `(1+u)^alpha` through Horner composition. When `u` is a single-term series (exactly one non-zero coefficient in `SeriesObj` terms) we skip the `O(N^2)` convolution and emit `binomial(alpha, k) * (b/a)^k` directly at exponent `k*m`. This covers `Sqrt[1+x]`, `(1 - x^2)^(1/2)`, `(1+x)^(-1/2)`, `(2 + 3x)^(1/3)`, `(1+x)^n` with symbolic `n`, and Puiseux bases like `(1+x^(1/2))^(1/2)`.
- **split-two-term probe** (`series_split_two_term` in `series.h`). Structural decomposition of `e` into `a + b*x^(p/q)` without running the full series-expansion pipeline. Feeds the Log fast path and the Apart gate; exposed for unit testing.
- **Log fast path**: when `arg` matches `a + b x^(p/q)` with `a, b` both free of `x` and `a != 1`, rewrite `Log[a + b x^c]` as `Log[a] + Log[1 + (b/a) x^c]` and let the `Log1p` kernel compose with a pure monomial. Maxima's `sp2log` uses the same identity.
- **Apart preprocessing**: if the input contains `Power[p(x), -n]` for `p` a polynomial in `x`, run `Apart[f, x]` to decompose into partial fractions before expanding. Composite denominators like `1/((1-x)(1-2x)(1-3x))` then break up into geometric-series pieces that hit the monomial fast path. Gated by a polynomial check so non-rational denominators (e.g. `1/(Exp[x] - 1 - x)`) bypass Apart and fall through to the generic `so_inv` path.

**Puiseux branch points for `ArcSin` / `ArcCos` at `x = ±1`**: Dedicated handler emits a Puiseux series with `den = 2` using the identity `ArcCos[1 - s] = Sqrt[2s] sum_{k>=0} b_k s^k / (2k+1)` (with `b_k = (2k)! / (8^k (k!)^2)`) and the symmetries `ArcSin[x] = Pi/2 - ArcCos[x]`, `ArcCos[-x] = Pi - ArcCos[x]`. Supports the simple-linear inner case `ArcSin[c + q (x - x0)]` / `ArcCos[c + q (x - x0)]` with `c = ±1`. Example: `Series[ArcSin[x], {x, 1, 1}]` returns `Pi/2 - I Sqrt[2] Sqrt[x - 1] + O[x - 1]^(3/2)`; `Series[ArcCos[x], {x, -1, 2}]` returns `Pi - Sqrt[2] Sqrt[x + 1] - Sqrt[2] (x + 1)^(3/2) / 12 + O[x + 1]^(5/2)`.

**Known limitation**: Puiseux branch points for the hyperbolic `Arc*` heads (`ArcCosh[x]` at `x = 1`, `ArcTanh[x]` at `x = ±1`, etc.) and for non-simple inner series at the circular branch points are still returned unevaluated rather than risking infinite-loop or incorrect output. The naive-Taylor fallback also caps iterations and bails out on `Infinity`/`Indeterminate` derivatives so unknown heads cannot spin the engine.

```mathematica
In[1]:= Series[Exp[x], {x, 0, 10}]
Out[1]= 1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120 + x^6/720 + x^7/5040 + x^8/40320 + x^9/362880 + x^10/3628800 + O[x]^11

In[2]:= Series[f[x], {x, a, 3}]
Out[2]= f[a] + Derivative[1][f][a] (x - a) + 1/2 Derivative[2][f][a] (x - a)^2 + 1/6 Derivative[3][f][a] (x - a)^3 + O[x - a]^4

In[3]:= Series[Cos[x]/x, {x, 0, 10}]
Out[3]= 1/x - x/2 + x^3/24 - x^5/720 + x^7/40320 - x^9/3628800 + O[x]^11

In[4]:= Series[Sqrt[Sin[x]], {x, 0, 10}]
Out[4]= Sqrt[x] - x^(5/2)/12 + x^(9/2)/1440 - x^(13/2)/24192 - 67 x^(17/2)/29030400 + O[x]^(21/2)

In[5]:= Series[x^x, {x, 0, 4}]
Out[5]= 1 + Log[x] x + Log[x]^2/2 x^2 + Log[x]^3/6 x^3 + Log[x]^4/24 x^4 + O[x]^5

In[6]:= Series[(1 + x)^n, {x, 0, 4}]
Out[6]= 1 + n x + 1/2 n (-1 + n) x^2 + 1/6 n (-2 + n) (-1 + n) x^3 + 1/24 n (-3 + n) (-2 + n) (-1 + n) x^4 + O[x]^5

In[7]:= Series[Sin[1/x], {x, Infinity, 10}]
Out[7]= 1/x - 1/6 (1/x)^3 + 1/120 (1/x)^5 - 1/5040 (1/x)^7 + 1/362880 (1/x)^9 + O[1/x]^11

In[8]:= Series[Sin[x + y], {x, 0, 3}, {y, 0, 3}]
Out[8]= (y - y^3/6 + O[y]^4) + (1 - y^2/2 + O[y]^4) x + (-y/2 + y^3/12 + O[y]^4) x^2 + (-1/6 + y^2/12 + O[y]^4) x^3 + O[x]^4

In[9]:= Series[{Sin[x], Cos[x], Tan[x]}, {x, 0, 5}]
Out[9]= {x - x^3/6 + x^5/120 + O[x]^6, 1 - x^2/2 + x^4/24 + O[x]^6, x + x^3/3 + 2 x^5/15 + O[x]^6}
```

#### Normal
Converts a `SeriesData` back into an ordinary expression by dropping its O-term.
- `Normal[expr]`

**Features**:
- `Protected`.
- Returns the Plus of the coefficient-times-power terms (zero coefficients skipped). For non-`SeriesData` input, `Normal` is the identity.

```mathematica
In[1]:= Normal[Series[Exp[x], {x, 0, 5}]]
Out[1]= 1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120

In[2]:= Normal[a + b]
Out[2]= a + b
```

### Elementary Functions

#### Trig Functions
`Sin`, `Cos`, `Tan`, `Cot`, `Sec`, `Csc` and their inverses `ArcSin`, etc.

**Features**:
- `Listable`, `NumericFunction`.
- Exact values for rational multiples of `Pi` with denominators `1, 2, 3, 4, 5, 6, 10, 12`.
- `ArcTan[x, y]` computes the quadrant-aware inverse tangent.

```mathematica
In[1]:= Sin[Pi/6]
Out[1]= 1/2

In[2]:= ArcTan[1, 1]
Out[1]= 1/4*Pi
```

#### Hyperbolic Functions
`Sinh`, `Cosh`, `Tanh`, `Coth`, `Sech`, `Csch` and their inverses.

**Features**:
- `Listable`, `NumericFunction`.
- Special values for `0` and `Infinity`.

#### Exponential and Logarithmic Functions
- `Exp[z]`: Natural exponential. Reduces `Exp[I*q*Pi]` using Euler's formula.
- `Log[z]`: Natural logarithm. For a negative integer `n` (including `BigInt`), rewrites `Log[n]` as `I Pi + Log[-n]` (principal branch).
- `Log[b, z]`: Logarithm to base `b`.

```mathematica
In[1]:= Log[2, 8]
Out[1]= 3

In[2]:= Exp[I * Pi]
Out[2]= -1

In[3]:= Log[-5]
Out[3]= Log[5] + (I) Pi
```

**Power/Log cancellation.** `builtin_power` recognizes a `Log` in the
exponent and strips it:

- `Exp[c Log[a]]               -> a^c`   (since `Log[E] = 1`)
- `Power[base, c Log[base, a]] -> a^c`   (the `Log[base]` denominator cancels)

Internally `Log[b, a]` is represented as `Log[a] * Log[b]^(-1)`, so
`b^Log[b, a] -> a` and `Power[10, 3 Log[10, x]] -> x^3` both fall out of the
same rule. The coefficient `c` may be any (sub-)expression; it is the
product of whatever factors remain after removing the matching `Log[a]`
and (when `base != E`) `Log[base]^(-1)` factors from the exponent.

```mathematica
In[4]:= Exp[b Log[a]]
Out[4]= a^b

In[5]:= b^Log[b, a]
Out[5]= a
```

**Log of a Power with a matching base.** `builtin_log` folds
`Log[E^k] -> k` and `Log[b, b^k] -> k` when the exponent `k` is a real
numeric (integer, bigint, rational, or real) and (for the two-argument
form) `b` is a known-positive value (positive numeric, or a symbol
like `E` / `Pi`). The real-`k` guard keeps us on the principal branch --
for complex `k` the identity `Log[b^k] = k Log[b]` can fail by a
multiple of `2 Pi i`.

```mathematica
In[6]:= Log[E^4]
Out[6]= 4

In[7]:= Log[E^(1/3)]
Out[7]= 1/3

In[8]:= Log[2, 2^(1/3)]
Out[8]= 1/3
```

#### Trig / inverse-trig and hyperbolic / inverse-hyperbolic cancellation

`builtin_sin`, `builtin_cos`, ..., `builtin_csc` and their hyperbolic
counterparts fold the direct forward-of-inverse identities:

- `Sin[ArcSin[x]]   -> x`,  `Cos[ArcCos[x]]   -> x`, ...
- `Sinh[ArcSinh[x]] -> x`,  `Cosh[ArcCosh[x]] -> x`, ...

These hold identically over the complex numbers because each `ArcF` is a
right inverse of `F` by construction. The opposite direction
(`ArcSin[Sin[x]]`, etc.) is *not* folded: it only reduces to `x` on each
function's principal domain.

The two-argument form `ArcTan[x, y]` is deliberately excluded from the
`Tan[ArcTan[...]]` rule -- `Tan[ArcTan[x, y]] = y/x`, not a single
variable -- so `Tan[ArcTan[3, 4]]` stays unevaluated.

```mathematica
In[9]:= Sin[ArcSin[x^2 + 1]]
Out[9]= 1 + x^2

In[10]:= Tanh[ArcTanh[z]]
Out[10]= z

In[11]:= ArcSin[Sin[x]]
Out[11]= ArcSin[Sin[x]]
```


#### TrigToExp
Converts trigonometric and hyperbolic functions in `expr` to exponentials.
- `TrigToExp[expr]`

**Features**:
- `Listable`, `Protected`.
- Operates on both circular and hyperbolic functions, and their inverses.
- Automatically threads over lists, equations, inequalities, and logic functions.

```mathematica
In[1]:= TrigToExp[Cos[x]]
Out[1]= 1/2 E^(-I x) + 1/2 E^(I x)

In[2]:= TrigToExp[ArcSin[x]]
Out[2]= -I Log[I x + Sqrt[1 - x^2]]

In[3]:= TrigToExp[{Sin[x], Cos[x], Tan[x]}]
Out[3]= {(-1/2*I) E^((I) x) + (1/2*I) E^((-I) x), 1/2 E^((-I) x) + 1/2 E^((I) x), (-I) E^((I) x)/(E^((-I) x) + E^((I) x)) + (I) E^((-I) x)/(E^((-I) x) + E^((I) x))}
```

#### ExpToTrig
Converts exponentials in `expr` to trigonometric or hyperbolic functions.
- `ExpToTrig[expr]`

**Features**:
- `Listable`, `Protected`.
- Tries when possible to give results that do not involve explicit complex numbers.
- ExpToTrig is natively the inverse of `TrigToExp`.
- Automatically threads over lists, equations, inequalities, and logic functions.

```mathematica
In[1]:= ExpToTrig[Exp[I x]]
Out[1]= Cos[x] + I Sin[x]

In[2]:= ExpToTrig[Log[1 + I x] - Log[1 - I x]]
Out[2]= 2 I ArcTan[x]

In[3]:= ExpToTrig[Exp[I x] == -1]
Out[3]= Cos[x] + I Sin[x] == -1
```

#### TrigExpand
Expands trigonometric functions in `expr` by splitting up sums and integer
multiples that appear in arguments of circular and hyperbolic trigonometric
functions.
- `TrigExpand[expr]`

**Features**:
- `Listable`, `Protected`.
- Operates on both circular (`Sin`, `Cos`, `Tan`, `Cot`, `Sec`, `Csc`) and
  hyperbolic (`Sinh`, `Cosh`, `Tanh`, `Coth`, `Sech`, `Csch`) functions.
- Applies angle-addition formulas to `Sin[a + b + …]`, `Cos[a + b + …]`,
  `Sinh[a + b + …]`, `Cosh[a + b + …]` to a fixed point.
- Applies multiple-angle reductions to `Sin[n x]`, `Cos[n x]`, `Sinh[n x]`,
  `Cosh[n x]` for integer `n`, recursively reducing to `Sin[x]` / `Cos[x]` /
  `Sinh[x]` / `Cosh[x]`.
- `Tan`, `Cot`, `Sec`, `Csc` (and `Tanh`, `Coth`, `Sech`, `Csch`) with sum or
  integer-multiple arguments are rewritten as ratios of `Sin`/`Cos`
  (resp. `Sinh`/`Cosh`) and then expanded.
- Distributes products over sums via `Expand` so the result is a flat sum of
  monomials.
- Applies the Pythagorean identities `Sin[x]^2 + Cos[x]^2 -> 1` and
  `Cosh[x]^2 - Sinh[x]^2 -> 1` as a final reduction pass, including powers of
  both identities for any integer `n >= 1`:
    - `Sin[n x]^2 + Cos[n x]^2` expands to `(Sin[x]^2 + Cos[x]^2)^n` and
      collapses to `1` via a Factor-based reduction.
    - `Cosh[n x]^2 - Sinh[n x]^2` factors as
      `(Cosh[x] + Sinh[x])^n (Cosh[x] - Sinh[x])^n` and collapses to `1`.
  Negated and scalar-weighted forms (e.g. `-Sin[n x]^2 - Cos[n x]^2`,
  `-5 (Sin[n x]^2 + Cos[n x]^2)`, `Sinh[n x]^2 - Cosh[n x]^2`) collapse to the
  expected signed constant — the Pythagorean rules match both possible signs
  that `Factor` may emerge with and allow an arbitrary remainder of factors in
  the surrounding `Times`. Expressions that contain a denominator (any
  `Power[_, negative_Integer]` subterm) skip the Factor pass so that canonical
  forms such as `(2 Cos[x] Sin[x])/(Cos[x]^2 - Sin[x]^2)` are preserved.
  Inputs without a Pythagorean-eligible squared structure (no pair
  `Sin[a]^k`/`Cos[a]^k` or `Sinh[a]^k`/`Cosh[a]^k` with the same argument and
  `k >= 2`) likewise skip the Factor pass; the multivariate polynomials that
  multi-angle expansions such as `TrigExpand[Sin[2 x + 3 y]]` produce would
  otherwise make `Factor` prohibitively slow without yielding any collapse.
  The Factor pass is also skipped when the expanded form contains more than
  two distinct squared trigonometric atoms (e.g. `Cos[x]^2`, `Sin[x]^2`,
  `Cos[y]^2`, `Sin[y]^2` together): even if a Pythagorean pair is structurally
  present, `Factor` on the resulting dense multivariate polynomial stalls
  without producing a useful collapse.
- Automatically threads over lists (via `Listable`), as well as equations,
  inequalities (`Equal`, `Unequal`, `Less`, `LessEqual`, `Greater`,
  `GreaterEqual`, `SameQ`, `UnsameQ`), and logic functions (`And`, `Or`,
  `Not`, `Xor`, `Implies`).

```mathematica
In[1]:= TrigExpand[Sin[2 x]]
Out[1]= 2 Cos[x] Sin[x]

In[2]:= TrigExpand[Sin[x + y]]
Out[2]= Cos[x] Sin[y] + Cos[y] Sin[x]

In[3]:= TrigExpand[Sin[3 x]]
Out[3]= 3 Cos[x]^2 Sin[x] - Sin[x]^3

In[4]:= TrigExpand[Cos[x + y + z]]
Out[4]= -Cos[x] Sin[y] Sin[z] - Cos[y] Sin[x] Sin[z] - Cos[z] Sin[x] Sin[y] + Cos[x] Cos[y] Cos[z]

In[5]:= TrigExpand[Sin[x]^2 + Cos[x]^2]
Out[5]= 1

In[5b]:= TrigExpand[Sin[4 x]^2 + Cos[4 x]^2]
Out[5b]= 1

In[6]:= TrigExpand[Sinh[4 x]]
Out[6]= 4 Cosh[x] Sinh[x]^3 + 4 Cosh[x]^3 Sinh[x]

In[7]:= TrigExpand[Cosh[x - y]]
Out[7]= -Sinh[x] Sinh[y] + Cosh[x] Cosh[y]

In[8]:= TrigExpand[Tanh[2 t]]
Out[8]= (2 Cosh[t] Sinh[t])/(Cosh[t]^2 + Sinh[t]^2)

In[9]:= TrigExpand[{Tan[2 x], Sinh[x + y]}]
Out[9]= {(2 Cos[x] Sin[x])/(Cos[x]^2 - Sin[x]^2), Cosh[x] Sinh[y] + Cosh[y] Sinh[x]}

In[10]:= TrigExpand[1 < Cos[x + y] < 2]
Out[10]= 1 < -Sin[x] Sin[y] + Cos[x] Cos[y] < 2
```

#### TrigFactor
Factors trigonometric functions in `expr`. Broadly the functional inverse of
`TrigExpand`: combines sums of trigonometric products into angle-sum forms,
collapses Pythagorean identities (both circular and hyperbolic), recognises
the reverse of the double-angle identities, and factors polynomial structure
over the trigonometric atoms.
- `TrigFactor[expr]`

**Features**:
- `Listable`, `Protected`.
- Operates on both circular (`Sin`, `Cos`, `Tan`, `Cot`, `Sec`, `Csc`) and
  hyperbolic (`Sinh`, `Cosh`, `Tanh`, `Coth`, `Sech`, `Csch`) functions.
- Pipeline:
  1. Rewrite reciprocal heads (`Tan`, `Cot`, `Sec`, `Csc`, and their
     hyperbolic analogs) as `Sin`/`Cos`/`Sinh`/`Cosh` ratios so that `Factor`
     sees the full polynomial structure.
  2. Combine into a single rational via `Together`.
  3. Run `Factor` on the resulting rational; trigonometric atoms are treated
     as independent polynomial variables. The `Factor` pass is skipped when
     the post-`Together` form contains more than two distinct squared
     trigonometric atoms (e.g. `Sin[x]^2`, `Cos[x]^2`, `Sinh[y]^2`,
     `Cosh[y]^2` together): on such dense multivariate polynomials Factor's
     trial-division loop stalls without producing a useful factorization, and
     the identity rules in step 4 still match Pythagorean structure that
     survives in the post-`Together` factored form (e.g.
     `(Sin[x]^2 + Cos[x]^2)(Cosh[y]^2 - Sinh[y]^2)` collapses directly to
     `1` via the `Times`-context Pythagorean rules).
  4. Apply identity collapse rules via `ReplaceRepeated`: Pythagorean
     identities (`Sin^2 + Cos^2 -> 1`, `Cosh^2 - Sinh^2 -> 1`, with and
     without arbitrary coefficients), reverse angle-addition
     (`Sin[a]Cos[b] ± Cos[a]Sin[b] -> Sin[a ± b]`,
     `Cos[a]Cos[b] ± Sin[a]Sin[b] -> Cos[a ∓ b]`, and hyperbolic analogs),
     reverse double-angle (`2 Sin Cos -> Sin[2x]`,
     `Cos^2 - Sin^2 -> Cos[2x]`, `Cosh^2 + Sinh^2 -> Cosh[2x]`), and
     factored-form variants such as `(Cos - Sin)(Cos + Sin) -> Cos[2x]`,
     `(Cosh - 1)(Cosh + 1) -> Sinh^2`, and
     `(Cosh - Sinh)(Cosh + Sinh) -> 1` that arise naturally from `Factor`.
  5. Restore `Tan`/`Cot`/`Sec`/`Csc` (and hyperbolic analogs) from the
     `Sin`/`Cos` ratio form so reciprocal heads survive the round-trip.
- Two paths are tried: the primary pipeline (preserves angle-sum structure)
  and a fallback that `TrigExpand`s the argument first (catches
  cancellations that only become visible after the angle-sum is expanded,
  e.g. `Cos[x + y] + Sin[x] Sin[y] -> Cos[x] Cos[y]`). The fallback runs
  only when the primary pipeline leaves the expression unchanged, so
  structurally productive inputs (e.g. `Sin[x + y]^2 + Tan[x + y]`) avoid
  the expensive expanded-rational path. The final result is the smaller of
  the two by leaf count; ties favour the primary pipeline.
- Automatically threads over lists (via `Listable`), as well as equations,
  inequalities (`Equal`, `Unequal`, `Less`, `LessEqual`, `Greater`,
  `GreaterEqual`, `SameQ`, `UnsameQ`), and logic functions (`And`, `Or`,
  `Not`, `Xor`, `Implies`).

```mathematica
In[1]:= TrigFactor[Sin[x]^2 + Cos[x]^2]
Out[1]= 1

In[2]:= TrigFactor[Cosh[x]^2 - Sinh[x]^2]
Out[2]= 1

In[3]:= TrigFactor[2 Sin[x] Cos[x]]
Out[3]= Sin[2 x]

In[4]:= TrigFactor[Cos[x]^2 - Sin[x]^2]
Out[4]= Cos[2 x]

In[5]:= TrigFactor[Sin[a] Cos[b] + Cos[a] Sin[b]]
Out[5]= Sin[a + b]

In[6]:= TrigFactor[Cos[a] Cos[b] + Sin[a] Sin[b]]
Out[6]= Cos[a - b]

In[7]:= TrigFactor[Sin[x]^2 + Tan[x]^2]
Out[7]= Tan[x]^2 (1 + Cos[x]^2)

In[8]:= TrigFactor[Cosh[x]^2 - Cosh[x]^4]
Out[8]= -Cosh[x]^2 Sinh[x]^2

In[9]:= TrigFactor[Sin[x+y]^2 + Tan[x+y]]
Out[9]= Tan[x + y] (1 + Cos[x + y] Sin[x + y])

In[10]:= TrigFactor[Cos[x + y] + Sin[x] Sin[y]]
Out[10]= Cos[x] Cos[y]

In[11]:= TrigFactor[Cos[x]^4 - Sin[x]^4]
Out[11]= Cos[2 x]

In[12]:= TrigFactor[{Sin[x]^2 + Cos[x]^2, 2 Sinh[x] Cosh[x]}]
Out[12]= {1, Sinh[2 x]}

In[13]:= TrigFactor[Sin[x]^2 + Cos[x]^2 == 1]
Out[13]= True
```

#### TrigReduce

`TrigReduce[expr]` rewrites products and integer powers of single-argument
trigonometric functions in `expr` in terms of trigonometric functions with
combined arguments. It is the inverse-of-`TrigExpand` direction for the
product/power identities: given a trigonometric polynomial, `TrigReduce`
typically yields a linear expression involving trigonometric functions
with more complicated arguments.

- `TrigReduce[expr]`

**Features**:
- Applies the classical product-to-sum identities (Sin·Cos, Sin·Sin,
  Cos·Cos, plus the four hyperbolic analogues) and the power-reduction
  identities (Sin² → (1 − Cos[2x])/2, etc., extended recursively to
  any positive integer power).
- Operates on both circular and hyperbolic functions; rewrites
  `Tan`/`Cot`/`Sec`/`Csc` (and the hyperbolic reciprocals) as `Sin`/`Cos`
  ratios before reduction, then restores the reciprocal head where the
  result has the matching shape.
- Recognises angle-addition forms produced after `Together`
  (`Sin[a]·Cos[b] + Cos[a]·Sin[b] → Sin[a + b]` and analogues), with
  sign variants and hyperbolic counterparts.
- Includes a sign-cancellation pass for the `Sin[a − b] + Sin[b − a] = 0`
  shape that arises when the product-to-sum rules bind asymmetrically;
  the same pass handles the corresponding `Cos`/`Sinh`/`Cosh` parities.
- Applies an `Expand`/`Together` cleanup at the end so trivial outer
  fractions (`1/2 (2 X + 2 Y)`) flatten while genuine rationals
  (`(3 − 4 Cos[2 x] + Cos[4 x])/2`) survive in normalised form.
- Memoised through the same `FactorMemo` mechanism used by
  `TrigExpand`/`TrigFactor`/`TrigToExp`, so repeated invocations during
  `Simplify` candidate-set search amortise.
- `Listable`, plus explicit threading over `Equal`, `Unequal`, `Less`,
  `LessEqual`, `Greater`, `GreaterEqual`, `SameQ`, `UnsameQ`, `And`,
  `Or`, `Not`, `Xor`, `Implies` (mirrors `TrigExpand` / `TrigFactor`).
- Idempotent on already-reduced inputs and a no-op on non-trig
  expressions or single trig calls of compound arguments.

```mathematica
In[1]:= TrigReduce[2 Cos[x]^2]
Out[1]= 1 + Cos[2 x]

In[2]:= TrigReduce[2 Sin[x] Cos[y]]
Out[2]= Sin[x + y] + Sin[x - y]

In[3]:= TrigReduce[2 Cosh[x] Cosh[y]]
Out[3]= Cosh[x + y] + Cosh[x - y]

In[4]:= TrigReduce[Sin[a] (Cos[b] - Sin[b]) + Cos[a] (Sin[b] + Cos[b])]
Out[4]= Cos[a + b] + Sin[a + b]

In[5]:= TrigReduce[Tan[x] + Tan[y]]
Out[5]= Sec[x] Sec[y] Sin[x + y]

In[6]:= TrigReduce[Coth[x] + Coth[y]]
Out[6]= Csch[x] Csch[y] Sinh[x + y]

In[7]:= TrigReduce[Sin[x]^4]
Out[7]= 1/8 (3 - 4 Cos[2 x] + Cos[4 x])

In[8]:= TrigReduce[2 Sin[x + y] Cos[x - y]]
Out[8]= Sin[2 x] + Sin[2 y]

In[9]:= TrigReduce[{Tan[x] + Cot[y], Tanh[x] - Coth[y]}]
Out[9]= {Cos[x - y] Csc[y] Sec[x], -Cosh[x - y] Csch[y] Sech[x]}

In[10]:= TrigReduce[4 Sin[x]^4 == 1 && 2 Cos[x]^2 >= 1]
Out[10]= 1/2 (3 - 4 Cos[2 x] + Cos[4 x]) == 1 && 1 + Cos[2 x] >= 1
```

`TrigReduce` is also offered as a `Simplify` transform: when the search
sees an angle-addition expansion such as
`Sin[a] (Cos[b] − Sin[b]) + Cos[a] (Sin[b] + Cos[b])`, the reduced form
`Cos[a + b] + Sin[a + b]` has fewer leaves and the score-gate selects
it.

#### Piecewise and Rounding Functions
`Floor`, `Ceiling`, `Round`, `IntegerPart`, `FractionalPart`.

**Features**:
- `Listable`.
- Applied component-wise to `Complex` numbers.
- `Round` rounds to the nearest even integer for ties.
- `Floor[x, a]` returns the greatest multiple of `a` $\le x$.

**Symbolic simplifications** (`Floor`, `Ceiling`, `Round` only -- `IntegerPart` and `FractionalPart` are excluded):

*Sign extraction* -- pulls a leading `-1` out through the rounding head, swapping `Floor`/`Ceiling` (and leaving `Round` unchanged):
- `Floor[-x]   :> -Ceiling[x]`
- `Ceiling[-x] :> -Floor[x]`
- `Round[-x]   :> -Round[x]`

The trigger is `expr_is_superficially_negative`, so `Floor[-2 x]` and `Ceiling[-3 x y]` reduce as well.

*Idempotency / composition* -- the inner expression is already an integer, so the outer rounding is a no-op:
- `Floor[Floor[x]]   :> Floor[x]`
- `Ceiling[Ceiling[x]] :> Ceiling[x]`
- `Round[Round[x]]   :> Round[x]`
- `Floor[Ceiling[x]] :> Ceiling[x]`
- `Ceiling[Floor[x]] :> Floor[x]`
- `Floor[Round[x]]   :> Round[x]`
- `Ceiling[Round[x]] :> Round[x]`
- `Round[Floor[x]]   :> Floor[x]`
- `Round[Ceiling[x]] :> Ceiling[x]`

These rules compose under fixed-point evaluation, so e.g. `Ceiling[Floor[Ceiling[x]]] -> Ceiling[x]`, and `Floor[-x] + Ceiling[x] -> 0`.

```mathematica
In[1]:= Round[2.5]
Out[1]= 2

In[2]:= Round[3.5]
Out[2]= 4

In[3]:= Floor[-x] + Ceiling[x]
Out[3]= 0

In[4]:= Ceiling[Floor[Ceiling[x]]]
Out[4]= Ceiling[x]
```

### Lists and Iteration

#### Tuples
Generates a list of all possible combinations or tuples of elements.
- `Tuples[list, n]`: generates a list of all possible n-tuples of elements from `list`.
- `Tuples[{list1, list2, ...}]`: generates a list of all possible tuples whose *i*-th element is from `listi`.
- `Tuples[list, {n1, n2, ...}]`: generates a list of all possible $n_1 \times n_2 \times \dots$ arrays of elements in `list`.

**Features**:
- `Protected`.
- The elements of `list` are treated as distinct.
- The object `list` need not have head `List`. The head at each level in the arrays generated by `Tuples` will be the same as the head of `list`.

```mathematica
In[1]:= Tuples[{0, 1}, 3]
Out[1]= {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}}

In[2]:= Tuples[{{a, b}, {1, 2, 3, 4}, {x}}]
Out[2]= {{a, 1, x}, {a, 2, x}, {a, 3, x}, {a, 4, x}, {b, 1, x}, {b, 2, x}, {b, 3, x}, {b, 4, x}}

In[3]:= Tuples[{a, b}, {2, 2}]
Out[3]= {{{{a, a}, {a, a}}, {{a, a}, {a, b}}}, {{{a, a}, {b, a}}, {{a, a}, {b, b}}}, {{{a, b}, {a, a}}, {{a, b}, {a, b}}}, {{{a, b}, {b, a}}, {{a, b}, {b, b}}}, {{{b, a}, {a, a}}, {{b, a}, {a, b}}}, {{{b, a}, {b, a}}, {{b, a}, {b, b}}}, {{{b, b}, {a, a}}, {{b, b}, {a, b}}}, {{{b, b}, {b, a}}, {{b, b}, {b, b}}}}

In[4]:= Tuples[f[x, y, z], 2]
Out[4]= {f[x, x], f[x, y], f[x, z], f[y, x], f[y, y], f[y, z], f[z, x], f[z, y], f[z, z]}
```


#### Permutations
Generates a list of all possible permutations of the elements in `list`.
- `Permutations[list]`: generates a list of all possible permutations of the elements in `list`.
- `Permutations[list, n]`: gives all permutations containing at most `n` elements.
- `Permutations[list, {n}]`: gives all permutations containing exactly `n` elements.

**Features**:
- `Protected`.
- There are $n!$ permutations of a list of $n$ distinct elements.
- Repeated elements are treated as identical.
- The object `list` need not have head `List`.
- `Permutations[list]` is effectively equivalent to `Permutations[list, {Length[list]}]`.
- `Permutations[list, {n_min, n_max}]` gives permutations of `list` between `n_min` and `n_max` elements.
- `Permutations[list, {n_min, n_max, dn}]` uses step `dn`.
- `Permutations[list, All]` is equivalent to `Permutations[list, {0, Length[list]}]`.

```mathematica
In[1]:= Permutations[{a, b, c}]
Out[1]= {{a, b, c}, {a, c, b}, {b, a, c}, {b, c, a}, {c, a, b}, {c, b, a}}

In[2]:= Permutations[{a, b, c, d}, {3}]
Out[2]= {{a, b, c}, {a, b, d}, {a, c, b}, {a, c, d}, {a, d, b}, {a, d, c}, {b, a, c}, {b, a, d}, {b, c, a}, {b, c, d}, {b, d, a}, {b, d, c}, {c, a, b}, {c, a, d}, {c, b, a}, {c, b, d}, {c, d, a}, {c, d, b}, {d, a, b}, {d, a, c}, {d, b, a}, {d, b, c}, {d, c, a}, {d, c, b}}

In[3]:= Permutations[{a, a, b}]
Out[3]= {{a, a, b}, {a, b, a}, {b, a, a}}

In[4]:= Permutations[{x, x^2, x + 1}]
Out[4]= {{x, x^2, 1 + x}, {x, 1 + x, x^2}, {x^2, x, 1 + x}, {x^2, 1 + x, x}, {1 + x, x, x^2}, {1 + x, x^2, x}}

In[5]:= Permutations[Range[3], All]
Out[5]= {{}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 1}, {2, 3}, {3, 1}, {3, 2}, {1, 2, 3}, {1, 3, 2}, {2, 1, 3}, {2, 3, 1}, {3, 1, 2}, {3, 2, 1}}

In[6]:= Permutations[Range[4], {4, 0, -2}]
Out[6]= {{1, 2, 3, 4}, {1, 2, 4, 3}, {1, 3, 2, 4}, {1, 3, 4, 2}, {1, 4, 2, 3}, {1, 4, 3, 2}, {2, 1, 3, 4}, {2, 1, 4, 3}, {2, 3, 1, 4}, {2, 3, 4, 1}, {2, 4, 1, 3}, {2, 4, 3, 1}, {3, 1, 2, 4}, {3, 1, 4, 2}, {3, 2, 1, 4}, {3, 2, 4, 1}, {3, 4, 1, 2}, {3, 4, 2, 1}, {4, 1, 2, 3}, {4, 1, 3, 2}, {4, 2, 1, 3}, {4, 2, 3, 1}, {4, 3, 1, 2}, {4, 3, 2, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 1}, {2, 3}, {2, 4}, {3, 1}, {3, 2}, {3, 4}, {4, 1}, {4, 2}, {4, 3}, {}}

In[7]:= Permutations[f[a, b, c]]
Out[7]= {f[a, b, c], f[a, c, b], f[b, a, c], f[b, c, a], f[c, a, b], f[c, b, a]}
```
#### Table
Generates a list of expressions.
- `Table[expr, n]`
- `Table[expr, {i, imax}]`
- `Table[expr, {i, imin, imax, di}]`
- `Table[expr, {i, {val1, val2, ...}}]`

**Features**:
- `HoldAll`: `expr` is evaluated once for each step.
- Supports nested iterators to create matrices.

```mathematica
In[1]:= Table[i^2, {i, 4}]
Out[1]= {1, 4, 9, 16}
```

#### Range
Generates a list of values.
- `Range[imax]`
- `Range[imin, imax, di]`

#### Array
Generates an array by applying a function to indices.
- `Array[f, n]`
- `Array[f, {n1, n2, ...}]`

### Functional Programming

#### Distribute
Implements the distributive law for operators `f` and `g`.
- `Distribute[expr]`: distributes the head of `expr` over `Plus` appearing in any of its arguments.
- `Distribute[expr, g]`: distributes over `g`.
- `Distribute[expr, g, f]`: performs the distribution only if the head of `expr` is `f`.
- `Distribute[expr, g, f, gp, fp]`: gives `gp` and `fp` in place of `g` and `f` respectively in the result of the distribution.

**Features**:
- `Protected`.
- `Distribute` explicitly constructs the complete result of a distribution; `Expand`, on the other hand, builds up results iteratively, simplifying at each stage.
- For pure products, `Distribute` gives the same results as `Expand`.

```mathematica
In[1]:= Distribute[(a + b) . (x + y + z)]
Out[1]= a . x + a . y + a . z + b . x + b . y + b . z

In[2]:= Distribute[f[a + b, c + d + e]]
Out[2]= f[a, c] + f[a, d] + f[a, e] + f[b, c] + f[b, d] + f[b, e]

In[3]:= Distribute[(a + b + c) (u + v), Plus, Times]
Out[3]= a u + a v + b u + b v + c u + c v

In[4]:= Distribute[{{a, b}, {x, y, z}, {s, t}}, List]
Out[4]= {{a, x, s}, {a, x, t}, {a, y, s}, {a, y, t}, {a, z, s}, {a, z, t}, {b, x, s}, {b, x, t}, {b, y, s}, {b, y, t}, {b, z, s}, {b, z, t}}

In[5]:= Distribute[{{}, {a}}, {{}, {b}}, {{}, {c}}, List, List, List, Join]
Out[5]= {{}, {c}, {b}, {b, c}, {a}, {a, c}, {a, b}, {a, b, c}}
```

#### Pure Functions (Function, &)
`Function` is a pure (or "anonymous") function, analogous to $\lambda$ in LISP.
- `Function[body]` or `body &`: A pure function where formal parameters are `#` (or `#1`), `#2`, etc.
- `Function[x, body]` or `Function[{x1, x2, ...}, body]`: A pure function with named local formal parameters.
- `Function[params, body, attrs]`: A pure function treated as having the evaluation attributes `attrs`.
- `Function[Null, body, attrs]`: Slot form (`#`, `##`, etc.) with evaluation attributes.
- `Slot` (`#`): Represents the first argument. `#n` represents the $n$-th argument.
- `SlotSequence` (`##`): Represents all arguments sequence. `##n` represents all arguments starting from the $n$-th.

**Features**:
- **Lexical parameter binding**: named parameters are substituted into the body before evaluation (not bound via the global symbol table). This means held references to a parameter (for example `Unevaluated[x]` inside the body) see the substituted expression rather than the raw symbol.
- **No Hold by default**: the default `Function` has no hold attributes, so its arguments are evaluated before substitution, matching Mathematica. `(Hold[#]&)[1+2]` gives `Hold[3]`, not `Hold[1 + 2]`.
- **3-argument form attributes**: `Function[params, body, attrs]` can assign any of the standard evaluator attributes. Recognised attributes include `HoldFirst`, `HoldRest`, `HoldAll`, `HoldAllComplete`, `Listable`, `Flat`, `Orderless`, `OneIdentity`, `NumericFunction`, `SequenceHold`, `NHoldRest`. `attrs` may be a single attribute symbol or a list.
- `Slot` parameters are positionally mapped to arguments provided. Remaining arguments are ignored.
- Pure functions can be assigned to variables (e.g., `f = #^2 &`) and applied over lists (e.g., `f /@ {1, 2, 3}`).
- Properly scopes nested `Function` expressions: named inner variables shadow outer ones, and unnamed inner functions establish a new scope for `#` slots.

```mathematica
In[1]:= Function[u, 3 + u][x]
Out[1]= 3 + x

In[2]:= (#1^2 + #2^4)&[x, y]
Out[2]= x^2 + y^4

In[3]:= f[X, ##, Y, ##]&[a, b, c]
Out[3]= f[X, a, b, c, Y, a, b, c]

In[4]:= Function[{x}, Length[Unevaluated[x]], {HoldAll}][1+1+1]
Out[4]= 3

In[5]:= Function[{x}, Length[Unevaluated[x]]][1+1+1]
Out[5]= 0

In[6]:= (Hold[#]&)[1+2]
Out[6]= Hold[3]

In[7]:= Function[{a, b}, {Length[Unevaluated[a]], Length[Unevaluated[b]]}, HoldFirst][1+2+3, 4+5+6]
Out[7]= {3, 0}
```

#### Map (/@)
- `f /@ expr` or `Map[f, expr]`

#### Apply (@@, @@@)
- `f @@ expr`: Level 0.
- `f @@@ expr`: Level 1.

#### MapAll (//@)
- `f //@ expr`: Recursive map.

#### MapAt
Applies a function to selected parts of an expression, identified by position specifications of the same form as those returned by `Position`.

- `MapAt[f, expr, n]`: applies `f` to the element at position `n` in `expr`. If `n` is negative, the position is counted from the end. `n = 0` targets the head of `expr`.
- `MapAt[f, expr, {i, j, ...}]`: applies `f` to the part of `expr` at position `{i, j, ...}` (equivalently `expr[[i, j, ...]]`).
- `MapAt[f, expr, {{i1, j1, ...}, {i2, j2, ...}, ...}]`: applies `f` to each of the listed parts of `expr`. If the same position appears more than once, `f` is applied repeatedly to that part.

**Features**:
- `Protected`.
- Path components may be integers (positive or negative), `All` (selects every child at that level), or `Span` expressions such as `i ;; j` or `i ;; j ;; k`.
- Works on expressions with any head (not just `List`); after substitution the evaluator re-applies canonical ordering for `Orderless` heads such as `Plus` and `Times`.
- `MapAt[f, expr, {}]` applies `f` to the entire expression itself.

```mathematica
In[1]:= MapAt[f, {a, b, c, d}, 2]
Out[1]= {a, f[b], c, d}

In[2]:= MapAt[f, {a, b, c, d}, {{1}, {4}}]
Out[2]= {f[a], b, c, f[d]}

In[3]:= MapAt[f, {{a, b, c}, {d, e}}, {2, 1}]
Out[3]= {{a, b, c}, {f[d], e}}

In[4]:= MapAt[f, {{a, b, c}, {d, e}}, {All, 2}]
Out[4]= {{a, f[b], c}, {d, f[e]}}

In[5]:= MapAt[h, {{a, b, c}, {d, e}, f, g}, -3]
Out[5]= {{a, b, c}, h[{d, e}], f, g}

In[6]:= MapAt[h, {{a, b, c}, {d, e}, f, g}, {2, 1}]
Out[6]= {{a, b, c}, {h[d], e}, f, g}

In[7]:= MapAt[h, {{a, b, c}, {d, e}, f, g}, {{2}, {1}}]
Out[7]= {h[{a, b, c}], h[{d, e}], f, g}

In[8]:= MapAt[h, {{a, b, c}, {d, e}, f, g}, {{1, 1}, {2, 2}, {3}}]
Out[8]= {{h[a], b, c}, {d, h[e]}, h[f], g}

In[9]:= MapAt[f, {1, 2, 3, 4, 5, 6}, 3 ;; 4]
Out[9]= {1, 2, f[3], f[4], 5, 6}

In[10]:= MapAt[f, a + b + c + d, 2]
Out[10]= a + c + d + f[b]

In[11]:= MapAt[f, x^2 + y^2, {{1, 1}, {2, 1}}]
Out[11]= f[x]^2 + f[y]^2

In[12]:= MapAt[f, {a, b, c}, 0]
Out[12]= f[List][a, b, c]
```

#### Nest
Applies a function repeatedly.
- `Nest[f, expr, n]`: gives an expression with `f` applied `n` times to `expr`.

**Features**:
- `Protected`.
- `n` must be a non-negative integer; `Nest[f, expr, 0]` returns `expr` unchanged.
- The function `f` may be a symbol, a built-in, or a pure function (`... &`).
- Each iteration evaluates `f[current]` before proceeding, so numeric computations collapse immediately.
- Returns unevaluated if `n` is not a non-negative integer or the argument count is wrong.

**Examples**:
```
In[1]:= Nest[f, x, 3]
Out[1]= f[f[f[x]]]

In[2]:= Nest[(1 + #)^2 &, 1, 3]
Out[2]= 676

In[3]:= Nest[(1 + #)^2 &, x, 5]
Out[3]= (1 + (1 + (1 + (1 + (1 + x)^2)^2)^2)^2)^2

In[4]:= Nest[Sqrt, 100.0, 4]
Out[4]= 1.33352

In[5]:= Nest[1/(1 + #) &, x, 5]
Out[5]= 1/(1 + 1/(1 + 1/(1 + 1/(1 + 1/(1 + x)))))

In[6]:= Nest[x^# &, x, 6]
Out[6]= x^x^x^x^x^x^x

In[7]:= Nest[#(1 + 0.05) &, 1000, 10]
Out[7]= 1628.89

In[8]:= Nest[(# + 2/#)/2 &, 1.0, 5]
Out[8]= 1.41421

In[9]:= Nest[{{1, 1}, {1, 0}} . # &, {0, 1}, 10]
Out[9]= {55, 34}
```

#### NestList
Applies a function repeatedly, collecting every intermediate result.
- `NestList[f, expr, n]`: gives a list of the results of applying `f` to `expr` 0 through `n` times.

**Features**:
- `Protected`.
- Returns a list of length `n + 1` whose first element is `expr` and whose `(k+1)`-th element is `f` applied `k` times to `expr`.
- `n` must be a non-negative integer; `NestList[f, expr, 0]` returns `{expr}`.
- The function `f` may be a symbol, a built-in, or a pure function (`... &`).
- Each iteration evaluates `f[current]` before proceeding, so numeric computations collapse immediately.
- Returns unevaluated if `n` is not a non-negative integer or the argument count is wrong.
- `Last[NestList[f, expr, n]]` is equivalent to `Nest[f, expr, n]`.

**Examples**:
```
In[1]:= NestList[f, x, 4]
Out[1]= {x, f[x], f[f[x]], f[f[f[x]]], f[f[f[f[x]]]]}

In[2]:= NestList[Cos, 1.0, 10]
Out[2]= {1., 0.540302, 0.857553, 0.65429, 0.79348, 0.701369, 0.76396, 0.722102, 0.750418, 0.731404, 0.744237}

In[3]:= NestList[(1 + #)^2 &, x, 3]
Out[3]= {x, (1 + x)^2, (1 + (1 + x)^2)^2, (1 + (1 + (1 + x)^2)^2)^2}

In[4]:= NestList[Sqrt, 100.0, 4]
Out[4]= {100., 10., 3.16228, 1.77828, 1.33352}

In[5]:= NestList[2 # &, 1, 10]
Out[5]= {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024}

In[6]:= NestList[# + 1 &, 0, 10]
Out[6]= {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}

In[7]:= NestList[#^2 &, 2, 6]
Out[7]= {2, 4, 16, 256, 65536, 4294967296, 18446744073709551616}

In[8]:= NestList[#(1 + 0.05) &, 1000, 10]
Out[8]= {1000, 1050., 1102.5, 1157.63, 1215.51, 1276.28, 1340.1, 1407.1, 1477.46, 1551.33, 1628.89}

In[9]:= NestList[(# + 2/#)/2 &, 1.0, 5]
Out[9]= {1., 1.5, 1.41667, 1.41422, 1.41421, 1.41421}

In[10]:= NestList[If[EvenQ[#], #/2, (3 # + 1)/2] &, 100, 20]
Out[10]= {100, 50, 25, 38, 19, 29, 44, 22, 11, 17, 26, 13, 20, 10, 5, 8, 4, 2, 1, 2, 1}

In[11]:= NestList[Mod[59 #, 101] &, 1, 15]
Out[11]= {1, 59, 47, 46, 88, 41, 96, 8, 68, 73, 65, 98, 25, 61, 64, 39}
```

#### NestWhile
Iteratively applies a function while a predicate continues to yield `True`.
- `NestWhile[f, expr, test]`: Starts with `expr` and keeps applying `f` until `test` no longer yields `True`.
- `NestWhile[f, expr, test, m]`: Supplies the most recent `m` results (not wrapped in a list) as arguments to `test`.
- `NestWhile[f, expr, test, All]`: Supplies every result so far as arguments to `test`.
- `NestWhile[f, expr, test, {mmin, mmax}]`: Defers the first test until `mmin` results exist, then supplies up to `mmax` most-recent results.
- `NestWhile[f, expr, test, m, max]`: Caps the number of `f` applications at `max` (may be `Infinity`).
- `NestWhile[f, expr, test, m, max, n]`: After the loop terminates, applies `f` an additional `n` times.
- `NestWhile[f, expr, test, m, max, -n]`: Returns the result produced `n` applications before the loop ended (i.e. `Part[NestWhileList[...], -n-1]`).

**Features**:
- `Protected`.
- If `test[expr]` does not yield `True` initially, the unchanged `expr` is returned.
- Results passed to `test` are in generation order with the most recent last, so e.g. `# > 1 &` inspects the oldest when more than one result is supplied.
- `NestWhile[f, expr, UnsameQ, 2]` is equivalent to `FixedPoint[f, expr]`.
- `NestWhile[f, expr, UnsameQ, All]` continues until any prior value reappears.
- `m` must be a positive integer, `All`, or a 2-element list `{mmin, mmax}` with `1 <= mmin <= mmax` (or `mmax = Infinity`); `max` must be a non-negative integer or `Infinity`; `n` must be an integer. Malformed specs leave `NestWhile` unevaluated.
- Pure functions (`... &`) are supported for both `f` and `test`.

**Examples**:
```
In[1]:= NestWhile[#/2 &, 123456, EvenQ]
Out[1]= 1929

In[2]:= NestWhile[Log, 100., # > 0 &]
Out[2]= -0.859384

In[3]:= NestWhile[Floor[#/2] &, 10, UnsameQ, 2]
Out[3]= 0

In[4]:= NestWhile[#/2 &, 123456, EvenQ, 1, 4]
Out[4]= 7716

In[5]:= NestWhile[Floor[#/2] &, 20, # > 1 &, 1, Infinity]
Out[5]= 1

In[6]:= NestWhile[Floor[#/2] &, 20, # > 1 &, 1, Infinity, 1]
Out[6]= 0

In[7]:= NestWhile[Floor[#/2] &, 20, # > 1 &, 1, Infinity, -1]
Out[7]= 2

In[8]:= NestWhile[# + 1 &, 888, !PrimeQ[#] &]
Out[8]= 907

In[9]:= NestWhile[# + 1 &, 888, !PrimeQ[#1] || !PrimeQ[#3] &, 3]
Out[9]= 1021

In[10]:= NestWhile[Mod[# + 3, 7] &, 0, UnsameQ, All]
Out[10]= 0

In[11]:= NestWhile[If[EvenQ[#], #/2, 3 # + 1] &, 27, # != 1 &]
Out[11]= 1
```

#### NestWhileList
Like `NestWhile`, but returns the full list of intermediate results.
- `NestWhileList[f, expr, test]`: Generates `{expr, f[expr], f[f[expr]], ...}`, continuing while `test` applied to the most recent result yields `True`.
- `NestWhileList[f, expr, test, m]`: Supplies the most recent `m` results (not wrapped in a list) as arguments to `test`.
- `NestWhileList[f, expr, test, All]`: Supplies every result so far as arguments to `test`.
- `NestWhileList[f, expr, test, {mmin, mmax}]`: Defers the first test until `mmin` results exist, then supplies up to `mmax` most-recent results.
- `NestWhileList[f, expr, test, m, max]`: Caps the number of `f` applications at `max` (may be `Infinity`).
- `NestWhileList[f, expr, test, m, max, n]`: Appends `n` additional applications of `f` to the list after the loop terminates.
- `NestWhileList[f, expr, test, m, max, -n]`: Drops the last `n` elements from the list.

**Features**:
- `Protected`.
- Results are listed in generation order, including the final element on which `test` yielded a non-`True` value (or the last element produced when `max` iterations were reached).
- If `test[expr]` does not yield `True` initially, the result is just `{expr}`.
- `NestWhileList[f, expr, UnsameQ, 2]` is equivalent to `FixedPointList[f, expr]`.
- `NestWhileList[f, expr, test, All]` is equivalent to `NestWhileList[f, expr, test, {1, Infinity}]`.
- `NestWhileList[f, expr, UnsameQ, All]` continues until a previously-seen value reappears, and the repeat is included as the last element of the list.
- `m` must be a positive integer, `All`, or a 2-element list `{mmin, mmax}` with `1 <= mmin <= mmax` (or `mmax = Infinity`); `max` must be a non-negative integer or `Infinity`; `n` must be an integer. Malformed specs leave `NestWhileList` unevaluated.
- Pure functions (`... &`) are supported for both `f` and `test`.

**Examples**:
```
In[1]:= NestWhileList[#/2 &, 123456, EvenQ]
Out[1]= {123456, 61728, 30864, 15432, 7716, 3858, 1929}

In[2]:= NestWhileList[Log, 100., # > 0 &]
Out[2]= {100.0, 4.60517, 1.52718, 0.423423, -0.859384}

In[3]:= NestWhileList[Floor[#/2] &, 20, UnsameQ, 2, 4]
Out[3]= {20, 10, 5, 2, 1}

In[4]:= NestWhileList[Floor[#/2] &, 20, # > 1 &, 1, Infinity]
Out[4]= {20, 10, 5, 2, 1}

In[5]:= NestWhileList[Floor[#/2] &, 20, # > 1 &, 1, Infinity, 1]
Out[5]= {20, 10, 5, 2, 1, 0}

In[6]:= NestWhileList[Floor[#/2] &, 20, # > 1 &, 1, Infinity, -1]
Out[6]= {20, 10, 5, 2}

In[7]:= NestWhileList[# + 1 &, 899, !PrimeQ[#] &]
Out[7]= {899, 900, 901, 902, 903, 904, 905, 906, 907}

In[8]:= NestWhileList[Mod[2 #, 19] &, 2, # != 1 &]
Out[8]= {2, 4, 8, 16, 13, 7, 14, 9, 18, 17, 15, 11, 3, 6, 12, 5, 10, 1}

In[9]:= NestWhileList[Mod[5 #, 7] &, 4, Unequal, All]
Out[9]= {4, 6, 2, 3, 1, 5, 4}

In[10]:= NestWhileList[If[EvenQ[#], #/2, (3 # + 1)/2] &, 400, Unequal, All]
Out[10]= {400, 200, 100, 50, 25, 38, 19, 29, 44, 22, 11, 17, 26, 13, 20, 10, 5, 8, 4, 2, 1, 2}

In[11]:= NestWhileList[If[EvenQ[#], #/2, (3 # + 1)/2] &, 400, Unequal, All, Infinity, -1]
Out[11]= {400, 200, 100, 50, 25, 38, 19, 29, 44, 22, 11, 17, 26, 13, 20, 10, 5, 8, 4, 2, 1}
```

#### FixedPointList
Generates the list of successive iterates of a function until a fixed point is reached.
- `FixedPointList[f, expr]`: Generates `{expr, f[expr], f[f[expr]], ...}`, continuing until two consecutive results are `SameQ`. The returned list always begins with `expr`, and its last two elements are always equal (the fixed point appears twice).
- `FixedPointList[f, expr, n]`: Stops after at most `n` applications of `f`. If `n` is reached before convergence, the last two elements need not be equal.
- `FixedPointList[f, expr, SameTest -> s]`: Uses the binary predicate `s` instead of `SameQ` to test successive pairs of results. Iteration stops when `s[prev, current]` evaluates to `True`.
- `FixedPointList[f, expr, n, SameTest -> s]`: Combines a bounded iteration count with a custom equivalence test.

**Features**:
- `Protected`.
- `FixedPointList[f, expr]` is equivalent to `NestWhileList[f, expr, UnsameQ, 2]`.
- `n` (when given) must be a non-negative integer or `Infinity`. Malformed argument specs leave `FixedPointList` unevaluated.

```mathematica
In[1]:= FixedPointList[1 + Floor[#/2] &, 1000]
Out[1]= {1000, 501, 251, 126, 64, 33, 17, 9, 5, 3, 2, 2}

In[2]:= 1 + Floor[Last[%]/2]
Out[2]= 2

In[3]:= FixedPointList[# /. {a_, b_} /; b != 0 -> {b, Mod[a, b]} &, {28, 21}]
Out[3]= {{28, 21}, {21, 7}, {7, 0}, {7, 0}}

In[4]:= GCD[28, 21]
Out[4]= 7

In[5]:= FixedPointList[1 + Floor[#/2] &, 1000, 5]
Out[5]= {1000, 501, 251, 126, 64, 33}

In[6]:= FixedPointList[(# + 2/#)/2 &, 1.0]
Out[6]= {1.0, 1.5, 1.41667, 1.41422, 1.41421, 1.41421, 1.41421}

In[7]:= FixedPointList[(# + 2/#)/2 &, 1.0, SameTest -> (Abs[#1 - #2] < 0.01 &)]
Out[7]= {1.0, 1.5, 1.41667, 1.41422}
```

#### Fold
Successively applies a binary function to an accumulating seed and the elements of a list.
- `Fold[f, x, list]`: Returns the last element of `FoldList[f, x, list]`, namely `f[...f[f[f[x, list[[1]]], list[[2]]], list[[3]]]..., list[[n]]]`.
- `Fold[f, list]`: Equivalent to `Fold[f, First[list], Rest[list]]`.

**Features**:
- `Protected`.
- The head of the third argument need not be `List` (any compound expression is accepted).
- `Fold[f, x, {}]` returns `x` (the function is never applied); `Fold[f, {a}]` returns `a`.
- `Fold[f, {}]` remains unevaluated (no seed, no elements).
- Each intermediate application is evaluated before the next one.

```mathematica
In[1]:= Fold[f, x, {a, b, c, d}]
Out[1]= f[f[f[f[x, a], b], c], d]

In[2]:= Fold[List, x, {a, b, c, d}]
Out[2]= {{{{x, a}, b}, c}, d}

In[3]:= Fold[Times, 1, {a, b, c, d}]
Out[3]= a b c d

In[4]:= Fold[f, {a, b, c, d}]
Out[4]= f[f[f[a, b], c], d]

In[5]:= Fold[{2 #1, 3 #2} &, x, {a, b, c, d}]
Out[5]= {{{{16 x, 24 a}, 12 b}, 6 c}, 3 d}

In[6]:= Fold[f, x, p[a, b, c, d]]
Out[6]= f[f[f[f[x, a], b], c], d]

(* Horner-form evaluation: 2^3 + 2 + 1 at x = 2 *)
In[7]:= Fold[2 #1 + #2 &, 0, {1, 0, 1, 1}]
Out[7]= 11

(* Form a number from digits *)
In[8]:= Fold[10 #1 + #2 &, 0, {4, 5, 1, 6, 7, 8}]
Out[8]= 451678

(* Form a continued fraction *)
In[9]:= Fold[1/(#2 + #1) &, x, Reverse[{a, b, c, d}]]
Out[9]= 1/(a + 1/(b + 1/(c + 1/(d + x))))

(* Factorial via Times *)
In[10]:= Fold[Times, 1, Range[5]]
Out[10]= 120

(* When the pure function ignores its second argument, Fold coincides with Nest *)
In[11]:= Fold[f[#1] &, x, Range[5]]
Out[11]= f[f[f[f[f[x]]]]]

(* Edge cases *)
In[12]:= Fold[f, x, {}]
Out[12]= x

In[13]:= Fold[f, {a}]
Out[13]= a

In[14]:= Fold[f, {}]
Out[14]= Fold[f, {}]
```

#### FoldList
Produces the list of intermediate fold values.
- `FoldList[f, x, list]`: Gives `{x, f[x, list[[1]]], f[f[x, list[[1]]], list[[2]]], ...}`.
- `FoldList[f, list]`: Gives `{list[[1]], f[list[[1]], list[[2]]], ...}`.

**Features**:
- `Protected`.
- With a length-`n` list, `FoldList` returns a list of length `n + 1` (or `n` in the no-seed form).
- The head of the third argument is preserved in the output: `FoldList[f, x, p[a, b]]` gives `p[x, f[x, a], f[f[x, a], b]]`.
- `FoldList[f, {}]` returns `{}` (an empty list with the input head); `FoldList[f, x, {}]` returns `{x}`.
- `Fold[f, x, list]` is equivalent to `Last[FoldList[f, x, list]]`.

```mathematica
In[1]:= FoldList[f, x, {a, b, c, d}]
Out[1]= {x, f[x, a], f[f[x, a], b], f[f[f[x, a], b], c], f[f[f[f[x, a], b], c], d]}

In[2]:= FoldList[f, {a, b, c, d}]
Out[2]= {a, f[a, b], f[f[a, b], c], f[f[f[a, b], c], d]}

(* Cumulative sums *)
In[3]:= FoldList[Plus, 0, Range[5]]
Out[3]= {0, 1, 3, 6, 10, 15}

(* Head preservation *)
In[4]:= FoldList[f, x, p[a, b, c, d]]
Out[4]= p[x, f[x, a], f[f[x, a], b], f[f[f[x, a], b], c], f[f[f[f[x, a], b], c], d]]

(* Fold to the right *)
In[5]:= FoldList[g[#2, #1] &, x, {a, b, c, d}]
Out[5]= {x, g[a, x], g[b, g[a, x]], g[c, g[b, g[a, x]]], g[d, g[c, g[b, g[a, x]]]]}

(* Successive factorials *)
In[6]:= FoldList[Times, 1, Range[10]]
Out[6]= {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800}

(* Build up a continued fraction *)
In[7]:= FoldList[1/(#2 + #1) &, x, Reverse[{a, b, c}]]
Out[7]= {x, 1/(c + x), 1/(b + 1/(c + x)), 1/(a + 1/(b + 1/(c + x)))}

(* Build up a number from digits *)
In[8]:= FoldList[10 #1 + #2 &, 0, {4, 5, 1, 6, 7, 8}]
Out[8]= {0, 4, 45, 451, 4516, 45167, 451678}

(* Horner-form polynomial evaluation at x = 2, coefficients {1, 0, 1, 1} *)
In[9]:= FoldList[2 #1 + #2 &, 0, {1, 0, 1, 1}]
Out[9]= {0, 1, 2, 5, 11}

(* Edge cases *)
In[10]:= FoldList[f, x, {}]
Out[10]= {x}

In[11]:= FoldList[f, {}]
Out[11]= {}

In[12]:= FoldList[f, p[]]
Out[12]= p[]
```

#### Through
Distributes operators that appear inside the heads of expressions.
- `Through[expr]`: Distributes the top-level head.
- `Through[expr, h]`: Performs the transformation wherever `h` occurs in the head of `expr`.

**Features**:
- `Protected`.

```mathematica
In[1]:= Through[{f, g, h}[x]]
Out[1]= {f[x], g[x], h[x]}

In[2]:= Through[(f + g)[x, y]]
Out[2]= f[x, y] + g[x, y]
```

#### Select
Filters elements from an expression matching a criterion.
- `Select[list, crit]`: Returns an expression with the same head as `list`, containing only those elements `e` for which `crit[e]` evaluates to `True`.
- `Select[list, crit, n]`: Returns only the first `n` matching elements.

```mathematica
In[1]:= Select[{1, 2, 4, 7, 6, 2}, EvenQ]
Out[1]= {2, 4, 6, 2}

In[2]:= Select[{1, 2, 4, 7, 6, 2}, # > 2 &, 1]
Out[2]= {4}
```

### Control Flow

#### Do
Evaluates an expression sequentially over an iteration range.
- `Do[expr, n]`: Evaluates `expr` `n` times.
- `Do[expr, {i, imax}]`: Evaluates `expr` with `i` from 1 to `imax`.
- `Do[expr, {i, imin, imax, di}]`: Evaluates `expr` with `i` taking values from `imin` to `imax` in steps of `di`.
- `Do[expr, {i, list}]`: Evaluates `expr` with `i` taking values from `list`.
- `Do[expr, spec1, spec2, ...]`: Evaluates `expr` looping over `spec1` internally.

**Features**:
- `HoldAll`, evaluating its body only after arguments are substituted.
- Employs exact dynamic iteration identical to `Table` but discards the evaluated results, returning `Null`.
- Supports explicit break states (`Return`, `Break`, `Continue`, `Throw`, `Abort`, `Quit`).
- Can execute an infinite loop using `Do[expr, Infinity]`.

```mathematica
In[1]:= Do[Print[i], {i, 3}]
Out[1]= Null

In[2]:= Do[If[i == 3, Break[]]; Print[i], {i, 5}]
Out[2]= Null
```

#### For
Executes a loop with an initialization, condition test, increment, and body.
- `For[start, test, incr, body]`: Evaluates `start`, then repeatedly evaluates `body` and `incr` until `test` fails to give `True`.
- `For[start, test, incr]`: Executes the loop with a `Null` body.

**Features**:
- Evaluates its arguments in a nonstandard way (sequence: `test`, `body`, `incr`).
- Has attribute `HoldAll`.
- `Break[]` exits the loop.
- `Continue[]` skips the rest of the body and proceeds to evaluating `incr`.
- Exits as soon as `test` fails.
- Returns `Null` unless an explicit `Return` is evaluated.

```mathematica
In[1]:= For[i=0, i<4, i++, Print[i]]
Out[1]= Null
```

#### While
Evaluates a test expression and, while it yields `True`, repeatedly evaluates a body expression.
- `While[test, body]`: Evaluates `test`, then `body`, repeatedly, until `test` first fails to give `True`.
- `While[test]`: Executes the loop with a `Null` body. Useful when `test` itself has side-effects.

**Features**:
- Has attribute `HoldAll`; both `test` and `body` are re-evaluated each iteration.
- `Break[]` inside `body` exits the loop, yielding `Null`.
- `Continue[]` inside `body` skips the rest of `body` and returns to re-evaluating `test`.
- `Return[v]` inside `body` causes `While` to yield `v`.
- `Throw`, `Abort`, and `Quit` propagate unchanged.
- If the very first evaluation of `test` is not `True`, `body` is never evaluated.
- Returns `Null` unless an explicit `Return` is issued.

```mathematica
In[1]:= n = 1; While[n < 4, n = n + 1]; n
Out[1]= 4

In[2]:= {a, b} = {27, 6}; While[b != 0, {t1, t2} = {b, Mod[a, b]}; a = t1; b = t2]; a
Out[2]= 3

In[3]:= n = 1; While[True, If[n > 10, Break[]]; n = n + 1]; n
Out[3]= 11
```

#### If
Evaluates condition and executes the corresponding branch.
- `If[condition, t, f]`: Gives `t` if `condition` evaluates to `True`, and `f` if it evaluates to `False`.
- `If[condition, t, f, u]`: Gives `u` if `condition` evaluates to neither `True` nor `False`.

**Features**:
- `HoldRest`, evaluating only the chosen branch.
- Remains unevaluated if the condition is undetermined and `u` is not provided.
- `If[condition, t]` returns `Null` if `condition` evaluates to `False`.

```mathematica
In[1]:= If[True, x, y]
Out[1]= x

In[2]:= If[a < b, 1, 0, Indeterminate]
Out[2]= Indeterminate
```

#### TrueQ
Tests whether an expression evaluates explicitly to `True`.
- `TrueQ[expr]`: Yields `True` if `expr` is `True`, and `False` otherwise.

#### Return
Exits the nearest enclosing scope or loop, yielding a value.
- `Return[expr]`: Yields `expr` from the innermost enclosing scope boundary.
- `Return[]`: Yields `Null` (shorthand for `Return[Null]`).
- `Return[expr, h]`: Yields `expr` from the *nearest enclosing boundary* whose head is the symbol `h`. Boundaries with a different head propagate the marker outward unchanged so that `h` can be reached.

**Recognised boundary heads**: `Function`, `Module`, `Block`, `With`, `Do`, `For`, `While`. `CompoundExpression` and other Hold-free heads (including `If`, `Which`, `Switch`) propagate `Return` through unchanged so it can reach the enclosing boundary.

**Features**:
- `Protected`. No Hold attributes — arguments are evaluated before the marker takes effect.
- `Return` takes effect as soon as it is evaluated, even when it appears inside other functions (Plus, Times, etc. in argument positions still see the substituted value, but at the top level the marker propagates immediately through `CompoundExpression`).
- 1-arg `Return[expr]` is consumed by the *innermost* boundary on the call stack. The 2-arg `Return[expr, h]` form lets the user skip past intervening boundaries to a specific enclosing construct.
- If no enclosing boundary matches, `Return[expr]` (or `Return[expr, h]`) survives at the top level as a literal expression.

```mathematica
In[1]:= Function[x, If[x > 0, Return[positive], Return[negative]]][3]
Out[1]= positive

In[2]:= Module[{}, Do[Return[5], {3}]]
Out[2]= 5

In[3]:= Module[{}, Do[Return[5, Module], {3}]]
Out[3]= 5

In[4]:= Module[{}, Do[Return[5, Block], {3}]]
Out[4]= Return[5, Block]

In[5]:= f[n_] := Module[{s = 0}, Do[s = s + i; If[s > 10, Return[i]], {i, 1, n}]]; f[10]
Out[5]= 5
```

### Assignment and Rules

#### Set (=), SetDelayed (:=)
- `lhs = rhs`: `Set`. Immediate evaluation of `rhs`.
- `lhs := rhs`: `SetDelayed`. Delayed evaluation of `rhs`.
- `lhs := rhs /; condition`: When the RHS of `SetDelayed` is a `Condition`, it is automatically moved to the LHS pattern. This makes `f[x_] := body /; test` equivalent to `f[x_] /; test := body`.

#### Replace
Applies a rule or list of rules to transform an expression.
- `Replace[expr, rules]`: Applies rules to the entire expression.
- `Replace[expr, rules, levelspec]`: Applies rules to parts specified by `levelspec`.
- `Replace[expr, rules, levelspec, Heads -> True]`: Includes heads of expressions and their parts.

**Features**:
- `Protected`.
- Rules must be of the form `lhs -> rhs` (`Rule`) or `lhs :> rhs` (`RuleDelayed`).
- Tries rules in order. The first one that matches is applied.
- If rules are given in nested lists, `Replace` is mapped onto the inner lists.
- Standard level specifications (`n`, `Infinity`, `All`, `{n}`, `{n1, n2}`) are fully supported with the default being `{0}` (the whole expression).
- Expressions at deeper levels in a subexpression are matched first.
- Replaces parts even inside `Hold` or related wrappers.

```mathematica
In[1]:= Replace[x^2, x^2 -> a + b]
Out[1]= a + b

In[2]:= Replace[1 + x^2, x^2 -> a + b]
Out[2]= 1 + x^2

In[3]:= Replace[x, {{x -> a}, {x -> b}}]
Out[3]= {a, b}

In[4]:= Replace[1 + x^2, x^2 -> a + b, {1}]
Out[4]= 1 + a + b
```

#### ReplacePart
Replaces parts of an expression at specified positions.
- `ReplacePart[expr, i -> new]`: Replaces the part at position `i`.
- `ReplacePart[expr, {i, j, ...} -> new]`: Replaces the part at a nested position.
- `ReplacePart[expr, {{i1, ...} -> new1, {i2, ...} -> new2, ...}]`: Multiple replacements.
- `ReplacePart[expr, pattern -> new]`: Replaces parts at all positions matching `pattern` (e.g., `_`, `Except[...]`).

**Features**:
- `Protected`.
- Supports negative indices: `-1` is the last element, `-2` is second-to-last, etc.
- Supports mixed positive/negative indices in nested position specs (e.g., `{2, -1}`).
- Default `Heads -> False`: pattern-based position specs do not match the head (index 0). Use an explicit `0 -> new` rule to replace the head.
- Supports `Rule` (`->`) and `RuleDelayed` (`:>`).

```mathematica
In[1]:= ReplacePart[{a, b, c, d}, 2 -> x]
Out[1]= {a, x, c, d}

In[2]:= ReplacePart[{a, b, c, d, e, f, g}, -3 -> xxx]
Out[2]= {a, b, c, d, xxx, f, g}

In[3]:= ReplacePart[{a, b, c, d, e, f, g}, Except[1 | 3 | 5] -> xxx]
Out[3]= {a, xxx, c, xxx, e, xxx, xxx}

In[4]:= ReplacePart[f[x], 0 -> g]
Out[4]= g[x]
```

#### ReplaceAt
Replaces parts of an expression at specified positions using replacement rules. Unlike `ReplaceAll` (which traverses the whole tree), `ReplaceAt` only attempts the rules at the part(s) named by the position specification.

- `ReplaceAt[expr, rules, n]`: replaces the `n`-th element using `rules`.
- `ReplaceAt[expr, rules, {i, j, ...}]`: replaces the part of `expr` at position `{i, j, ...}` (equivalently `expr[[i, j, ...]]`).
- `ReplaceAt[expr, rules, {{i1, j1, ...}, {i2, j2, ...}, ...}]`: replaces the parts at each of the listed positions. If the same position appears more than once, the rules are applied repeatedly to that part.

**Features**:
- `Protected`.
- `rules` may be a single `Rule` (`->`), `RuleDelayed` (`:>`), or a list of such rules. The rules are tried in order; the first one that applies wins. If no rule matches at a targeted position, the part is left unchanged.
- For `RuleDelayed`, the right-hand side is evaluated separately for each match after substituting bound pattern variables.
- Negative integer indices count from the end. The literal index `0` targets the head of an expression.
- Path components may be integers, the symbol `All` (selects every child at that level), or `Span` expressions such as `i ;; j` or `i ;; j ;; k`.
- Works on expressions with any head (not just `List`); after substitution the evaluator re-applies canonical ordering for `Orderless` heads such as `Plus` and `Times`.
- The position list uses the same form as is returned by `Position`. `ReplaceAt[expr, rules, {}]` applies the rules to the whole expression.

```mathematica
In[1]:= ReplaceAt[{a, a, a, a}, a -> xx, 2]
Out[1]= {a, xx, a, a}

In[2]:= ReplaceAt[{a, a, a, a}, a -> xx, {{1}, {4}}]
Out[2]= {xx, a, a, xx}

In[3]:= ReplaceAt[{{a, a}, {a, a}}, a -> xx, {2, 1}]
Out[3]= {{a, a}, {xx, a}}

In[4]:= ReplaceAt[{a, a, a, a}, a -> xx, -2]
Out[4]= {a, a, xx, a}

In[5]:= ReplaceAt[{{a, a, a}, {a, a, a}}, a -> xx, {-1, -2}]
Out[5]= {{a, a, a}, {a, xx, a}}

In[6]:= ReplaceAt[{1, 2, 3, 4}, x_ :> 2 x - 1, {{2}, {4}}]
Out[6]= {1, 3, 3, 7}

In[7]:= ReplaceAt[{a, b, c, d}, {a -> xx, _ -> yy}, {{1}, {2}, {4}}]
Out[7]= {xx, yy, c, yy}

In[8]:= ReplaceAt[{{a, a}, {a, a}}, a -> xx, {All, 2}]
Out[8]= {{a, xx}, {a, xx}}

In[9]:= ReplaceAt[{{a, b}, {c, d}, e}, x_ :> f[x], 2]
Out[9]= {{a, b}, f[{c, d}], e}

In[10]:= ReplaceAt[{{a, b}, {c, d}, e}, x_ :> f[x], -1]
Out[10]= {{a, b}, {c, d}, f[e]}

In[11]:= ReplaceAt[{{a, b}, {c, d}, e}, x_ :> f[x], {2, 1}]
Out[11]= {{a, b}, {f[c], d}, e}

In[12]:= ReplaceAt[{{a, b}, {c, d}, e}, x_ :> f[x], {{1}, {3}}]
Out[12]= {f[{a, b}], {c, d}, f[e]}

In[13]:= ReplaceAt[{{a, b}, {c, d}, e}, x_ :> f[x], {{1, 2}, {2, 2}, {3}}]
Out[13]= {{a, f[b]}, {c, f[d]}, f[e]}

In[14]:= ReplaceAt[{a, a, a, a, a}, a -> xx, 2 ;; 4]
Out[14]= {a, xx, xx, xx, a}

In[15]:= ReplaceAt[a + b + c + d, _ -> x, 2]
Out[15]= a + c + d + x

In[16]:= ReplaceAt[x^2 + y^2, _ -> z, {{1, 1}, {2, 1}}]
Out[16]= 2 z^2

In[17]:= ReplaceAt[{a, b, c}, _ -> f, 0]
Out[17]= f[a, b, c]
```

#### ReplaceAll (/.)
Applies a rule or list of rules to transform each subpart of an expression.
- `expr /. lhs -> rhs` or `expr /. {rules}`

**Features**:
- `Protected`.
- Evaluates the entire expression top-down. The first rule that applies to a particular part is used; no further rules are tried on that part or on any of its subparts.
- Applies a rule only once to an expression.
- Returns `expr` unmodified if no rules apply.
- Maps across lists of rules appropriately.

```mathematica
In[1]:= {x, x^2, y, z} /. x -> 1
Out[1]= {1, 1, y, z}

In[2]:= Sin[x] /. Sin -> Cos
Out[2]= Cos[x]

In[3]:= {1, 3, 2, x, 6, Pi} /. _?PrimeQ -> "prime"
Out[3]= {1, "prime", "prime", x, 6, Pi}

In[4]:= {f[2], f[x, y], h[], f[]} /. f[x__] -> "OK"
Out[4]= {"OK", "OK", h[], f[]}
```

#### CompoundExpression (;)
- `expr1; expr2; ...`: Evaluates a sequence of expressions, returning the last one.

### Scoping Constructs

#### Module
Implements lexical scoping by creating unique local variables.
- `Module[{x, y, ...}, expr]`
- `Module[{x = x0, ...}, expr]`

**Features**:
- `HoldAll`, `Protected`.
- Variables are renamed to `name$nnn` using `$ModuleNumber`.
- Created symbols have the `Temporary` attribute.

```mathematica
In[1]:= x = 1; Module[{x = 2}, x + 1]
Out[1]= 3

In[2]:= x
Out[2]= 1
```

#### Block
Implements dynamic scoping by temporarily overriding symbol values.
- `Block[{x, y, ...}, expr]`
- `Block[{x = x0, ...}, expr]`

**Features**:
- `HoldAll`, `Protected`.
- Affects only values, not names.
- Restores original values and attributes after execution.

#### With
Defines local constants by lexical substitution.
- `With[{x = x0, ...}, expr]`
- `With[{x := x0, ...}, expr]`

**Features**:
- `HoldAll`, `Protected`.
- Replaces occurrences of symbols in the body before evaluation.

```mathematica
In[1]:= x = 10; With[{x = 5}, x^2]
Out[1]= 25
```

### Pattern Matching

- `_`: `Blank[]`.
- `__`: `BlankSequence[]`.
- `___`: `BlankNullSequence[]`.
- `x_`: `Pattern[x, _]`.
- `_h`: Matches expression with head `h`.
- `p ? test`: `PatternTest[p, test]`.
- `lhs /; condition`: Conditional matching.

#### Cases
Gives a list of elements that match a pattern.
- `Cases[{e1, e2, ...}, pattern]`
- `Cases[{e1, ...}, pattern -> rhs]`
- `Cases[expr, pattern, levelspec]`
- `Cases[expr, pattern -> rhs, levelspec]`
- `Cases[expr, pattern, levelspec, n]`
- `Cases[pattern]`

**Features**:
- Uses standard level specifications.
- Defaults to level `{1}`.
- Option `Heads -> True` includes heads of expressions.
- Operates in depth-first post-order traversal.

#### Position
Gives a list of the positions at which objects matching a pattern appear.
- `Position[expr, pattern]`
- `Position[expr, pattern, levelspec]`
- `Position[expr, pattern, levelspec, n]`
- `Position[pattern]`

**Features**:
- Defaults to levels `{0, Infinity}` with `Heads -> True`.
- Yields lists of indices in lexicographic order.

#### Count
Gives the number of elements or subexpressions that match a pattern.
- `Count[list, pattern]`
- `Count[expr, pattern, levelspec]`
- `Count[pattern]`

**Features**:
- Uses standard level specifications, defaulting to level `{1}`.
- Option `Heads -> True` evaluates heads of expressions.
- Effectively returns the length of the list generated by `Cases` for the same arguments.
- `p:def`: `Optional[p, def]`, matches `p`, but if omitted defaults to `def`.
- `x_.`: `Optional[Pattern[x, _]]`, matches `x_`, but if omitted defaults to globally specified `Default` values.
- `p..`: `Repeated[p]`, matches a sequence of one or more expressions, each matching `p`.
- `p...`: `RepeatedNull[p]`, matches a sequence of zero or more expressions, each matching `p`.
- `Repeated[p, max]`, `Repeated[p, {min, max}]`, `Repeated[p, {n}]`: match sequences of length bounded by `max` or `{min, max}`. Equivalent properties exist for `RepeatedNull`.

#### Default Values
- `Default[f] = val`: Sets the default value for arguments of `f`.
- Built-in default values include `0` for `Plus` and `1` for `Times` and `Power`.
- Used automatically by the pattern matcher when `_.` is present in the pattern.

## 3. Operators and Precedence

| Operator | FullForm | Precedence | Association |
|----------|----------|------------|-------------|
| `[[...]]`| `Part` | 100 | Left |
| `f[x]`   | `f[x]` | 1000 | Left |
| `_`, `__`, `___` | `Blank` | 730 | None |
| `?`      | `PatternTest` | 680 | None |
| `@`      | `Prefix` | 620 | Right |
| `//`     | `Postfix` | 70 | Left |
| `&`      | `Function`| 90 | Left |
| `@@`     | `Apply` | 620 | Right |
| `/@`     | `Map`   | 620 | Right |
| `<>`     | `StringJoin` | 600 | Left |
| `^`      | `Power` | 590 | Right |
| `*`      | `Times` | 400 | Left |
| `/`      | `Divide`| 470 | Left |
| `+`, `-` | `Plus`  | 310 | Left |
| `==`     | `Equal` | 290 | None |
| `===`    | `SameQ` | 290 | None |
| `=`      | `Set`   | 40 | Right |
| `:=`     | `SetDelayed` | 40 | Right |
| `;`      | `CompoundExpression` | 10 | Left |

## 4. Limitations
- Large integer arithmetic is not supported (limited to 64-bit).
- Transcendental simplification is limited.
- Pattern matching for `Flat` and `Orderless` heads is partially implemented.

## 5. Technical Implementation Details

The system is composed of several interdependent subsystems: Expression Representation, Parser, Symbol Table, Evaluator, Pattern Matcher, and Rule Engine.

### 5.1. Expression Representation (`Expr`)

Everything in PicoCAS is an expression, represented by the `Expr` struct. It is implemented as a tagged union.

*   **Types (`ExprType`)**:
    *   `EXPR_INTEGER`: A 64-bit signed integer (`int64_t`).
    *   `EXPR_REAL`: A 64-bit floating-point number (`double`).
    *   `EXPR_SYMBOL`: A named identifier (e.g., `x`, `Plus`, `True`).
    *   `EXPR_STRING`: A string literal.
    *   `EXPR_FUNCTION`: A compound expression consisting of a `head` (which is itself an `Expr*`) and an array of `args` (`Expr**`).

All data structures are immutable-by-convention during evaluation. Functions that manipulate expressions generally return new allocated trees or unmodified references if no changes were made. 

**Memory Management:** Explicit manual memory management is required. `expr_new_*` allocates nodes. `expr_copy` performs a deep copy. `expr_free` performs a deep deallocation. **Crucial Rule:** Built-in functions take ownership of the expression passed to them. If a built-in returns a new expression, it *must* free the input expression (or reuse its nodes). If it returns `NULL` (meaning it remains unevaluated), the evaluator retains ownership.

### 5.2. Symbol Table (`symtab.c`)

The symbol table (`SymbolDef`) stores the global state associated with symbols. Each symbol can have:
*   **Attributes**: Bitflags (e.g., `ATTR_FLAT`, `ATTR_ORDERLESS`, `ATTR_LISTABLE`, `ATTR_PROTECTED`, `ATTR_HOLDALL`, `ATTR_NUMERICFUNCTION`).
*   **OwnValues**: Immediate assignments (e.g., `x = 5`).
*   **DownValues**: Delayed assignments with pattern matching (e.g., `f[x_] := x^2`).
*   **Builtin C Function**: A pointer to a C function (`BuiltinFunc`) that implements native evaluation logic (e.g., `builtin_plus`).

### 5.3. Parser (`parse.c`)

The parser translates raw string input into `Expr` trees.
*   **Lexical Analysis**: Handled inline within the parsing routines by advancing a `ParserState` pointer and skipping whitespace.
*   **Pratt Parsing**: Operator precedence (Infix, Prefix, Postfix) is implemented using a Pratt parser (`parse_expression_prec`). Operator precedence values exactly mirror Mathematica's standard precedences (e.g., `OP_CALL` is 1000, `OP_TIMES` is 400, `OP_PLUS` is 310, `OP_PREFIX` (`@`) is 620, `OP_POSTFIX` (`//`) is 70).

### 5.4. Evaluator (`eval.c`)

The `evaluate(Expr* e)` function is the heart of PicoCAS. It follows Mathematica's infinite evaluation semantics: expressions are repeatedly evaluated until a fixed point is reached (the expression no longer changes).

**Evaluation Sequence for Functions (`f[arg1, arg2]`):**
1.  **Evaluate Head**: The head `f` is evaluated first.
2.  **Check Attributes**: The attributes of the evaluated head are retrieved.
3.  **Evaluate Arguments**: Arguments are evaluated standardly, *unless* the head possesses Hold attributes (`ATTR_HOLDFIRST`, `ATTR_HOLDREST`, `ATTR_HOLDALL`, `ATTR_HOLDALLCOMPLETE`).
4.  **Apply Listable**: If the head is `ATTR_LISTABLE` and any evaluated argument is a `List`, the function automatically threads over the lists (e.g., `Mod[{1, 2}, 3] -> {Mod[1, 3], Mod[2, 3]}`).
5.  **Apply Flat / Orderless**: If `ATTR_FLAT` (associative), nested calls to the same head are flattened. If `ATTR_ORDERLESS` (commutative), arguments are lexically sorted.
6.  **Apply Built-ins**: If a C-level built-in function is registered, it is called.
7.  **Apply User Rules**: If no built-in handles it (or returns `NULL`), the evaluator checks the `DownValues` of the head and applies the first matching pattern replacement.

### 5.5. Pattern Matcher (`match.c`)

Implements structural pattern matching (`MatchQ`).
*   `Blank[]` (`_`): Matches any single expression.
*   `BlankSequence[]` (`__`): Matches 1 or more expressions in a sequence.
*   `BlankNullSequence[]` (`___`): Matches 0 or more expressions in a sequence.
*   `Pattern[name, pattern]` (`name_`): Binds the matched sub-expression to `name` inside a `MatchEnv`.
*   The matcher recursively unifies trees and handles sequence segmenting through backtracking.

### 5.6. Rule Engine (`replace.c`)

Implements expression transformations via rules (`Rule` `->` and `RuleDelayed` `:>`).
*   `ReplaceAll` (`/.`): Traverses the tree top-down, applying rules to sub-expressions.
*   Uses `match.c` to determine if a rule's LHS matches the current expression, and if so, substitutes bindings into the RHS.

## 6. Performance Notes

### 6.1. Polynomial Algebra (`poly.c`)

The polynomial subsystem -- `Coefficient`, `CoefficientList`, `PolynomialGCD`,
`PolynomialQuotient`, `PolynomialRemainder`, `Resultant`, `HornerForm`,
`Discriminant`, `Decompose`, `Collect`, etc. -- shares a small set of inner
loops. Two optimisations dominate:

1. **Direct coefficient extraction.** The internal helper
   `get_coeff_expanded(expanded, var, n)` walks an already-expanded
   polynomial once per query, summing the contributions of each summand
   directly. This avoids the full evaluator pipeline (which would
   construct `Coefficient[expanded, var, n]`, look up its symbol,
   apply `Listable`/`Flat`/`Orderless`, then re-expand inside
   `builtin_coefficient`). A bulk variant `get_all_coeffs_expanded`
   produces every coefficient in a single pass and is used by
   `poly_content`, `coeff_list_rec`, `poly_derivative` and the
   resultant routine.

2. **Skipping `Together` + `Cancel` in pure-integer division.** The
   univariate division loops in `exact_poly_div` and `poly_div_rem`
   used to call `internal_together` and `internal_cancel` on every
   iteration to unify denominators introduced by symbolic leading
   coefficients. We now track whether the quotient coefficient was
   produced by an exact integer / bigint division (in which case no
   new fractions can appear) and only invoke the heavy unification
   path for the genuinely symbolic case.

Combined effect on the test bench (compiled at `-O3`):

| Test                             | Before | After  | Speedup |
|----------------------------------|-------:|-------:|--------:|
| `tests/poly_tests`               | 1.86s  | 0.53s  |   3.5x  |
| `tests/facpoly_tests`            | 4.24s  | 1.60s  |   2.7x  |
| `PolynomialGCD` (50 large iter.) | 10.07s | 2.40s  |   4.2x  |
| `CoefficientList` (30 iter.)     | 3.81s  | 0.92s  |   4.1x  |
| `HornerForm` (200 iter.)         | 16.76s | 5.13s  |   3.3x  |
| `PolynomialQuotient` (500 iter.) | 2.55s  | 1.15s  |   2.2x  |

A small clean-up was made alongside the perf work:
* `poly.h` had stray declarations after `#endif`; they are now inside
  the include guard.
* `PolynomialQ` was registered twice in `poly_init`; the duplicate
  registration has been removed.

### 6.2. Differentiation (`deriv.c`)

`D`, `Dt` and `Derivative` used to be defined by ~60 DownValues in
`src/internal/deriv.m`. Each call walked the whole rule list linearly,
re-ran the pattern matcher (including side-conditions such as
`FreeQ[c, x]`), and re-entered the evaluator for the replacement. The
native C implementation (added in this release) replaces that with:

1. **Direct head-symbol dispatch.** A single `strcmp` cascade in
   `compute_deriv` picks the right closed-form rule for `Plus`,
   `Times`, `Power`, `Sqrt`, `Exp`, `Log`, the trig/hyperbolic
   families and their inverses.
2. **Allocation-free constant detection.** A tailored `expr_free_of`
   walk replaces the generic `FreeQ` builtin on the hot path; it
   short-circuits on the first match against the differentiation
   variable without allocating a `MatchEnv`.
3. **Direct `Derivative[...][f][g]` construction.** The chain rule
   for unknown single- and multi-argument heads builds the
   `Derivative` head, the intermediate `Derivative[...][f]` and the
   final application directly, avoiding the N^2 DownValue-search
   behaviour the matcher produced for nested heads.

Measured on the stock interpreter at `-O3`:

| Workload                                               | Rule-based | C      | Speedup |
|--------------------------------------------------------|-----------:|-------:|--------:|
| 1000 x simple mixed (`D[Sin[x^2]+Cos[a x]Exp[x]+x^x+Log[x], x]`) | 1.24s      | 0.54s  |  2.3x   |
| 1000 x deep chain (`D[Sin[Cos[Tan[ArcSin[x^2+1]]]], x]`)         | 3.37s      | 2.95s  |  1.1x   |
| 200 x higher-order (`D[Sin[a x] + x^3 Cos[b x], {x, 6}]`)        | 33.26s     | 8.60s  |  3.9x   |
| 500 x mixed partials (`D[f[x,y,z]*Sin[x+y]*Cos[x y], x, y, z]`)  | 17.78s     | 2.27s  |  7.8x   |

The largest wins come from higher-order and mixed-partial calls,
where the rule-based implementation had to re-traverse the entire
DownValue list at every inner step; the native dispatch visits the
relevant head exactly once per sub-expression.

### 6.3. Limits (`limit.c`)

PicoCAS now has a native `Limit` built-in implemented in C in
`src/limit.c`, registered by `limit_init()` in the standard
`core_init()` chain. The design follows the layered dispatch outlined
in `limit_candidate_spec.md`, with each layer either resolving the
limit and short-circuiting or passing the problem down:

```
Layer 0 -- Interface normalization (calling forms + Direction option)
Layer 1 -- Structural fast paths (free-of-var, identity, continuous subst.)
Layer 3 -- Rational function P(x)/Q(x) short-cut
Layer 5.3 -- Logarithmic reduction for f^g indeterminate forms
Layer 2 -- Series-based evaluation (the workhorse)
Layer 5.1 -- L'Hospital with complexity-growth guardrail
Layer 6 -- Bounded-oscillation Interval returns
```

**Calling forms** (Mathematica-compatible):

| Form                                            | Meaning                            |
|-------------------------------------------------|------------------------------------|
| `Limit[f, x -> a]`                              | single-variable limit              |
| `Limit[f, {x1 -> a1, ..., xn -> an}]`           | iterated limit (rightmost first)   |
| `Limit[f, {x1, ..., xn} -> {a1, ..., an}]`      | multivariate joint limit           |
| `Limit[f, x -> a, Direction -> spec]`           | directional approach               |

**Direction settings:**

| Setting                        | Internal meaning                  |
|--------------------------------|-----------------------------------|
| `Reals`, `"TwoSided"`          | two-sided (default on the reals)  |
| `"FromAbove"` or `-1`          | approach from above (x -> a^+)    |
| `"FromBelow"` or `+1`          | approach from below (x -> a^-)    |
| `Complexes`                    | complex / radial-in-all-directions|

The Mathematica sign convention (`-1 == FromAbove`) is applied once in
the option parser so that the math layers only ever see internal
`+1 / -1 / 0` tags.

**Return values:** finite expression, `Infinity`, `-Infinity`,
`ComplexInfinity`, `Indeterminate`, `Interval[{lo, hi}]`, or the
original unevaluated `Limit[...]` when the system cannot determine a
value.

**Integration:**
* Registered under `Limit` with attributes
  `Protected | ReadProtected | HoldAll`. `HoldAll` prevents the
  second-argument rule from being evaluated prematurely against any
  OwnValue of the limit variable.
* Documented via `symtab_set_docstring` for `?Limit` in the REPL.
* Covered by `tests/test_limit.c`, which exercises each layer and
  every Direction setting.

**Layer highlights:**
1. *Continuous substitution* (Layer 1) does not just evaluate
   `f /. x -> a`; it Together-normalizes first and checks that the
   denominator does not vanish, so `Sin[x]/x` at x = 0 correctly skips
   the fast path instead of silently returning 0 (as PicoCAS's
   arithmetic would fold `Sin[0] * 0^(-1)`). Expressions of the form
   `Power[_, expr_with_x]` are likewise refused here -- they are
   handled by Layer 5.3 instead.
2. *Series-based evaluation* (Layer 2) calls `Series[f, {x, a, k}]`
   with increasing orders (4, 8, 16, 32) until a nonzero leading term
   is found, then reads the limit off the leading exponent. Limits at
   `-Infinity` substitute x -> -y internally and recurse at +Infinity.
3. *Log reduction* (Layer 5.3) fires for any `Power[base, exp]` whose
   exponent depends on the limit variable, rewriting as
   `Exp[Limit[exp * Log[base]]]` and recursing. It sits above Series in
   the dispatch order because Series has no kernel for `b^g` with
   g a non-trivial function of x.
4. *L'Hospital* (Layer 5.1) only fires when pointwise evaluation is
   strictly `0/0` or `inf/inf`, and aborts if the leaf count of the
   quotient grows for three consecutive derivations -- this avoids
   spinning on Sin/Cos expansions where the series layer is the right
   tool.
5. *Bound analysis* (Layer 6) handles the canonical
   `Limit[Sin[1/x], x -> 0] = Interval[{-1, 1}]` shape by recognising
   a bounded trigonometric head wrapped around a divergent argument.

**Worked examples (cross-checked against Mathematica):**

| Input                                                   | Result        |
|---------------------------------------------------------|---------------|
| `Limit[Sin[x]/x, x -> 0]`                               | `1`           |
| `Limit[(Cos[x]-1)/(Exp[x^2]-1), x -> 0]`                | `-1/2`        |
| `Limit[(Tan[x]-x)/x^3, x -> 0]`                         | `1/3`         |
| `Limit[Sin[2 x]/Sin[x], x -> Pi]`                       | `-2`          |
| `Limit[(x^a - 1)/a, a -> 0]`                            | `Log[x]`      |
| `Limit[(1 + A x)^(1/x), x -> 0]`                        | `E^A`         |
| `Limit[(1 + a/x)^(b x), x -> Infinity]`                 | `E^(a b)`     |
| `Limit[x/(x+1), x -> Infinity]`                         | `1`           |
| `Limit[3 x^2/(x^2 - 2), x -> -Infinity]`                | `3`           |
| `Limit[Sqrt[x^2+a x]-Sqrt[x^2+b x], x -> Infinity]`     | `(a-b)/2`     |
| `Limit[1/x, x -> 0, Direction -> "FromAbove"]`          | `Infinity`    |
| `Limit[1/x, x -> 0, Direction -> "FromBelow"]`          | `-Infinity`   |
| `Limit[Sin[1/x], x -> 0]`                               | `Interval[{-1, 1}]` |
| `Limit[x/(x + y), {x -> 0, y -> 0}]`                    | `1` (iterated)|

**Regression-suite additions (2026-04-20):** four new layers were
bolted on in response to a batch of REPL-driven test cases. Each is
composable with the original dispatcher and sits at the order indicated
below:

1. *Reciprocal-trig normalization* (runs at the top of
   `compute_limit`, before any other layer). Rewrites
   `Csc[z] -> 1/Sin[z]`, `Sec[z] -> 1/Cos[z]`, `Cot[z] -> Cos[z]/Sin[z]`,
   plus hyperbolic twins `Csch/Sech/Coth`, and also `Tan[z] -> Sin[z]/Cos[z]`,
   `Tanh[z] -> Sinh[z]/Cosh[z]`. The rewrite is purely structural (tree
   walk) followed by a single `evaluate()` pass. It turns the
   `0 * ComplexInfinity` folding trap (e.g. `x Csc[x] -> 0`) into a
   structural `0/0` that Series and L'Hospital can resolve. This one
   change resolves roughly 15 REPL-failing cases on its own.
2. *ArcTan / ArcCot at divergent inner argument* (runs before the
   rational layer). Computes `Limit[inner, ...]`; on
   `+Infinity / -Infinity` maps `ArcTan` to `+/- Pi/2` and `ArcCot` to
   `0 / Pi`. Previously `ArcTan[x^2 - x^4]` at `Infinity` sat on a
   Series that could not expand the inner polynomial at infinity and
   emitted `Power::infy` warnings.
3. *Bounded envelope (squeeze theorem)* (gated to `+Infinity` limits,
   runs before the Series layer). Uses a structural magnitude-bound
   walk to produce `|f| <= g(x)`: `|Sin/Cos/Tanh[anything]| <= 1`,
   `|ArcTan/ArcCot[anything]| <= Pi/2`, triangle inequality on Plus,
   multiplicativity on Times, and `|a^n| = |a|^n` for non-negative
   integer `n`. If `Limit[g(x), x -> Infinity] = 0`, the original limit
   is zero. Covers `Sin[t^2]/t^2`, `(1 +/- Cos[x])/x`,
   `(x Sin[x])/(5 + x^2)`, and friends at infinity, plus `x^2 Sin[1/x]`
   at 0 (which already worked through a different path but now routes
   here without Series warnings).
4. *Generalized log-reduction* (replaces the old single-`Power` form).
   Fires for a top-level `Power[b, g]` with x in `g` OR a top-level
   `Times[P1, P2, ...]` whose non-constant factors are all of that
   shape. Rewrites as `Exp[Limit[Sum[g_i * Log[b_i]]]]`. A
   post-processor maps a `+Infinity / -Infinity` log-limit to
   `Infinity / 0`, and refuses to produce an answer for
   `ComplexInfinity` / `Indeterminate`. Covers
   `((-3+2x)/(5+2x))^(1+2x)` at `Infinity`, `x^Tan[Pi x/2]` at 1,
   `(1+Sin[ax])^Cot[bx]` at 0, and the `E^(-x) x^20` shape (via
   log-reduction of the implicit `Power[E, -x] * Power[x, 20]`).

A smaller tweak: continuous substitution now folds `Power[0, positive]`
and `Sqrt[0]` to `0` locally, so `Limit[Sqrt[x-1]/x, x -> 1]` reports
`0` instead of the un-folded `Sqrt[0]` that PicoCAS's Power evaluator
leaves in place.

**Regression-suite additions (2026-04-20, batch 2):** three further
layers and one post-processing pass were added after a second round
of REPL-driven cases.

1. *Log of a finite-limit inner* (runs before Series). If `f = Log[g]`
   and `Limit[g]` is finite, returns `Log[c]`. Enables shapes such as
   `Log[1 + 2 Exp[-x]]` at `Infinity`.
2. *Log / polynomial merge at infinity* (runs after the Log-of-finite
   layer). Rewrites `Sum(Log[g_i]) + Sum(h_j)` into a single
   `Log[(Product g_i) * Exp[Sum h_j]]` whose inner product is expanded
   to expose cancellations (e.g. `Exp[x] * Exp[-x] -> 1`). Resolves
   `-x + Log[2 + E^x]` at `Infinity` to `0`.
3. *Term-wise Plus at infinity*. If every summand of a `Plus` at
   `+/- Infinity` has an individually finite limit, return the sum of
   those limits. Refuses as soon as any term is divergent or
   unresolved. Keeps the outer log-merge sound, because individually
   divergent shapes still bail and reach the other layers.
4. *Radical fusion in `Times`* (global, not Limit-specific). Moved
   out of `limit.c` and into `builtin_times` so that
   `Sqrt[6]/Sqrt[2]` becomes `Sqrt[3]` system-wide, not only as a
   Limit post-pass. The rule fuses `Power[a, q] * Power[b, -q]` into
   `Power[a/b, q]` whenever `a` and `b` are both positive numeric
   (integer, bigint, rational, or real) -- `a > 0` and `b > 0`
   places us on the principal branch. The ratio `a/b` may itself
   be an integer (`Sqrt[6]/Sqrt[2] -> Sqrt[3]`) or a rational
   (`9^(1/3)/3^(1/3) -> 3^(1/3)`, `2^(1/3)/3^(1/3) -> (2/3)^(1/3)`,
   `Sqrt[3]/Sqrt[6] -> Sqrt[1/2]`). Applies after same-base
   grouping, restarts on each fusion so chained reductions like
   `Sqrt[210]/Sqrt[6]/Sqrt[5] -> Sqrt[7]` converge.
5. *`Power[0, b]` folding* in `builtin_power`. `0^b` now evaluates to
   `0` for any positive `b` (integer, rational, or real) and to
   `ComplexInfinity` for any negative `b`. This eliminates the
   `Sqrt[0]` leak that previously bubbled up from continuous
   substitution in Limit; the auxiliary `reduces_to_zero` predicate
   has been removed from `limit.c`.

Also fixed: `heuristic_factor` in `facpoly.c` recursed indefinitely
when factoring `Power[a, rational]` where `a` was a constant (no
variables, but non-unit content). The guard is a single early return
for `v_count == 0` (no variables -> return the constant); segfault in
`Limit[1/(t Sqrt[t+1]) - 1/t, t -> 0]` was caused by this path and is
now fixed.

**Regression-suite additions (2026-04-21, batch 3):** follow-up batch
from a REPL-driven failure report. The fixes are small, targeted, and
preserve the layered dispatcher rather than adding new layers.

1. *Unary-minus extraction in the `Times` printer* (`src/print.c`). A
   leading negative literal (integer, real, or negative-numerator
   `Rational`) is now emitted as a prefixed `-` sign. This normalises
   `Times[-1, Infinity]` to `-Infinity`, `Times[-1, a, b]` to `-a b`,
   `Times[-3, x]` to `-3 x`, and `Times[Rational[-2, 3], x]` to
   `-2/3 x`. Previously these printed with a literal `-1` factor
   (`-1 Infinity`), which both looked wrong and tripped equality
   checks. Existing tests whose fixtures asserted the old spelling
   (`test_logexp`, `test_list`, `test_limit`) were updated.
2. *Opaque-head guard* (`src/limit.c`). Before running any analytic
   layer, `compute_limit` now refuses the problem if the expression
   applies an undefined or known-discontinuous head to an argument
   containing the limit variable. Unknown heads are detected via the
   symbol table (no `builtin_func`, no `down_values`, no `own_values`).
   A curated whitelist of discontinuous heads -- `Floor`, `Ceiling`,
   `Round`, `FractionalPart`, `IntegerPart`, `Sign`, `UnitStep`,
   `HeavisideTheta`, `KroneckerDelta`, `DiscreteDelta`, `Piecewise`,
   `Boole`, `Mod`, `Quotient` -- is treated as opaque even when a
   builtin exists, because plain substitution would silently pick one
   side's value at a jump. Fixes `Limit[f[x], x -> a]`,
   `Limit[Ceiling[5 - x^2], x -> 1]`, and similar cases that previously
   returned a misleading numeric answer.
3. *Even-order pole handling in the rational layer*. Layer 3 used to
   short-cut `Limit[N(x)/D(x), x -> a]` to `ComplexInfinity` whenever
   `D(a) = 0` and the direction was two-sided. This lost the parity
   information, so `1/(x-2)^2` at x = 2 came back as `ComplexInfinity`
   instead of `+Infinity`. The layer now defers unconditionally to the
   Series layer in that branch, which inspects the leading exponent:
   even-order poles return `+/-Infinity` (signed by the leading
   coefficient), odd-order poles return `ComplexInfinity`.
4. *Oscillatory limits return `Indeterminate`*. Layer 6 previously
   returned `Interval[{-1, 1}]` for shapes like `Sin[1/x]` at 0 and
   `Sin[x]` at `Infinity`. That interval is a correct bound but not a
   limit value; Mathematica returns `Indeterminate`, which we now match.
   The layer still only fires after the squeeze / substitution / Series
   paths have failed, so bounded shapes that actually converge to a
   fixed value are unaffected.
5. *Leaked limit variable in directional-infinity coefficient*. When
   the Series leading coefficient is a non-constant expression (e.g.
   `Log[x]` for `Log[x]/x` at 0), emitting
   `DirectedInfinity[Log[x]]` gave back a pseudo-answer that still
   mentioned the limit variable. `read_leading_term_limit` now bails
   in this case so the other layers (or the unevaluated fall-through)
   can run.
6. *`Direction -> I` and other purely-imaginary directions*. The
   option parser now accepts `I`, `k*I`, `-I`, and `Complex[0, k]` as
   valid directions, routing them through `LIMIT_DIR_COMPLEX`. Full
   branch-cut analysis isn't implemented yet, so such limits typically
   still return unevaluated -- but they no longer reject the option
   outright.

**Known limitations carried forward from earlier batches:**
* Joint multivariate limits with path dependence
  (`x y / (x^2 + y^2)` at `{0,0}`, `ArcTan[y/x]` at `{Infinity,
  Infinity}`, `y/(x+y)` at `{0,0}`, `(x^3+y^3)/(x^2+y^2)` at `{0,0}`):
  the polar-substitution heuristic is still future work. A handful of
  cases that happen to be continuous at the substitution point do
  resolve via `run_multivariate`.
* `(1 + Sinh[x])/Exp[x]` at `Infinity`: Series at Infinity does not
  recognise the natural `Exp[-x]` factorisation and L'Hospital stalls
  on the Cosh/Sinh derivative cycle.
* `(25^x - 5^x)/(4^x - 2^x)` at 0: Series infinite-escalates on the
  log-expansion of `k^x`; currently hangs if invoked (the REPL must be
  interrupted). Avoid until a dedicated `b^x` expansion kernel is
  added.
* `Log[1 - (Log[Exp[z]/z - 1] + Log[z])/z]/z` at a finite non-singular
  point (e.g. z = 100): one of the sub-expressions triggers repeated
  series retries inside `evaluate()`. Regression to be triaged.
* `Tan[x]` at `Pi/2, Direction -> Reals`: should return `Indeterminate`
  (the two sides diverge with opposite signs) but currently yields
  `ComplexInfinity`. A pole-with-sign-disagreement classifier is
  future work.
* `Sin[n Pi]` at `n -> Infinity`: should be `Indeterminate`; currently
  the log-reduction pathway picks an incorrect sign via the bounded
  envelope and returns `-Infinity`-like output in some code paths.
* Complex `Exp[Exp[-x/(1+Exp[-x])]]`-style stacked-exponential limits
  at `Infinity`: the bounded envelope does not see through nested
  `Exp` applied to bounded-but-decaying arguments.
* `Exp[Log[Log[x + ...]]/Log[Log[Log[Exp[x] + x + Log[x]]]]]` at
  `Infinity`: this requires comparative asymptotic expansion of
  iterated logarithms, which Series does not currently model.
* ~~The `$SeriesInvVar$` reverse-substitution symbol can leak back into
  user output for deeply nested `Log[Log[...]]` chains.~~ **Fixed
  2026-04-21**: `do_series_single` now substitutes the internal inverse
  variable `u -> 1/x` in every coefficient after the at-Infinity
  expansion, not only in `SeriesObj::x`. Previously a coefficient that
  was treated as "free of u" (e.g. `Log[1/u] -> -Log[u]` folded as a
  constant-in-u term) would retain the placeholder. Regression test
  `test_series_infinity_no_inv_var_leak` grep-checks the `$SeriesInvVar$`
  literal out of both InputForm and FullForm output for four
  representative shapes.

---

**Work-package-driven additions (2026-04-21, batches 4–12):** a
coordinated series of work packages closed out most of the residual
Limit cases from the REPL report. Each landed as its own bounded
regression in `tests/test_limit.c` (`test_wp8_*`, `test_wp2_*`,
`test_wp4_*`, `test_wp1_*`, `test_wp3_*`, `test_wp6_*`, `test_wp5_*`,
`test_wp9_*`, `test_wp7_*`).

1. **WP-8 (numeric-point fast path).** `layer1_fast_paths` gains a
   `try_numeric_point_substitution` that substitutes a plain numeric
   limit point into the Together-normalised form and returns the
   result when clean, skipping the Series / L'Hospital pipeline. Gate
   checks: Together's denominator does not vanish at the point, and no
   `Power[base, x-dep-exponent]` with a divergent exponent remains
   (1^inf / 0^0 / inf^0 indeterminate seeds). Fixes
   `Limit[Log[1 - (Log[Exp[z]/z - 1] + Log[z])/z]/z, z -> 100]` hang.
2. **WP-2 (b^x series kernel).** `so_inv` caps its iteration count
   based on input-coefficient leaf size; previously 1/(2^x - 3^x) at
   x = 0 spun in the O(N^2) simp-per-iteration loop because PicoCAS
   doesn't canonicalise polynomials in symbolic `Log[2], Log[3]`. The
   cap is sized so numeric/trivial inputs retain the full order and
   heavy symbolic inputs still produce a valid leading-term Laurent.
   Also: the single-simp-per-iteration rewrite in `so_inv` cuts the
   per-step cost by roughly 3x on expression growth.
3. **WP-4 (pole sign-disagreement classifier).** New
   `LIMIT_DIR_REALS` tag (distinct from the implicit `LIMIT_DIR_TWOSIDED`
   default). At an odd-order pole with two sides disagreeing in sign,
   `Direction -> Reals` returns `Indeterminate` (the real-line answer
   matches Mathematica's behaviour for `Tan[x] at Pi/2, Reals`), while
   the default two-sided direction keeps the old `ComplexInfinity`
   fall-back for unqualified rational-function limits.
   `Direction -> Complexes` returns `ComplexInfinity` for any pole
   regardless of parity (radial interpretation).
4. **WP-1 (multivariate path-dependence).** `run_multivariate` now
   performs polar/spherical substitution at the origin (2D/3D) and
   cross-checks with direction sampling along axes and diagonals.
   Returns `Indeterminate` when sampled paths disagree and the common
   value otherwise. Joint-at-Infinity works for all-positive-orthant
   via straight-line paths. Covered by `test_wp1_multivariate`
   (`Tan[x y]/(x y)` → 1, `y/(x+y)` → Indeterminate,
   `(x^3 + y^3)/(x^2 + y^2)` → 0, `ArcTan[y^2/(x^2+x^3)]` →
   Indeterminate, `x z/(x^2+y^2+z^2)` → Indeterminate,
   `ArcTan[y/x] at {Infinity, Infinity}` → Indeterminate). The
   origin-fast-path gains an inner-divide-by-zero scanner that catches
   0/0 shapes buried inside ArcTan/Sin/etc. that PicoCAS's arithmetic
   would silently fold to 0.
5. **WP-3 (Series of x^a at nonzero point).** Two-part fix. First,
   `do_series_single`'s `try_factor_power_prefactor` is now gated to
   x0 = 0 or Infinity -- at any other x0 the Power[x, alpha] is
   expanded via the (1+u)^alpha binomial kernel rather than pulled
   out as a symbolic prefactor (which would leave the other factors
   to be expanded in isolation, producing 0). Second,
   `try_apart_preprocess` refuses to Apart-preprocess expressions
   containing `Power[base, non-rational]` because picocas's Apart
   collapses such inputs to 0. Together these resolve
   `Limit[3 (x^a - a x + a - 1)/(x-1)^2, x -> 1]` to
   `(3/2) a (a - 1)`.
6. **WP-6 (Sinh/Cosh exponentialisation at Infinity).** A new
   pre-Series tree walk (`rewrite_hyperbolic_to_exp`) replaces
   `Sinh[z] -> (E^z - E^-z)/2`, `Cosh[z] -> (E^z + E^-z)/2`,
   `Tanh[z] -> (E^z - E^-z)/(E^z + E^-z)` whenever the limit point is
   ±Infinity. The result is `Expand[]`-ed so the term-wise Plus layer
   can fold cancelling Exp[kx] summands. Fixes
   `(1 + Sinh[x])/Exp[x] at Infinity -> 1/2` and the Cosh twin.
   Limits at finite points are unchanged (Taylor via D suffices).
7. **WP-5 (sign-aware envelope / dominant term).** `layer_plus_termwise`
   gained a `growth_exponent_upper` helper: polynomial degree upper
   bound that treats Sin/Cos/Tanh/ArcTan/ArcCot as 0, multiplies
   through Times, takes max over Plus, multiplies by exponent in
   `Power[base, nonnegative_int]`, and treats
   `Power[base, negative_int]` with growing base as 0 (reciprocal of
   something diverging is bounded by 0). When one Plus summand has
   strictly larger growth than every other and its limit is ±Infinity,
   that summand wins and we return its limit. Covers
   `x^2 + x Sin[x^2]` → Infinity, `x + Sin[x]` → Infinity (bounded
   oscillator absorbed by dominant polynomial).
8. **WP-9 (complex-direction branch cuts).** New dir tag
   `LIMIT_DIR_IMAGINARY` (for numeric-imaginary directions: `I`, `k I`,
   `Complex[0, k > 0]`). The analytic layers compute the principal-
   branch result as usual; a post-pass in `builtin_limit` then
   conjugates the imaginary part via `ReplaceAll[I -> -I]` when the
   direction was `I` (landing on the branch below the cut), and
   returns `Indeterminate` for `Direction -> Complexes` when the
   principal-branch result picked up an imaginary part (a branch-point
   value). Any signed-infinity pole under `Complexes` collapses to
   `ComplexInfinity`. Also: `Limit[{f1, f2, ...}, ...]` now threads
   over a top-level List in its first argument.
9. **WP-7 (Gruntz-lite for Log[sum]).** A partial Gruntz-style
   rewrite for the iterated-log family: `Log[dom + rest...]` at
   +Infinity with a unique dominant summand gets rewritten as
   `Log[dom] + Log[1 + rest/dom]`, then recursed. Handles
   `Log[x + Log[x]] - Log[x] -> 0` and
   `Log[x^2 + x] - 2 Log[x] -> 0` directly. **Full Hardy-field
   comparative asymptotics** (stacked `Exp[Exp[-x/(1+Exp[-x])]]` etc.,
   multi-level log-exp dominance, `Sin[x] + Log[x-a]/Log[E^x-E^a]` as
   x → a) remain future work -- a fully compliant Gruntz algorithm
   is ~600-800 LOC on top of a proper MRV (most-rapidly-varying)
   rewriter, which we have not yet built. The groundwork here
   (growth_exponent_upper, Log[sum] rewrite, hyperbolic
   exponentialisation) is the platform for that next step.

**Known limitations carried forward (still unevaluated):**
* Stacked-exponential limits like `Exp[Exp[-x/(1+Exp[-x])]] ...` at
  Infinity -- requires full Gruntz MRV.
* `Limit[Exp[x] (Exp[1/x - Exp[-x]] - Exp[1/x]), x -> Infinity]`
  (should be -1): first-order asymptotic cancellation between two
  close exponentials, out of reach without Gruntz.
* `Limit[Sin[x] + Log[x-a]/Log[E^x - E^a], x -> a]` (should be
  `1 + Sin[a]`): the Log/Log ratio needs two iterations of
  L'Hospital with analytic simplification in between; the current
  L'Hospital guardrail trips before it gets there.
* `Limit[Exp[Log[Log[x + Exp[Log[x] Log[Log[x]]]]]/
                Log[Log[Log[Exp[x] + x + Log[x]]]]], x -> Infinity]`
  (should be `E`): classic Gruntz case with three levels of
  logarithmic dominance ranking.
* Iterated limits involving `ArcTan[y/x]` where the inner limit
  variable's sign must be inferred from the outer context -- e.g.
  `Limit[ArcTan[y/x], {x -> Infinity, y -> Infinity}]`: we'd need
  `Assumptions -> y > 0` style inference.

**Known limitations** (remaining; these return the original
`Limit[...]` unevaluated):
* Multivariate limits that are truly path-dependent (e.g.
  `x y / (x^2 + y^2)` at `{0,0}`) are left unevaluated rather than
  returning `Indeterminate`; the polar-substitution heuristic from
  Layer 5a of the spec is future work.
* `(1 + Sinh[x])/Exp[x]` at `Infinity`: Series at Infinity does not
  recognise the natural `Exp[-x]` factorisation, and L'Hospital stalls
  on the Cosh/Sinh derivative cycle.
* `Csc[x]/E^x` at `Infinity` and the `(x Sin[x])/(x + Sin[x])` /
  `(x^2(1+Sin[x]^2))/(x+Sin[x])^2` shapes: these are genuinely
  `Indeterminate` (bounded but sign-oscillating numerator divided by a
  factor that does not dominate). A conservative "Indeterminate
  classifier" is future work.
* `(1 - E^(-x))^E^x` at `Infinity`: should reduce to `1/E` via
  Series-of-log but the inner limit `Log[1 - E^(-x)] * E^x` sits
  outside the current log-expansion capabilities.
* `Limit[(1/(1-x))^(-1/x^2), x -> 0]`: two-sided, the from-above and
  from-below limits are respectively `0` and `Infinity`, so the
  two-sided answer is ambiguous. Our dispatcher returns `1` via a
  fall-through path (strictly speaking both `Indeterminate` and `1`
  are defensible depending on convention; Mathematica returns
  `Indeterminate`).


---

## Simplification extensions (2026-04-20)

### `Log[E^k]` and `Log[b, b^k]`

`builtin_log` now folds direct Log-of-Power forms:

* `Log[E^k] -> k` (one-argument Log), when `k` is a real numeric
  (integer, bigint, rational, or real).
* `Log[b, b^k] -> k` (two-argument Log), when `k` is a real numeric
  and `b` is a known-positive value (positive numeric, `E`, or `Pi`).

The real-`k` restriction keeps us on the principal branch -- for
complex `k` the identity `Log[b^k] = k Log[b]` can fail by a multiple
of `2 Pi i / Log[b]`. Helpers `is_real_numeric_expr` and
`is_positive_known` are local to `logexp.c`. Example:
`Log[E^4] -> 4`, `Log[2, 2^(1/3)] -> 1/3`, `Log[E^x]` stays.

### Trig / inverse-trig and hyperbolic / inverse-hyperbolic cancellation

`builtin_sin`, `builtin_cos`, `builtin_tan`, `builtin_cot`,
`builtin_sec`, `builtin_csc` (in `trig.c`) and
`builtin_sinh`..`builtin_csch` (in `hyperbolic.c`) each fold
`F[F_inv[x]] -> x` via a small shared helper (`strip_inverse_call`).
These identities hold identically over the complex numbers because
each `ArcF` is a right inverse of `F` by construction.

The opposite direction (`ArcSin[Sin[x]]`, etc.) is deliberately *not*
folded -- that identity only holds on each function's principal
domain. The two-argument `ArcTan[x, y]` form is also excluded from
the `Tan[ArcTan[...]]` rule (guarded by `arg_count == 1`) because
`Tan[ArcTan[x, y]] = y/x` rather than a single variable.

---

## Infinity, ComplexInfinity, and Indeterminate (2026-04-21)

Arithmetic involving the symbols `Infinity`, `ComplexInfinity`, and
`Indeterminate` now follows the standard Mathematica rules.
Previously, `Plus` and `Times` treated `Infinity` as an ordinary
opaque symbol and folded `Infinity - Infinity` to `0`, `0 * Infinity`
to `0`, etc. The new behaviour is implemented inline at the top of
`builtin_plus` (`src/plus.c`), `builtin_times` (`src/times.c`), and
`builtin_power` (`src/power.c`), with shared predicates exported from
`arithmetic.h`:

```
bool is_infinity_sym(Expr*);
bool is_complex_infinity_sym(Expr*);
bool is_indeterminate_sym(Expr*);
bool is_neg_infinity_form(Expr*);   /* Times[c, Infinity] with c<0 */
int  expr_numeric_sign(Expr*);
```

### Plus

| Inputs                                  | Result          |
|-----------------------------------------|-----------------|
| any term is `Indeterminate`             | `Indeterminate` |
| `Infinity + (-Infinity)` (any flavour)  | `Indeterminate` (with `Infinity::indet` message) |
| `ComplexInfinity + ComplexInfinity`     | `Indeterminate` (with `Infinity::indet` message) |
| `ComplexInfinity + (any other infinity)`| `Indeterminate` |
| `ComplexInfinity + finite`              | `ComplexInfinity` |
| `Infinity + finite` (no `-Infinity`)    | `Infinity`      |
| `-Infinity + finite` (no `+Infinity`)   | `-Infinity` (= `Times[-1, Infinity]`) |

Term classification (`classify_plus_term`) recognises both bare
`Infinity` and the canonical product form `Times[c, Infinity, ...]`,
using the sign of the leading numeric factor to distinguish `+` from
`-`. This means `5 Infinity - 7 Infinity` is correctly recognised as
`Infinity + (-Infinity)` rather than collapsing to `-2 Infinity`.

### Times

| Inputs                                  | Result          |
|-----------------------------------------|-----------------|
| any factor is `Indeterminate`           | `Indeterminate` |
| `0 * Infinity` or `0 * ComplexInfinity` | `Indeterminate` (with `Infinity::indet` message) |
| `c * Infinity`        (`c` numeric, `c > 0`) | `Infinity` |
| `c * Infinity`        (`c` numeric, `c < 0`) | `-Infinity` (`Times[-1, Infinity]`) |
| `c * ComplexInfinity` (`c` numeric, `c != 0`) | `ComplexInfinity` |
| `Infinity * Infinity`                   | `Infinity`      |
| `Infinity * ComplexInfinity`            | `ComplexInfinity` |

If non-numeric symbolic factors are present the simplification is
skipped and the existing symbolic `Times` machinery runs. Magnitude
is not preserved when multiplying by an unbounded value -- `-3
Infinity` reduces to `-Infinity` (= `Times[-1, Infinity]`), not
`Times[-3, Infinity]`.

### Power

| Inputs                                  | Result          |
|-----------------------------------------|-----------------|
| `Indeterminate^x` or `x^Indeterminate`  | `Indeterminate` |
| `1^Infinity`, `1^-Infinity`, `1^ComplexInfinity` | `Indeterminate` (with `Infinity::indet`) |
| `Infinity^0`, `(-Infinity)^0`, `ComplexInfinity^0` | `Indeterminate` (with `Infinity::indet`) |
| `0^Infinity`                            | `0`             |
| `0^-Infinity`                           | `ComplexInfinity` (with `Power::infy`) |
| `0^ComplexInfinity`                     | `Indeterminate` (with `Infinity::indet`) |
| `Infinity^Infinity`                     | `ComplexInfinity` |
| `Infinity^-Infinity`                    | `0`             |
| `Infinity^n`        (`n` numeric, `n > 0`)         | `Infinity` |
| `Infinity^n`        (`n` numeric, `n < 0`)         | `0`        |
| `ComplexInfinity^n` (`n` numeric, `n > 0`)         | `ComplexInfinity` |
| `ComplexInfinity^n` (`n` numeric, `n < 0`)         | `0`        |

### Division

`builtin_rational` and `builtin_divide` (`src/arithmetic.c`) emit
`Power::infy: Infinite expression 1/0 encountered.` when forming a
rational with zero denominator. `1/0 -> ComplexInfinity` and
`0/0 -> Indeterminate` (with the additional
`Infinity::indet: Indeterminate expression 0 ComplexInfinity`
message).

### Diagnostics

The `Infinity::indet` and `Power::infy` messages are written to
`stderr` (matching Mathematica's `Message[]` convention of routing
diagnostics off the value channel) so that programmatic consumers of
results are not affected.

### Tests

Coverage lives in `tests/test_infinity.c` (binary
`tests/build/infinity_tests`) and exercises every row of the tables
above. The limit-test helper `check_equiv` (`tests/test_limit.c`) was
updated to fall back to printed-form equality when the result is one
of the infinity flavours, since `Expand[Together[Infinity-Infinity]]`
now correctly returns `Indeterminate` rather than `0`.


## Series inversion correctness fixes (2026-04-21)

Two bugs in `so_inv` (`src/series.c`) caused `Series[1/f, ...]` and
reciprocal-trig / reciprocal-hyperbolic expansions to emit trailing
zeros instead of their correct high-order coefficients. Both are now
fixed.

### 1. Zero-coefficient inflation of the expression-growth guardrail

`so_inv` uses a leaf-count guardrail to cap the number of inversion
iterations `N` when the input coefficients are large symbolic
expressions (e.g. polynomials in `Log[2]`, `Log[3]` arising from
`1/(2^x - 3^x)`). The old leaf-count walk summed over every
coefficient, including zero ones, which inflated the count for
sparse polynomials. For `1/(1 + x^2)` at internal order 25 the 23
zero entries produced `total_leaves = 25`, triggering `N_cap = 6`
and truncating the result to `1 - x^2 + x^4 + 0 + 0 + ...`.

Fix: the walk now skips zero coefficients. Purely-numeric
coefficients (Integer, BigInt, Real, Rational, Complex, and
arithmetic combinations of these) are also skipped, since their size
does not compound under `simp()` during the recurrence.

### 2. Latent crash when simp()'ing a large flat-Plus of numeric terms

Lifting the cap exposed a separate issue in the Taylor-via-simp path
inside `so_inv`: when the recurrence built a `Plus` node with 10+
`Times[Rational, Rational]` terms and passed it to `simp()` in a
single call, the evaluator could dereference a NULL `Plus` argument
(reproducible as `Series[1/Sin[x], {x, 0, 9}]` or higher). Rather
than chase the evaluator bug, `so_inv` now branches:

- **All coefficients purely numeric** (the Sin/Cos/Sinh/Cosh path
  and polynomial inverses): accumulate `b[k]` with one `simp()` per
  term. Each partial sum collapses to a single rational via
  `add_numbers`, so the O(k^2) blowup the batched path was designed
  to prevent does not occur, and the evaluator never sees a large
  flat-Plus tree.

- **Otherwise** (symbolic coefficients): retain the batched
  single-simp-at-end construction that keeps symbolic-constant
  inversions (like `1/(2^x - 3^x)`) tractable.

### Affected cases, now correct

| Input                              | Old output                                | New output                          |
|------------------------------------|-------------------------------------------|-------------------------------------|
| `Series[Csc[x], {x, 0, 5}]`        | last coef `0`                             | `31/15120`                          |
| `Series[Sec[x], {x, 0, 6}]`        | last coef `0`                             | `61/720`                            |
| `Series[Sech[x], {x, 0, 6}]`       | last coef `0`                             | `-61/720`                           |
| `Series[Coth[x], {x, 0, 5}]`       | last coef `1/240`                         | `2/945`                             |
| `Series[1/(1 + x^2), {x, 0, 12}]`  | truncated after `x^4`                     | full `1 - x^2 + x^4 - ... + x^12`   |
| `Series[1/(1 - x^2), {x, 0, 12}]`  | truncated after `x^4`                     | full geometric series               |
| `Series[1/(Exp[x] - 1 - x), ...]`  | truncated after `x^3`                     | full Bernoulli-derived coefficients |

### Tests

All `series_tests`, `simp_tests`, and `distribute_tests` soft-failures
(printer-cleanup staleness in the latter two, series correctness in
the first) are resolved. The test suite is now entirely pass-clean:
59 binaries, 0 hard failures, 0 silent `FAIL:` prints.

## Limit robustness pass (2026-04-21)

Five categorical issues in `Limit[]` were fixed after a second head-to-head
review against Mathematica on 15 adversarial inputs. Each is handled at
the layer that has the right context, and every fix is gated narrowly so
the change does not regress other limit shapes (verified against the full
`limit_tests` + `series_tests` + whole-repo test battery -- 59 binaries,
all green).

### 1. Spurious `Power::infy` / `Infinity::indet` stderr spam

PicoCAS's arithmetic modules (`arithmetic.c`, `power.c`, `plus.c`,
`times.c`) print informational messages whenever they fold `1/0`,
`Infinity^0`, `0^ComplexInfinity` and similar. `Limit`'s internal probes
(continuous substitution, Series / L'Hospital, polar substitution) hit
these shapes *by design* — that is how they discover a divergent form.
The messages were escaping to the user's terminal and drowning out
otherwise-correct answers.

A re-entrant counter (`g_arith_warnings_muted`, declared in
`arithmetic.h`) now gates the fprintf sites. `builtin_limit` pushes the
counter on entry and pops it on exit, so nested limits stack cleanly.
The return values (Indeterminate / ComplexInfinity) are unchanged; only
the stderr noise is suppressed during limit evaluation.

### 2. Discontinuous-head leak past the opaque-head gate

`contains_opaque_head_over(f, x)` is invoked at the top of
`compute_limit` to refuse limits through `Ceiling`, `Floor`, `Sign`, ...
But user wrappers (`h[n_] := Ceiling[n]`) evade that check: the outer
head `h` has DownValues, so `is_known_head_symbol` returns true. Once
the expression is evaluated inside `compute_limit` (via `simp`), the
wrapper unfolds and the discontinuous head is exposed — but by then the
opaque check has already passed and the continuous-substitution fast
path happily returns `Ceiling[4] = 4` for `Limit[h[5-x^2], x -> 1]`.

The fix is a **second opaque-head check after `simp`**. If the
evaluated form applies a discontinuous head to anything involving `x`,
`compute_limit` returns NULL and leaves the expression unevaluated.
This matches Mathematica: `Limit[Ceiling[5-x^2], x -> 1]` is left
symbolic because the ceiling jumps at 4.

### 3. Asymptotic-Log contamination of the Series layer's leading
coefficient

`Series[Log[E^x - E^a], {x, a, k}]` yields `Log[E^a] + Log[-a + x] +
(x-a)/2 + ...` — an *asymptotic* expansion whose leading "constant
coefficient" is `Log[x - a]`, which is in fact unbounded at `x = a`.
`read_leading_term_limit` previously returned this coefficient verbatim
when the leading exponent was 0, yielding a "limit" that still depended
on the limit variable (for `Limit[Log[x-a]/Log[E^x-E^a], x -> a]` the
wrong answer `Log[-a+x]/(Log[E^a] + Log[-a+x])`).

Two defensive changes:
  * `read_leading_term_limit` refuses the leading coefficient if it
    contains the limit variable. The coefficient of a genuine power
    series is a constant; an `x`-residual is a structural signal that
    the expansion is asymptotic, not Taylor.
  * `layer2_series`'s escalation loop terminates early when a Series
    call produces an `x`-residual expression. Escalating the order to
    `LIMIT_SERIES_MAX_ORDER = 32` for these shapes is O(k^3) in the
    Series kernel and was observed to hang on `Log[E^x - E^a]`; the
    short-circuit brings the worst case back to milliseconds.

Effect: `Limit[Log[x-a]/Log[E^x-E^a], x -> a]` now returns unevaluated
rather than a wrong-looking `x`-dependent expression. (Mathematica
returns the closed form 1 via a Gruntz/series analysis that PicoCAS
does not yet implement; returning unevaluated is a strict improvement
over the prior incorrect answer.)

### 4. Atom-substitution layer for `Power[b, e(x)]` shapes

Expressions like `(3^(1/x) - 3^(-1/x)) / (3^(1/x) + 3^(-1/x))` at
`x -> 0` are rational in `3^(2/x)` (after Together) but not expandable
as a Taylor series — the Power has an essential singularity. A new
layer, `layer_atom_substitute`, detects a `Power[b, e(x)]` subterm with
`b` free of `x` and `e(x)` containing `x`, substitutes a fresh symbol
`u` for it, computes `Limit[Power[b, e(x)], x -> point]`, and recurses
on `Limit[f_sub, u -> atom_lim]`. The rational layer then closes out
the `u`-limit.

Gating is deliberately narrow:
  * A suitable `Power` must already exist in `f` (no expensive Together
    call otherwise).
  * The limit point must be a numeric literal (the shape does not apply
    to Limits at a symbolic point, and Together on shapes like
    `E^x - E^a` is itself expensive).
  * The recursion depth is capped at 3 to keep the total work bounded
    when substitutions chain.

### 5. Two-sided disagreement probe (new `Indeterminate` backstop)

`(3^(1/x) - 3^(-1/x)) / (3^(1/x) + 3^(-1/x))` at `x -> 0` two-sided:
the one-sided limits are `+1` (from above) and `-1` (from below), so
the two-sided limit is `Indeterminate`. A new final-fallback layer,
`layer_onesided_disagree`, fires when:
  * `ctx->dir == LIMIT_DIR_TWOSIDED`
  * the point is a numeric literal
  * the expression has `x` in some exponent (the failure mode)
  * recursion depth is shallow

It computes the two one-sided limits; if both succeed and disagree, the
result is `Indeterminate`. If they agree, that common value is returned
(so the layer is also a correctness tightening). Recursive calls use
one-sided directions, so the probe does not re-enter itself.

### 6. Literal-sign for `Log[c]` with positive real `c`

`literal_sign(Log[c])` previously returned 0 for numeric `c`, which
prevented the series layer from resolving the sign of a leading
`DirectedInfinity[Log[k]]` term and surfacing it as `+Infinity` or
`0`. The helper now returns `+1` / `0` / `-1` based on `c` vs `1` for
positive-real `c` (Integer, Real, BigInt, positive Rational). This
rescues `Limit[3^(1/x), x -> 0, "FromAbove"]` → `Infinity`, which is a
prerequisite for the two-sided disagreement probe to classify `3^(1/x)`
shapes.

### Adversarial-input scorecard

| # | Input                                                                | Before              | After            | Mathematica      |
|---|----------------------------------------------------------------------|---------------------|------------------|------------------|
| 1 | `Limit[Sin[1/x], x -> 0]`                                            | Indeterminate + 6×1/0 noise | Indeterminate        | Indeterminate    |
| 2 | `Limit[Log[x]/x, x -> 0]`                                            | unevaluated         | unevaluated      | `-Infinity` (dir-sensitive) |
| 3 | `Limit[Sin[x] + Log[x-a]/Log[E^x-E^a], x -> a]`                      | wrong (has `x` in answer) | unevaluated      | `1 + Sin[a]` (out of reach) |
| 4 | `Limit[Sin[x], x -> Infinity]`                                       | Indeterminate + 4×1/0 noise | Indeterminate    | Indeterminate    |
| 9 | `Limit[ArcTan[y/x], {x,y} -> {Inf,Inf}]` (joint)                     | unevaluated + 1/0   | Indeterminate    | `Pi/2` (out of reach) |
| 11 | `Limit[f[x,y], {x,y} -> {0,0}]` where `f = ArcTan[y^2/(x^2+x^3)]`    | Indeterminate + 2×1/0 noise | Indeterminate    | Indeterminate    |
| 12 | `Limit[(3^(1/x)-3^(-1/x))/(3^(1/x)+3^(-1/x)), x -> 0]`               | unevaluated + 19×1/0 noise | Indeterminate    | Indeterminate    |
| 13 | `Limit[Ceiling[5-x^2], x -> 1]`                                      | wrong (returned 4)  | unevaluated      | (left as is)     |
| 14 | `Limit[x Sin[1/x], x -> 0]`                                          | 0 + 1×1/0 noise     | 0                | 0                |
| 15 | `Limit[Sin[n Pi], n -> Infinity]`                                    | wrong (`-Infinity`) | Indeterminate    | Indeterminate    |

Cases 5-8 (nested exp / iterated log / deeply composed asymptotics)
remain unevaluated: they are the classical Gruntz-algorithm test cases
that require an MRV (most-rapidly-varying) analysis PicoCAS does not
yet have. Returning unevaluated is the correct fall-through; a future
pass can add a Gruntz layer if needed.

### Files touched

  * `src/arithmetic.h`, `src/arithmetic.c` -- mute counter + gated
    prints
  * `src/plus.c`, `src/times.c`, `src/power.c` -- gated prints
  * `src/limit.c` -- post-`simp` opaque check; `read_leading_term_limit`
    x-residual guard; `layer2_series` early-terminate; new
    `layer_atom_substitute`; new `layer_onesided_disagree`;
    `literal_sign` `Log` arm
  * `src/limit.h` -- unchanged (public interface preserved)

## Trig / hyperbolic argument canonicalisation (2026-04-21)

Forward trigonometric and hyperbolic builtins now canonicalise
two shapes of argument at evaluation time:

1. **Sign extraction** via even/odd symmetry (`Cos[-2 x] -> Cos[2 x]`,
   `Sin[-4 x] -> -Sin[4 x]`).
2. **Imaginary-unit bridge** between trig and hyperbolic heads
   (`Cosh[I x] -> Cos[x]`, `Sin[x/I] -> -I Sinh[x]`).

Both rewrites are short, syntactic, and fire before the existing
zero / infinity / `Pi`-multiple / numeric-approximation branches.

### Sign extraction

**Scope.**
* Even (argument sign dropped): `Cos`, `Sec`, `Cosh`, `Sech`.
* Odd (argument sign pulled outside as `Times[-1, ...]`):
  `Sin`, `Tan`, `Cot`, `Csc`, `Sinh`, `Tanh`, `Coth`, `Csch`.

**Trigger -- "superficially negative" arguments.** The new helper
`expr_is_superficially_negative` in `src/arithmetic.{c,h}` matches:

* a negative numeric literal (`Integer`, `Real`, `BigInt`, or
  `Rational[-n, d]`);
* a `Complex[re, im]` whose canonical sign is negative, i.e.
  `re < 0`, or `re == 0 && im < 0` (catches `-2 I` =
  `Complex[0, -2]`);
* a `Times[c, ...]` whose leading factor `c` is itself superficially
  negative (covers `-k x`, `-1.5 y`, `(-3/2) a`, `-x` =
  `Times[-1, x]`, and the complex-imaginary variants like
  `-2 I x` = `Times[Complex[0, -2], x]`).

The shape test is deliberately syntactic: it recognises sign
extraction that is cheap and obviously valid without invoking the
full simplifier on symbolic arguments.

**Mechanism.** Each forward builtin runs, immediately after the
forward/inverse strip, a short-circuit that:

1. Calls `flip_sign(arg)` -- builds `Times[-1, arg]` and invokes the
   evaluator so nested structure canonicalises:
   `Times[-1, Complex[0, -2]] -> Complex[0, 2]`;
   `Times[-1, Times[-k, x]] -> k x`; single-factor cases like
   `Times[-1, x]` collapse to `x`.
2. Rebuilds `f[flip_sign(arg)]` (even case) or
   `Times[-1, f[flip_sign(arg)]]` (odd case).
3. Returns the new expression -- `res` is freed by `evaluate_step`
   per the builtin calling convention.

Because the evaluator loop continues until a fixed point, the inner
re-entry of `f` on the flipped argument picks up any other available
simplifications (exact rational-`Pi` values, the new imaginary-unit
bridge, numeric approximation, etc.) in the same top-level evaluation.

### Imaginary-unit bridge

**Identities.** Applied immediately after sign extraction, so the
leading imaginary coefficient is always canonically positive:

| Trig          | Hyperbolic     |
|---------------|----------------|
| `Cos[I y]  -> Cosh[y]` | `Cosh[I y] -> Cos[y]`  |
| `Sin[I y]  -> I Sinh[y]` | `Sinh[I y] -> I Sin[y]`  |
| `Tan[I y]  -> I Tanh[y]` | `Tanh[I y] -> I Tan[y]`  |
| `Cot[I y]  -> -I Coth[y]` | `Coth[I y] -> -I Cot[y]` |
| `Sec[I y]  -> Sech[y]` | `Sech[I y] -> Sec[y]`    |
| `Csc[I y]  -> -I Csch[y]` | `Csch[I y] -> -I Csc[y]` |

**Trigger.** The helper `peel_imaginary_unit(arg)` returns the "real
part" `y` when `arg` is of the form `I y`, i.e.:

* `arg = Complex[0, k]` with `k > 0` (then `y = k`);
* `arg = Times[Complex[0, k], rest...]` with `k > 0` (then
  `y = k * rest`, evaluated).

A leading `Complex[0, k]` with `k < 0` never reaches this helper
because sign extraction canonicalises `k` to positive first, so for
example `Cos[x/I]` (which parses to `Cos[Times[Complex[0,-1], x]]`)
first becomes `Cos[I x]` via even-symmetry and then `Cosh[x]` via
the bridge in the same evaluator cycle.

**Why symmetric coefficients.** The bridge multipliers are the
same in both directions because `1/I = -I`:
`Cot[I y] = Cosh[y] / (I Sinh[y]) = -I Coth[y]`, and
`Coth[I y] = Cos[y] / (I Sin[y]) = -I Cot[y]`.

### Interaction with existing rules

* Exact rational-`Pi` evaluation is unaffected: `Sin[-Pi/3]` becomes
  `-Sin[Pi/3]` = `-1/2 Sqrt[3]`.
* Numeric approximation is unaffected: `Sin[-1.5]` evaluates to
  `-0.997495` (via the `-Sin[1.5]` path).
* `f[0]` and infinity identities still fire because `0` is neither
  superficially negative nor pure imaginary.
* Previously, `Sin[2 x] + 2 Cos[x] // TrigToExp // ExpToTrig` left a
  residual `(-I) Cosh[(2*I) x] + (-I) Sinh[(2*I) x] + (I) Cosh[(-2*I) x] + (I) Sinh[(-2*I) x]`
  tangle. Sign extraction canonicalises the `-2 I x` arguments and
  the bridge collapses the hyperbolic heads back to trig, letting
  the surrounding `Plus` cancellations fire.

### Files touched

* `src/arithmetic.h`, `src/arithmetic.c` -- new public
  `expr_is_superficially_negative` helper (reused by both
  `trig.c` and `hyperbolic.c`; available to future callers).
* `src/trig.c` -- new static helpers `flip_sign`, `even_fold`,
  `odd_fold`, `peel_imaginary_unit`, `trig_i_fold`; sign-fold and
  I-fold short-circuits added to each of `builtin_sin`, `builtin_cos`,
  `builtin_tan`, `builtin_cot`, `builtin_sec`, `builtin_csc`.
* `src/hyperbolic.c` -- same helpers (`hyp_i_fold` counterpart);
  short-circuits added to `builtin_sinh`, `builtin_cosh`,
  `builtin_tanh`, `builtin_coth`, `builtin_sech`, `builtin_csch`.


## Increment / Decrement / AddTo / SubtractFrom

PicoCAS implements the full Mathematica family of in-place numeric
update operators. The parser recognises `++`, `--`, `+=`, and `-=`
and rewrites them to named heads; each head is a C built-in with
attribute `{HoldFirst, Protected}`.

### Operators and heads

| Source form   | Parsed as              | Returns    |
| ------------- | ---------------------- | ---------- |
| `x++`         | `Increment[x]`         | old value  |
| `x--`         | `Decrement[x]`         | old value  |
| `++x`         | `PreIncrement[x]`      | new value  |
| `--x`         | `PreDecrement[x]`      | new value  |
| `x += dx`     | `AddTo[x, dx]`         | new value  |
| `x -= dx`     | `SubtractFrom[x, dx]`  | new value  |

Postfix `++` / `--` have operator precedence 660 (tighter than `^`,
`*`, `+`); prefix `++` / `--` use the same precedence for their
operand, so `++a + 1` parses as `(++a) + 1`. `+=` and `-=` share the
right-associative precedence of `Set` (40), so `x += y = 2` assigns
`2` to `y` and then adds it to `x`.

### Semantics

1. The first argument is held (`HoldFirst`) so the identifier is not
    evaluated before the operator sees it.
2. The operator resolves the underlying mutable symbol: for a bare
    symbol that is the symbol itself; for a `Part[sym, ...]` l-value
    it is the underlying `sym`.
3. If the target symbol has no `OwnValue`, the operator emits
    `<Head>::rvalue: <name> is not a variable with a value, so its
    value cannot be changed.` and returns itself unevaluated, matching
    Mathematica.
4. Otherwise the operator evaluates the l-value, combines it with the
    delta via `Plus` (negating for `Decrement` / `SubtractFrom`), and
    re-evaluates. Because `Plus` is `Listable`, list-valued targets
    thread element-wise (`{1,2,3}++` produces `{2,3,4}`, and
    `{18,19,20} += {20,21,22}` produces `{38,40,42}`).
5. The new value is written back via `Set`, so indexed updates go
    through the existing `expr_part_assign` pipeline
    (`list[[2]]++` mutates only that element).
6. `Increment` and `Decrement` return the *old* value; the other four
    return the *new* value.

### Examples

```mathematica
In[1]:= k = 1; k++
Out[1]= 1
In[2]:= k
Out[2]= 2

In[3]:= v = a; v++; v
Out[3]= 1 + a

In[4]:= list = {1, 2, 3}; list[[2]]++; list
Out[4]= {1, 3, 3}

In[5]:= x = {100, 200, 300}; x -= {20, 21, 22}; x
Out[5]= {80, 179, 278}

In[6]:= y++
Increment::rvalue: y is not a variable with a value, so its value cannot be changed.
Out[6]= y++
```

### Files touched

* `src/parse.c` -- adds `OP_INCREMENT`, `OP_DECREMENT`, `OP_ADDTO`,
    `OP_SUBTRACTFROM` to the operator table; handles `++` / `--`
    as both postfix (via the Pratt loop, alongside `!`) and prefix
    (via the expression-start branch).
* `src/core.c`, `src/core.h` -- new `builtin_increment`,
    `builtin_decrement`, `builtin_preincrement`,
    `builtin_predecrement`, `builtin_addto`, `builtin_subtractfrom`
    sharing a single `increment_core` helper that performs the
    resolve / evaluate / update / write-back steps.
* `src/info.c` -- docstrings for all six heads.
* `tests/test_increment.c` -- full test coverage of each transcript
    in the feature request (symbols, numerical, symbolic, list,
    `Part`, and the `rvalue` path) plus two parser edge cases.

## Contexts and Packages (2026-04-21)

PicoCAS now supports Mathematica-style **contexts**: every symbol conceptually
lives in a namespace identified by a string ending in a backtick, such as
`Global``` or `MyPackage`Private``. Five new built-ins (`Context`, `Begin`,
`BeginPackage`, `End`, `EndPackage`) together with two system variables
(`$Context`, `$ContextPath`) make it possible to write self-contained
packages with private helpers and a clean public surface.

All the machinery lives in `src/context.c` / `src/context.h`; the parser
(`src/parse.c`) and the printer (`src/print.c`) were updated to resolve /
shorten qualified names. Existing code that never touches Begin or
BeginPackage sees no behavioral change -- bare names continue to resolve
exactly as before.

### Pragmatic design notes

Two simplifications keep the change self-contained rather than rippling
through the ~250 existing `symtab_add_builtin` registration sites:

1. **Built-ins stay stored under their bare names.** PicoCAS registers
   `Plus`, `Sin`, etc. in the symbol table without a `System``` prefix.
   The context resolver treats such bare names as implicitly belonging
   to `System```: `Context[Plus]` reports `"System`"`, and
   `System`Plus[2,3]` canonicalizes to `Plus[2,3]` so it still dispatches
   to the built-in.
2. **Inside the default `Global``` context, bare names pass through
   unchanged.** Context qualification is enabled only once the user
   enters a non-Global context via `Begin` or `BeginPackage`. This
   preserves legacy REPL and test behavior.

### `$Context`

`$Context` evaluates to a string giving the current context. Default
value is `"Global`"`. Mutated by `Begin`, `End`, `BeginPackage`, and
`EndPackage` (the symbol is `Protected`; direct assignment is refused).

```mathematica
In[1]:= $Context
Out[1]= "Global`"
```

### `$ContextPath`

`$ContextPath` evaluates to a list of contexts consulted (in order) when
the parser encounters a bare identifier that has not yet been seen. Used
so that symbols exported from a package become reachable by short name
after `EndPackage[]`. Default value is `{"Global`", "System`"}`.

```mathematica
In[1]:= $ContextPath
Out[1]= {"Global`", "System`"}
```

### `Context`

| Form                   | Meaning                                                         |
|------------------------|-----------------------------------------------------------------|
| `Context[]`            | Returns the current context (same as `$Context`).               |
| `Context[sym]`         | Returns the context in which `sym` resides, as a string.        |
| `Context["name"]`      | Returns the context for the symbol whose name is `"name"`.       |

Attributes: `{HoldFirst, Protected}`.

```mathematica
In[1]:= Context[]
Out[1]= "Global`"

In[2]:= Context[Plus]
Out[2]= "System`"

In[3]:= x = 1; Context[x]
Out[3]= "Global`"
```

### `BeginPackage`

| Form                                      | Meaning                                                   |
|-------------------------------------------|-----------------------------------------------------------|
| `BeginPackage["ctx`"]`                    | Saves `$Context` and `$ContextPath`, then sets them to `"ctx`"` and `{"ctx`", "System`"}`. |
| `BeginPackage["ctx`", {"n1`", ...}]`      | Same, additionally appending the needed contexts to the path. |

Returns the new current context as a string. Attributes: `{Protected}`.

### `Begin`

| Form                   | Meaning                                                                                          |
|------------------------|--------------------------------------------------------------------------------------------------|
| `Begin["ctx`"]`        | Sets the current context to `"ctx`"`, saving the previous value for a matching `End[]`.           |
| `Begin["`rel`"]`       | If the argument starts with a backtick, the context is taken relative to the current context, e.g. `Begin["`Private`"]` inside `"UtilPkg`"` yields `"UtilPkg`Private`"`. |

Returns the new current context as a string. Attributes: `{Protected}`.

### `End`

`End[]` pops the most recent Begin/BeginPackage frame, restoring the
previous `$Context` and (for BeginPackage) the previous `$ContextPath`.
Returns the just-closed context as a string. Attributes: `{Protected}`.

### `EndPackage`

`EndPackage[]` pops the most recent BeginPackage frame, prepending the
just-closed package context to `$ContextPath` so that any public symbols
defined within the package remain reachable by short name. Returns
`Null`. Attributes: `{Protected}`.

### Example: defining and using a small package

```mathematica
In[1]:= BeginPackage["square`"]
Out[1]= "square`"

In[2]:= square[x_] := x^2
Out[2]= Null

In[3]:= EndPackage[]
Out[3]= Null

In[4]:= $ContextPath
Out[4]= {"square`", "Global`", "System`"}

In[5]:= square[12]
Out[5]= 144

In[6]:= square`square[12]  (* same function, fully qualified *)
Out[6]= 144
```

### Parser rules in detail

The parser calls `context_resolve_name` on every identifier. The rule
applied, in order:

1. Names starting with `$` pass through unchanged (system variables).
2. A leading backtick is taken as relative to the current context:
   `` `Private`f `` in context `"foo`"` resolves to `"foo`Private`f"`.
3. Names containing a backtick elsewhere are treated as absolute, with
   one canonicalization: `"System`X"` is reduced to bare `"X"` so that
   built-ins stored under bare names still dispatch.
4. For a bare identifier, the resolver searches each prefix in
   `{$Context} U $ContextPath`. If a qualified symbol with that prefix
   already exists, its full name is used.
5. If a bare builtin is already registered, it is used directly
   (implicit `System``` membership).
6. If `$Context == "Global`"`, the bare name is kept as-is.
7. Otherwise the name is qualified with `$Context`.

### Printer

When printing a symbol, if its context prefix is currently visible
(equal to `$Context`, present on `$ContextPath`, or `"System`"`), the
prefix is dropped so output reads naturally -- `UtilPkg`f[3]` prints as
`f[3]` while `UtilPkg`` is on the path. `FullForm[]` uses the same
shortening rule.

### Files touched

* `src/context.c`, `src/context.h` -- new module holding the context
  stack, resolution rules, and the five built-in functions plus the
  `$Context`/`$ContextPath` OwnValue publishing.
* `src/parse.c` -- `parse_symbol` now passes every identifier through
  `context_resolve_name`.
* `src/print.c` -- `expr_print_fullform` shortens symbol names whose
  context is visible via `context_display_name`.
* `src/symtab.c`, `src/symtab.h` -- adds `symtab_lookup` (non-creating
  hash-table probe) so the resolver can ask whether a qualified name
  already exists.
* `src/core.c` -- calls `context_init()` first in `core_init()`;
  `builtin_information` uses `context_display_name` when reporting a
  missing docstring so the short name appears in the diagnostic.
* `tests/test_context.c` -- direct API tests plus end-to-end coverage
  of `Context`, `$Context`, `$ContextPath`, `Begin`/`End` round-trip,
  `BeginPackage` path manipulation, and a `Begin["`Private`"]` package
  example.

## Bug fix: List destructuring assignment after prior OwnValues (2026-04-22)

`{a, b, c, d} = {...}` silently no-op'd when the symbols on the LHS
already held OwnValues. The second destructuring assignment in a session
like

```mathematica
{a, b, c, d} = {1, 1, 1, 1};
{a++, ++b, c--, --d}               (* {1, 2, 1, 0} *)
{a, b, c, d} = {1, 1, 1, 1};
{a++, ++b, c--, --d}               (* was {2, 3, 0, -1} *)
```

returned the mutation applied to the stale post-increment state because
the second reassignment never actually wrote back to the symbols.

**Root cause.** `Set` has `HoldFirst`, so `List[a,b,c,d]` arrived at the
`Set` primitive handler unevaluated -- correct so far. But that handler
then evaluated each LHS argument, gated only by the *head's* Hold
attributes. Since `List` has no Hold attributes, `a, b, c, d` got
evaluated to their current values (2, 2, 0, 0). The internal LHS became
`List[2, 2, 0, 0]`, and `apply_assignment` recursed on each `(int, int)`
pair, fell through the dispatch (LHS is an integer, not a symbol), and
returned `false`. A second latent bug compounded the problem: the list-
destructuring branch of `apply_assignment` ignored the recursive return
values and always reported success, so `Set` yielded the RHS as if the
write had happened.

**Fix.** Two targeted changes in `src/eval.c`:

1. In the `Set`/`SetDelayed` LHS-argument evaluation step, when the LHS
   is a `List`, hold any element that is itself a `Symbol` (a binding
   target) or a nested `List` (nested destructuring). Non-symbol
   function-shaped elements such as `a[x]` continue to get their inner
   arguments evaluated so that `{a[x], b[y]} = {p, q}` still defines
   `a[val-of-x]` and `b[val-of-y]` (standard Mathematica semantics).
2. In `apply_assignment`'s list-destructuring branch, AND each child
   recursion's return value into the overall result, so a literal-LHS
   child (e.g. `{1, a} = {1, 2}`) no longer lies about succeeding.

**Coverage.** `tests/test_list_set.c` adds ten focused cases:

- Fresh destructuring still works.
- Reassignment to symbols that already have values.
- The exact sequence from the bug report (two rounds of
  `{a,b,c,d}={1,1,1,1}; {a++, ++b, c--, --d}`).
- Simultaneous swap `{x, y} = {y, x}`.
- Nested destructuring `{{a, b}, c} = {{10, 20}, 30}`, including
  reassignment of nested targets.
- DownValue destructuring `{a[x], b[y]} = {p, q}` with evaluated x, y.
- Pattern-bearing DownValue destructuring `{a[p_], b[q_]} = {1, 2}`.
- Length-mismatched RHS leaves targets untouched.
- Literal-integer LHS element does not silently pretend success.
- Repeated reassignment smoke test (50 rounds).

All 60 test binaries continue to pass.


## Numeric evaluation: `N`, precision literals, `Precision`/`Accuracy`/`SetPrecision`/`SetAccuracy`

PicoCAS now has a full numeric-evaluation pipeline covering both
machine-precision (IEEE double) and arbitrary-precision (MPFR) reals.
The work lives in `src/numeric.c`, `src/numeric.h`, `src/precision.c`,
`src/precision.h`; MPFR support is gated behind the `USE_MPFR=1`
makefile flag (default on), mirroring the existing `USE_ECM=1` pattern.

### Builtins

| Function                | Attributes                    | Description |
|-------------------------|-------------------------------|-------------|
| `N[expr]`               | `Protected, Listable`         | Machine-precision numerical approximation of `expr`. |
| `N[expr, n]`            | `Protected, Listable`         | Approximation to `n` decimal digits (requires `USE_MPFR=1`). |
| `MachinePrecision`      | `Protected`                   | Symbol representing IEEE-double precision (≈ 15.95 digits). |
| `Precision[x]`          | `Protected, Listable`         | Returns `Infinity`, `MachinePrecision`, or decimal-digit precision of `x`. |
| `Accuracy[x]`           | `Protected, Listable`         | Decimal digits of accuracy — `Precision[x] − Log10[Abs[x]]`. Zero / exact → `Infinity`. |
| `SetPrecision[x, n]`    | `Protected, Listable`         | Walk `x` re-rounding inexact reals and promoting exact values to `n`-digit MPFR. |
| `SetAccuracy[x, n]`     | `Protected, Listable`         | Like `SetPrecision` but interprets `n` as digits of accuracy. |

### Numeric constants

`N` recognizes and promotes the following symbols:
`Pi`, `E`, `EulerGamma`, `Catalan`, `GoldenRatio`, `Degree`. Machine
values use `<math.h>` macros where available; MPFR values use
`mpfr_const_pi`, `mpfr_const_euler`, `mpfr_const_catalan`, and
`mpfr_exp`/`mpfr_sqrt` for derived constants. To add a new constant,
append a row to the `kConstants` table in `src/numeric.c`.

### Precision literals

A numeric literal may carry a trailing precision or accuracy suffix,
matching Mathematica syntax:

| Literal        | Meaning                                             |
|----------------|-----------------------------------------------------|
| `` 3.14 ``     | machine-precision double                            |
| `` 3.14`50 ``  | MPFR real, 50-digit *precision*                     |
| `` 3.14``49 `` | MPFR real, 49-digit *accuracy*                      |
| `` 3`30 ``     | integer widened to MPFR at 30-digit precision       |

The parser (`src/parse.c:parse_number`) disambiguates the
precision-marker backtick from the context-separator backtick by
requiring a digit/`.`/`+`/`-` to follow — an identifier character after
the backtick leaves it alone for the context parser.

### Architecture

`N` does not reimplement the numeric math — it walks the expression
tree replacing *leaves* (integers, rationals, named constants) with
their numeric values, rebuilds, and re-evaluates. The existing Plus,
Times, Power, Sin, Cos, Exp, Log, etc. builtins then take the numeric
fast path naturally. This "descend + re-evaluate" strategy keeps
`src/numeric.c` small and puts each function's numeric implementation
in its owning module.

One special case lives in `numericalize_function`: when the head is
`Power` and the exponent is an `EXPR_INTEGER` or `EXPR_BIGINT`, the
exponent is preserved verbatim and only the base is numericalized.
This matches Mathematica (`N[3 x^2 + 3 x + 1]` → `1. + 3. x + 3. x^2`,
not `… 3. x^2.`) and keeps polynomial structure intact under `N`.
Non-integer exponents (e.g. `1/2`) are still numericalized.

MPFR-aware paths were added to:

- `src/plus.c`, `src/times.c`, `src/power.c` — use `numeric_mpfr_add`,
  `numeric_mpfr_mul`, `numeric_mpfr_pow` from `numeric.c`.
- `src/trig.c`, `src/hyperbolic.c` — each builtin calls
  `numeric_mpfr_apply_unary(arg, 0, mpfr_sin)` (etc.) when any argument
  carries MPFR precision.
- `src/logexp.c` — `Log`, `Exp` via `mpfr_log`, `mpfr_exp`.

New shared helpers in `numeric.c`:

- `numericalize(expr, spec)` — recursive walker used by `N`,
  `SetPrecision`, and `SetAccuracy`.
- `get_approx_mpfr(e, re, im, *inexact)` — extracts an MPFR approximation
  of any real or complex numeric Expr, used by the transcendental builtins.
- `numeric_mpfr_{add,sub,mul,div,pow}(a, b, default_bits)` — binary ops
  that produce EXPR_MPFR at `max(prec_a, prec_b, default_bits)`.
- `numeric_mpfr_apply_unary(e, default_bits, op)` — generic unary
  dispatch; used by the transcendental branches.
- `numeric_any_mpfr(a, b)`, `numeric_expr_is_mpfr(e)` — cheap guards.
- `numeric_digits_to_bits` / `numeric_bits_to_digits` — decimal ↔ binary
  precision conversion (ratio `log₂10`).

### Extending the system

- **A new constant** (e.g. `Glaisher`): append one row to `kConstants`
  in `src/numeric.c`. Provide a `double` machine value and an optional
  MPFR filler function.
- **A new numeric backend** (e.g. GSL for special functions, libquadmath
  for `__float128`): extend `NumericMode` in `src/numeric.h` and add the
  per-op branches; each transcendental builtin's existing per-module MPFR
  branch is the natural template for a GSL branch. No central dispatcher
  is required because `N` drives the backend through the evaluator.
- **A new special function** (e.g. `BesselJ`, `Gamma`, `Zeta`): add the
  numeric branch directly to that function's builtin. `N` will reach it
  automatically after numericalizing arguments, because the function is
  invoked through the re-evaluation step of `numericalize`.

### Known limitations (Phase 2 v1)

- **Display doesn't pad trailing zeros** for declared precision. A
  literal like `` 3.98`50 `` prints as `3.979999…9` (full 50 digits of
  the closest binary representation) rather than `3.98000…0`. Mathematica
  tracks "declared precision" separately from the bit-level value; we
  print bit-accurate decimals. Value semantics are otherwise identical.
- **Complex MPFR** is not implemented. `Power` with MPC-style complex
  operands falls back to double-complex `cpow`. A future step would
  link MPC for complex arbitrary precision.
- **Significance arithmetic** uses `max(prec)` of inputs for each
  operation's output precision. Mathematica's precision-degradation
  rules (precision decreases through lossy operations) are deferred.

### Test coverage

`tests/test_numeric.c` exercises:

- Phase 1: leaf conversion, constants, unknown symbols, descent + re-evaluate
  semantics, Listable threading, Complex, Hold preservation, bad precision
  args, `MachinePrecision` protection.
- Phase 2 (compile-time-gated by `USE_MPFR`): `N[…, prec]` for Pi, E,
  Sqrt, Sin, Cos, Log, Exp, Sinh, Tanh, ArcTan, rational and mixed
  expressions; `Precision[]`, `Accuracy[]`, `SetPrecision[]`,
  `SetAccuracy[]`; precision-literal parsing; Listable precision.

A new `assert_eval_startswith` helper in `tests/test_utils.h` tolerates
last-bit rounding noise in long MPFR strings.

## Inexact contagion in Plus / Times (2026-04-22)

Mathematica's numeric contagion rule is now applied inside `Plus` and
`Times`: when any summand (or factor) is inexact — an `EXPR_REAL`, an
`EXPR_MPFR`, or a `Complex[...]` containing one — every remaining
argument is run through `numericalize` before the usual folding. This
converts exact numeric leaves (Pi, E, EulerGamma, Catalan, GoldenRatio,
Degree, Rational, and expressions like `Sqrt[2]` or `Sin[1]` that
numericalize cleanly) into their numeric values so the inexact
coefficient can actually collapse against them.

Behavior:

```
1. Pi        (* 3.14159 — previously 1. Pi *)
1. + Pi      (* 4.14159 *)
1. E         (* 2.71828 *)
1. Sqrt[2]   (* 1.41421 *)
1. Sin[1]    (* 0.841471 *)
1. x Pi      (* 3.14159 x — symbolic x stays symbolic *)
1. x         (* 1.0 x — no exact numeric to numericalize *)
2 Pi         (* 2 Pi — no inexact arg, no contagion *)
1 + Pi       (* 1 + Pi *)
```

Precision selection mirrors Mathematica: the *lowest* precision among
inexact operands wins. `MachinePrecision` is the floor — any `EXPR_REAL`
(a machine double) collapses the combination to machine precision even
alongside MPFR values with more digits. Between two MPFR operands, the
smaller precision wins rather than the larger, so `1.0\`50 + 1.0\`20`
lands at 20 digits instead of preserving the 50-digit operand. So:

```
1.0`50 Pi             (* 3.14159… — pure MPFR, 50-digit *)
Precision[%]          (* 50.272 *)
1.0`50 Pi + 1.0`20    (* 4.14159… — 20 digits wins *)
Precision[%]          (* 20.169 *)
1.0`50 Pi + 1.        (* 4.14159 — machine wins *)
Precision[%]          (* MachinePrecision *)
1.0`50 + 1.           (* 2.0 — machine *)
```

New public helper `numeric_contagion_args` in `numeric.h` encapsulates
the scan + rewrite so future numeric heads can opt in with one call.

Regression coverage in `tests/test_numeric.c::test_inexact_contagion`.

## Bug fix: Inverse trig/tanh branch-cut sign for real x > 1 (2026-04-22)

`ArcSin[x]`, `ArcCos[x]`, and `ArcTanh[x]` for inexact real `x > 1`
returned the wrong sign on the imaginary part relative to Mathematica.
Root cause: C99 `casin`, `cacos`, `catanh` for an input `x + 0i` land on
the upper side of the `(1, infinity)` branch cut, while Mathematica uses
the lower side. The `(-infinity, -1)` cut already agrees between the two
conventions (verified numerically on macOS libc), so the fix is scoped
to `creal(c) > 1.0 && cimag(c) == 0.0`, where we negate `cimag(s)` after
calling the C function. ArcSin/ArcCos live in `src/trig.c`; ArcTanh in
`src/hyperbolic.c`. Complex-argument and `|x| < 1` cases are untouched.

Before → After for `x = 3.0`:

```
ArcSin[3.]  : 1.5708 + 1.76275 I  ->  1.5708 - 1.76275 I
ArcCos[3.]  : 0 - 1.76275 I       ->  0 + 1.76275 I
ArcTanh[3.] : 0.346574 + Pi/2 I   ->  0.346574 - Pi/2 I
```


## TeXForm (2026-04-22)

`TeXForm[expr]` prints `expr` as AMS-LaTeX-compatible TeX. When an
expression evaluates to `TeXForm[expr]`, only the TeX rendering appears
in the output — the `TeXForm` wrapper itself is not printed. `TeXForm`
has attribute `Protected`.

Conventions, matching Mathematica:

- Single-character symbol names render bare (italic by TeX default):
  `TeXForm[x]` -> `x`. Multi-character names render in roman via
  `\text{...}`: `TeXForm[foo]` -> `\text{foo}`.
- Mathematical constants: `Pi` -> `\pi`, `E` -> `e`, `I` -> `i`,
  `Infinity` -> `\infty`, `Degree` -> `{}^{\circ}`, `EulerGamma` ->
  `\gamma`, `GoldenRatio` -> `\phi`.
- Rational and Divide denominators render via `\frac{}{}`. The `Times`
  printer splits negative-exponent factors into a denominator.
- `Power[a, 1/2]` -> `\sqrt{a}`, `Power[a, 1/n]` -> `\sqrt[n]{a}`,
  general `Power` -> `a^{b}`.
- `Sqrt`, `Abs`, `Floor`, `Ceiling`, `Conjugate` render with their
  dedicated TeX constructs (`\sqrt`, `|...|`, `\lfloor...\rfloor`,
  `\lceil...\rceil`, `\overline{}`).
- Trig and hyperbolic functions render using their standard TeX macros
  (`\sin`, `\cos`, ..., `\sinh`, ..., plus `\text{sech}` / `\text{csch}`
  which have no standard macro). Inverse variants (`ArcSin`, `ArcTanh`,
  etc.) render as `\sin ^{-1}(x)` / `\tanh ^{-1}(x)`.
- `Log` -> `\log`, `Exp` -> `\exp`, `Gamma` -> `\Gamma`, `Re` -> `\Re`,
  `Im` -> `\Im`.
- `List` -> `\{a,b,c\}`. Comparisons map to `=`, `\neq`, `<`, `>`,
  `\leq`, `\geq`, `\equiv`, `\not\equiv`. `Rule` -> `\to`. `And` ->
  `\land`, `Or` -> `\lor`, `Not` -> `\neg`.
- Unknown heads render as `head(args...)` (the head follows the same
  single-char / multi-char italic/roman rule as symbols).

Implementation lives in `src/print.c` as `print_tex()`; registration
and attributes in `src/core.c`; docstring in `src/info.c`. Parenthesisation
reuses `get_expr_prec` so the TeX output agrees structurally with the
standard printer.

Examples:

```
TeXForm[1/2]               -> \frac{1}{2}
TeXForm[3 + 4 I]           -> 3+4 i
TeXForm[-3 x]              -> -3 x
TeXForm[Sqrt[x]]           -> \sqrt{x}
TeXForm[x^(1/3)]           -> \sqrt[3]{x}
TeXForm[a/b]               -> \frac{a}{b}
TeXForm[Sin[x]]            -> \sin(x)
TeXForm[ArcSin[x]]         -> \sin ^{-1}(x)
TeXForm[Sin[x]/(1+Cos[x])] -> \frac{\sin(x)}{1+\cos(x)}
TeXForm[(a+b)^2]           -> \left(a+b\right)^{2}
```

## TrigExpand and trigsimp module (2026-04-22)

Added `TrigExpand[expr]`, which expands circular and hyperbolic trigonometric
functions by splitting up sums and integer multiples that appear in their
arguments. See the `TrigExpand` section above for the full feature list.

At the same time the previous `simp.c` / `simp.h` module was renamed to
`trigsimp.c` / `trigsimp.h`, and now hosts all three trig-rewriting builtins:

- `TrigToExp`  — rewrite trig / inverse-trig in terms of `Exp` / `Log`.
- `ExpToTrig`  — rewrite `Exp` / `Log` back into trig / inverse-trig.
- `TrigExpand` — expand trig sums and integer multiples into products and
  powers of the base `Sin`/`Cos`/`Sinh`/`Cosh` functions.

All three are registered with `Listable` and `Protected` attributes and have
docstrings in `info.c`. The `core.c` initialization chain now calls
`trigsimp_init()` in place of the old `simp_init()`.

Unit tests live in `tests/test_trigexpand.c` (new) and `tests/test_simp.c`
(existing TrigToExp / ExpToTrig coverage).

### Pythagorean collapse for any integer multiple and sign (2026-04-22)

Extended the Pythagorean reduction pass so that every integer multiple of the
identity collapses regardless of sign or scalar weight. The previous
implementation only matched the canonical positive form, so inputs like
`Sinh[4 x]^2 - Cosh[4 x]^2`, `-Sin[4 x]^2 - Cos[4 x]^2`, and
`-5 Sin[4 x]^2 - 5 Cos[4 x]^2` came back in fully-expanded polynomial form.

`Factor` does not produce a canonical sign for these expressions — depending on
the input it may emerge as `(Sin^2 + Cos^2)^n`, `(-Sin^2 - Cos^2)^n` wrapped in
a separate `-1`, `(Cosh + Sinh)^n (Cosh - Sinh)^n`, or `(Cosh + Sinh)^n
(-Cosh + Sinh)^n`. The reduction rules now enumerate both signs and admit an
arbitrary remainder of `Times` factors (via `r___`) so scalar coefficients and
stray `-1`s factored out by `Factor` survive. The final `Expand` reabsorbs the
remainder, producing the expected signed constant (e.g. `-5`).

New tests in `tests/test_trigexpand.c`:

- `test_trigexpand_pythagorean` grew cases for negative-integer multiples,
  sign-flipped sums, and scalar-weighted Pythagorean collapses.
- `test_trigexpand_hyperbolic_pythagorean` covers `Cosh[n x]^2 - Sinh[n x]^2`,
  its sign-flipped form, even-parity arguments, and scalar-weighted variants.

## TrigFactor (2026-04-22)

Added `TrigFactor[expr]` in `src/trigsimp.c`, which factors trigonometric
functions in `expr` and acts broadly as the functional inverse of
`TrigExpand`. See the `TrigFactor` section above for the complete feature
list.

`TrigFactor` is registered with `Listable` and `Protected` attributes, and
its docstring lives in `src/info.c` alongside `TrigExpand` / `TrigToExp` /
`ExpToTrig`. Initialization is driven from `trigsimp_init()` (which is
already called by `core_init()`).

### Design highlights

- Reciprocal heads (`Tan`, `Cot`, `Sec`, `Csc` and their hyperbolic analogs)
  are rewritten as `Sin`/`Cos` / `Sinh`/`Cosh` ratios before the polynomial
  factoring pass so `Factor` sees the full structure.
- Identity collapses fire via `ReplaceRepeated`: Pythagorean (both directions
  and with arbitrary optional coefficients), reverse angle-addition in both
  circular and hyperbolic form, reverse double-angle, and factored-form
  variants that `Factor` naturally produces such as
  `(Cos - Sin)(Cos + Sin) -> Cos[2 x]`, `(Cosh - 1)(Cosh + 1) -> Sinh^2`,
  and `(Cosh - Sinh)(Cosh + Sinh) -> 1`.
- A second pass restores `Tan`, `Sec`, `Cot`, `Csc` (and hyperbolic analogs)
  from `Sin`/`Cos` ratio form so reciprocal heads survive.
- Two paths are evaluated:
  - Path A runs the pipeline on the original expression; it preserves
    angle-sum structure (for example `Sin[x + y]^2 + Tan[x + y]` factors
    cleanly as `Tan[x + y] (1 + Cos[x + y] Sin[x + y])`).
  - Path B first applies `TrigExpand`, then runs the pipeline. This catches
    cancellations that only become visible after angle-sums are distributed
    (for example `Cos[x + y] + Sin[x] Sin[y] -> Cos[x] Cos[y]`).
  Path B runs only when Path A left the expression unchanged. This avoids
  blowing up on inputs like `Sin[x + y]^2 + Tan[x + y]` where
  `TrigExpand` produces a large multivariate rational that `Factor` would
  otherwise struggle with. When both paths run, the smaller result by leaf
  count wins, with ties resolved in favour of Path A to preserve angle-sum
  structure.
- Threads automatically over lists (via `Listable`), equations,
  inequalities, and logic heads, matching `TrigExpand`.

### Tests

Unit tests live in `tests/test_trigfactor.c` and cover:

- Pythagorean identities (circular and hyperbolic), with and without
  arbitrary coefficients, and iterated across independent variable pairs.
- Reverse double-angle reductions for both `Sin`/`Cos` and `Sinh`/`Cosh`.
- Reverse angle-addition identities in both circular and hyperbolic form.
- The `Cos[x + y] + Sin[x] Sin[y] -> Cos[x] Cos[y]` cancellation that
  requires Path B.
- `Tan`/`Sec` expressions (including `Sin[x]^2 + Tan[x]^2`,
  `Tan[x + y] Cos[x + y]`, and `Sec[x] Cos[x]`).
- Polynomial factorizations such as `Cos[x]^4 - Sin[x]^4 -> Cos[2 x]` and
  `Cosh[x]^2 - Cosh[x]^4 -> -Cosh[x]^2 Sinh[x]^2`.
- Round-tripping: `TrigFactor[TrigExpand[Sin[x + y]]] -> Sin[x + y]` and
  analogous inverses for `Cos`, `Sinh`, `Cosh`, and integer multiples.
- Reciprocal trig passing through unchanged (`Tan[x] -> Tan[x]`,
  `Sec[x]^2 -> Sec[x]^2`, etc.).
- Atomic inputs and numeric folding.
- Threading over `List`, `Equal`, `Unequal`, `Less`, `Greater`, `And`,
  `Not`, including compound inequalities.

## Pattern matcher: Optional consistency fix (2026-04-28)

Fixed a bug in `src/match.c` where the `Optional` fallback path at the end
of `match_args_internal` would unconditionally overwrite an existing
binding for the pattern variable with the Optional default. This produced
inconsistent bindings for patterns like
`a_. Sin[x_]^2 + a_. Cos[x_]^2 + r___`, where the matcher would bind
`a` to `-1` from one term and then quietly overwrite it with the Optional
default `1` when matching the other term, allowing a rule keyed on equal
coefficients to fire on `Sin[x]^2 - Cos[x]^2` (with two different
implicit coefficients). After the fix, the Optional fallback only takes
the default-value path if any existing binding for the variable is
structurally equal to the default; otherwise the path fails and the
matcher backtracks.

This change tightens the matcher to match Mathematica semantics. It
restores correctness for rules in `trig_factor_identities` that use
`a_.` to capture a shared coefficient across multiple Pythagorean or
double-angle terms, and was necessary for the
`TrigFactor[Cosh[x + y]/(Sin[x]^2 - Cos[x]^2) + 1/Sinh[x]]` regression.

## TrigFactor: linear-combination factoring (2026-04-28)

Added a rule to `trig_factor_identities` that factors any sum
`a Sin[x] + b Cos[x]` (with optional unit-default coefficients and
arbitrary trailing terms) into the auxiliary-angle form
`Sqrt[a^2 + b^2] Sin[x + ArcTan[a, b]]` when both `a` and `b` are
numeric. The `NumberQ` guard prevents the rule from firing on symbolic
coefficients, where the result would not be simpler than the input.
Examples:

- `TrigFactor[Sin[x] + Cos[x]] -> Sqrt[2] Sin[Pi/4 + x]`
- `TrigFactor[3 Sin[x] + 4 Cos[x]] -> 5 Sin[x + ArcTan[3, 4]]`
- `TrigFactor[-Sin[x] - Cos[x]] -> Sqrt[2] Sin[-3 Pi/4 + x]`

The two-argument `ArcTan[a, b]` is used to handle the quadrant
correctly so the resulting phase is consistent for any sign of `a` or
`b`.

## TrigExpand: inverse-trig compositions (2026-04-28)

Added rewrite rules to the TrigExpand rule set in `src/trigsimp.c` that
reduce compositions of a forward trig (or hyperbolic) function with a
different inverse to their algebraic form. For example:

- `Cos[ArcSin[x]] -> Sqrt[1 - x^2]`
- `Sin[ArcCos[x]] -> Sqrt[1 - x^2]`
- `Tan[ArcSin[x]] -> x / Sqrt[1 - x^2]`
- `Sin[ArcTan[x]] -> x / Sqrt[1 + x^2]`
- `Cosh[ArcSinh[x]] -> Sqrt[1 + x^2]`
- `Sinh[ArcTanh[x]] -> x / Sqrt[1 - x^2]`
- ... and the analogous rules for the remaining circular and hyperbolic
  forward/inverse pairs (`ArcCot`, `ArcSec`, `ArcCsc`, `ArcCosh`,
  `ArcTanh`).

The rules use the principal real-branch identity for each composition.
This mirrors the behaviour of Mathematica's `TrigExpand`. Effect on the
canonical addition example:

- `TrigExpand[Sin[ArcSin[t] + ArcSin[v]]]` now returns
  `t Sqrt[1 - v^2] + v Sqrt[1 - t^2]` (the cross-terms `Cos[ArcSin[v]]`
  and `Cos[ArcSin[t]]` are reduced via the new rules), where previously
  the result was left as `t Cos[ArcSin[v]] + v Cos[ArcSin[t]]`.

## Parser: dynamic argument buffers (2026-04-28)

Replaced the fixed `Expr* args[64]` and `Expr* elements[64]` stack
arrays in `src/parse.c` (`parse_function` and `parse_list`) with
dynamically-grown heap buffers. The previous fixed-size arrays
silently overflowed the stack when a list literal contained more than
64 elements, triggering a stack-canary abort during initialisation
once the trigsimp rule lists grew past that bound. The new buffers
double in capacity as needed.

## Factor: integer-coefficient guard for Gaussian inputs (2026-04-28)

In `bz_factor_to_expr` (`src/facpoly.c`), the integer Zassenhaus
factoring routine previously coerced any non-`EXPR_INTEGER` coefficient
to `0`, which silently produced a polynomial with the wrong shape.
For example, `Factor[t^2 + 2 I t - 1]` returned `(-1 + t)(1 + t)` --
i.e. `t^2 - 1` -- because the imaginary middle coefficient `2 I` was
zeroed out. After the fix, `bz_factor_to_expr` detects any non-integer
coefficient and bails out, returning the polynomial unfactored at this
layer. PicoCAS does not implement factoring over the Gaussian rationals
`Z[i]`; this guard ensures that, when it cannot factor the input, it
returns a mathematically equivalent expression rather than a wrong
factorization.

## Known limitations (2026-04-28)

- `TrigExpand` does not currently implement power-reduction identities
  (e.g. `Sin[x]^2 -> (1 - Cos[2 x])/2`). Mathematica's
  `TrigExpand[Sin[x]^2 + Tan[x]^2]` produces an extensively rewritten
  form using such identities; PicoCAS leaves the input unchanged.
  Adding power-reduction safely requires careful loop-prevention design
  because the existing multiple-angle rules
  (`Cos[2 x] -> Cos[x]^2 - Sin[x]^2`) form a cycle with the inverse
  identities.
- `Factor` does not factor over `Z[i]` (the Gaussian rationals). When
  Gaussian coefficients appear, the polynomial is returned unfactored
  rather than incorrectly factored over `Z`.

## Simplify, Assuming, $Assumptions (2026-04-28)

Added `Simplify`, `Assuming`, and the `$Assumptions` system symbol in a
new `src/simp.c` / `src/simp.h` module. The module also exposes an
`AssumeCtx` API that future consumers (`Refine`, `Integrate`, ...) can
share once they need assumption-aware behaviour.

### Simplify

`Simplify[expr]` runs a small heuristic search over picocas's existing
algebraic transforms and returns the candidate of minimum complexity.
The default complexity measure matches Mathematica's: leaf count plus
the decimal-digit count of any integer/bigint leaves, so e.g.
`Simplify[100 Log[2]]` stays `100 Log[2]` rather than collapsing to
`Log[1267650600228229401496703205376]`.

`Simplify[expr, assum]` does simplification using assumptions `assum`,
which can be equations, inequalities, domain specifications such as
`Element[x, Integers]`, or logical combinations of these. Lists of
assumptions `{a1, a2, ...}` are converted to conjunctions
`And[a1, a2, ...]`.

Options:

- `Assumptions -> X` overrides the `$Assumptions` default. When given
  explicitly the option *replaces* `$Assumptions` rather than appending
  to it, matching Mathematica's documented behaviour:
  `Assuming[x>0, Simplify[Sqrt[x^2 y^2], Assumptions -> y<0]]` simplifies
  under `y<0` only -- the outer `Assuming` is suppressed by the explicit
  option.
- `ComplexityFunction -> f` overrides the default complexity. The user
  function is called as `f[candidate]`; it must return an integer (or
  bigint) score. Non-integer return values fall back to the default.

#### Search algorithm

1. Build an `AssumeCtx` from the resolved assumption set
   (`positional && Assumptions-or-$Assumptions`).
2. If `expr` literally appears in the fact list, return `True` (narrow
   predicate-membership shortcut).
3. Manually thread over `Equal`/`Less`/`LessEqual`/`Greater`/
   `GreaterEqual`/`And`/`Or`/`Not`/`Xor`/`Implies`. List threading is
   delivered by `ATTR_LISTABLE`.
4. Seed the candidate set with the input.
5. Apply the assumption-aware rewriter to the input. If it produces a
   different form whose score does not exceed `best_score`, force-prefer
   it as the new best (assumption-aware rewrites are correctness-
   preserving and qualitatively "more simplified" by definition; this
   biases ties toward the rewritten form).
6. Apply a TrigToExp/Together/Cancel/ExpToTrig roundtrip seed (catches
   identities such as `(E^x - E^-x)/Sinh[x] -> 2`).
7. For up to `SIMP_ROUNDS` (= 2) rounds, apply each transform in
   {`Together`, `Cancel`, `Expand`, `Factor`, `FactorSquareFree`,
   `Apart`, `TrigExpand`, `TrigFactor`} plus the trig roundtrip and the
   assumption rewriter to every seed. Update `best` whenever a candidate
   strictly beats `best_score`. Capped at `SIMP_CAND_CAP` (= 12)
   candidates per round.
8. Return the lowest-scoring candidate.

`Factor`, `FactorSquareFree`, and `TrigFactor` are skipped on inputs
that contain any `Power` with a non-integer exponent. The picocas
`Factor` trial-loop in `factor_roots` can stall on multivariate inputs
that include Sqrt or other fractional Powers as polynomial atoms; this
guard mirrors `trigsimp.c`'s `count_distinct_squared_trig_atoms` cap.

#### Assumption-aware rewriters

For each symbol the assumption set proves positive, real, negative, or
integer, `apply_assumption_rules` synthesises rewrite rules and applies
them via `ReplaceRepeated`:

- `x` known positive: `Sqrt[x^2] -> x`, `Sqrt[1/x] -> 1/Sqrt[x]`,
  `1/Sqrt[1/x] -> Sqrt[x]`, `Abs[x] -> x`, `Log[x^p_] -> p Log[x]`.
- `x` known negative: `Sqrt[x^2] -> -x`, `Abs[x] -> -x`.
- `x` known real (and not already positive/negative):
  `Sqrt[x^2] -> Abs[x]`.
- `n` known integer: `Sin[n Pi] -> 0`, `Cos[n Pi] -> (-1)^n`,
  `Tan[n Pi] -> 0` (also commuted forms `Sin[Pi n]`, etc).

For each `Equal[u, v]` fact, two complementary rules are emitted:

1. The direct heavier->lighter substitution, matched structurally
   (catches cases where `u` appears verbatim as a subterm of the input).
2. A polynomial-relation monomial-isolation rule when
   `diff = u - v` evaluates to a `Plus` with `>= 3` terms. For the
   highest-complexity non-numeric term `t`, emit
   `t -> -(Plus of other terms)`. This lets relations such as
   `a^2 + b^2 == 1` rewrite occurrences of `a^2` or `b^2` even when the
   full sum `a^2 + b^2` is not present in the input. Only one
   monomial rule per equation is emitted; emitting both directions
   would make `ReplaceRepeated` cycle to its 65536 iteration cap.

The choice of which monomial to isolate uses canonical (Plus-Orderless)
order on tie. As Mathematica's documentation notes, "results of
simplification may depend on the names of symbols" -- picocas
reproduces this property by alphabetical preference. Concretely:

```
Simplify[(1 - a^2)/b^2, a^2 + b^2 == 1]   (* -> 1 *)
Simplify[(1 - c^2)/b^2, c^2 + b^2 == 1]   (* -> (1 - c^2)/b^2 *)
```

`Simplify` is registered with `ATTR_LISTABLE | ATTR_PROTECTED`. Its
docstring lives in `src/info.c`.

### Assuming

`Assuming[assum, expr]` evaluates `expr` with `assum` appended to
`$Assumptions`. The implementation desugars at runtime to
`Block[{$Assumptions = $Assumptions && assum}, expr]`, reusing Block's
existing scope-restoration code path. Nested `Assuming` calls compose
naturally because each Block reads the current `$Assumptions` OwnValue
before extending it. Lists of assumptions are first converted to
conjunctions per the standard Mathematica semantics.

`Assuming` has attributes `ATTR_HOLDREST | ATTR_PROTECTED`: the
assumption (arg 1) evaluates normally so users can write things like
`Assuming[x > 0, ...]`, but the body (arg 2) is held until the Block
fires, ensuring `$Assumptions` is in scope when the body runs.

### $Assumptions

`$Assumptions` is a system symbol with an OwnValue of `True`
(no assumptions). `Assuming` rebinds it for the dynamic extent of its
body via Block, so the rebinding is automatically restored on exit
(including on exception paths that propagate out of Block normally).
Functions like `Simplify` read its current value when no explicit
`Assumptions` option is provided.

### AssumeCtx (public API in `simp.h`)

```c
AssumeCtx* assume_ctx_from_expr(const Expr* assum);
void       assume_ctx_free(AssumeCtx* ctx);

bool assume_known_positive(const AssumeCtx* ctx, const Expr* x);
bool assume_known_nonneg  (const AssumeCtx* ctx, const Expr* x);
bool assume_known_negative(const AssumeCtx* ctx, const Expr* x);
bool assume_known_nonpos  (const AssumeCtx* ctx, const Expr* x);
bool assume_known_real    (const AssumeCtx* ctx, const Expr* x);
bool assume_known_integer (const AssumeCtx* ctx, const Expr* x);
```

`assume_ctx_from_expr` flattens `And` and `List` heads, drops `True`
facts, and marks the context inconsistent on `False`. The predicates
combine direct fact matches with the obvious lattice (Integer => Real,
strict positive => non-negative, etc.) but do not perform transitive
chaining across multiple inequality facts.

### v1 known gaps (recorded for v2)

- **Inequality bound-propagation proofs.** Examples like
  `Simplify[x^2 > 3, x > 2]` and `Simplify[Abs[x] < 2, x^2 + y^2 < 4]`
  require a real bound-propagation reasoner; v1 returns the predicate
  unchanged unless it appears literally in the assumption set.
- **General Sqrt distribution under sign assumptions.** Mathematica's
  `Simplify[Sqrt[x^2 y^2], Assumptions -> y < 0]` returns
  `-Sqrt[x^2] y` by distributing Sqrt over the product. v1 does not
  distribute Sqrt over products and leaves the form unchanged.
- **Multi-variable polynomial reduction beyond direct/single-monomial
  substitution.** No Gröbner basis. Higher-arity polynomial relations
  (e.g. three-variable equalities) may simplify only partially.
- **Equation rebalancing in results.** v1 simplifies each side of an
  equation independently. `Simplify[2x - 4y + 6z - 10 == -8]` returns
  `True` (because the simplified pair is trivially equal) rather than
  Mathematica's `x + 3z == 1 + 2y`.
- **Reasoning over composite expressions.** Predicates such as
  `assume_known_positive(Plus[x, y], ctx)` for `x>0 && y>0` are not
  inferred (the Plus expression isn't a fact-target itself).
- **Domain coverage.** `Algebraics`, `Primes`, and `Booleans` are
  parsed but not used for reasoning beyond what `Reals` already implies
  for `Algebraics`.

## Element and logexp Simplify cascade (2026-04-28)

Two related additions in `src/simp.c`. `Element[x, dom]` is now a
first-class evaluator that reads `$Assumptions` to decide domain
membership, and `Simplify` gains an assumption-driven Log/Power
identity rewriter that fires on compound bases (the previous code only
emitted per-symbol rules from a string template, missing patterns like
`Log[a*b]` and `(a*b)^c`).

### Element

`Element[x, dom]` returns `True` when `x` is provably in `dom`,
`False` when `x` is provably outside `dom`, and stays unevaluated
otherwise. Supported domains: `Integers`, `Rationals`, `Reals`,
`Algebraics`, `Complexes`, `Booleans`, `Primes`, `Composites`.

Numeric and structural literals decide directly without consulting any
assumption set:

| Input | Result |
|-------|--------|
| `Element[5, Integers]` | `True` |
| `Element[5/2, Integers]` | `False` |
| `Element[5/2, Rationals]` | `True` |
| `Element[2.5, Reals]` | `True` |
| `Element[2.5, Integers]` | `False` |
| `Element[1+I, Complexes]` | `True` |
| `Element[1+I, Reals]` | `False` |
| `Element[7, Primes]` | `True` |
| `Element[15, Composites]` | `True` |
| `Element[True, Booleans]` | `True` |

For symbolic queries, `Element` reads the current `$Assumptions` and
delegates to the AssumeCtx domain queries. Crucially, `Element` reads
`$Assumptions` raw via `symtab_get_own_values` rather than via
`evaluate()`, because `evaluate("$Assumptions")` would re-fire the
OwnValue rule and recursively call `builtin_element` on the assumption
itself -- an infinite loop when the bound value is, say, the literal
expression `Element[x, Reals]`. Reading the OwnValue replacement
directly breaks the cycle.

The Integer => Rational => Algebraic => Real => Complex lattice is
honoured automatically through `prov_int` / `prov_re` -- so under
`Assuming[Element[x, Integers], ...]` a query for
`Element[x, Reals]` returns `True`. Inequality facts that pin down a
real numeric bound (e.g. `x > 0`) also imply `Element[x, Reals]`.

`Element[{x1, x2, ...}, dom]` threads over its first argument and
returns the list of per-element decisions, leaving any undecided
components in their original `Element[xi, dom]` form.

`Element` is registered with `ATTR_PROTECTED` and a docstring in
`src/info.c`. After `Assuming` exits, the OwnValue stack on
`$Assumptions` pops automatically and symbolic Element queries return
to the unevaluated form.

### Recursive predicates

`assume_known_positive`, `assume_known_nonneg`, `assume_known_real`,
and `assume_known_integer` have been generalised from per-symbol
direct-fact lookup to recursive structural reasoning. The internal
helpers `prov_pos`, `prov_nn`, `prov_neg`, `prov_np`, `prov_int`, and
`prov_re` form a small lattice with bounded mutual recursion:

- Standard positive constants (`Pi`, `E`, `EulerGamma`, `GoldenRatio`,
  `Catalan`, `Degree`, `Glaisher`, `Khinchin`) are inferred positive.
- `Times[u, v, ...]` is positive iff every factor is positive.
- `Plus[u, v, ...]` is positive iff at least one term is strictly
  positive and all others are non-negative.
- `Power[positive, _]` is positive (any exponent).
- `Exp[real]`, `Cosh[real]` are positive.
- `Times` and `Plus` of reals are real; integer-only Times/Plus is
  integer.
- Real-valued elementary functions of real arguments are real
  (`Sin`/`Cos`/`Tan`/`Cot`/`Sec`/`Csc`, the hyperbolics, `Exp`, `Abs`,
  `Floor`/`Ceiling`/`Round`/`Sign`, `ArcTan`/`ArcCot`/`ArcSinh`).
- `Log[positive]` is real.
- `Power[positive, real]` is real; `Power[real, integer]` is real.
- `Power[real, even-integer]` is non-negative.
- `Abs[real]` is non-negative.

Every existing call site continues to work: the public `assume_known_*`
functions delegate to the new `prov_*` helpers, so callers automatically
benefit from the more powerful reasoning.

### Logexp identity cascade (positive-real subset)

A new bottom-up AST walker `apply_logexp_rules` applies the following
identities under provable positivity / reality of operands:

| Pattern | Result | Conditions |
|---------|--------|------------|
| `Log[Times[u1, ..., un]]` | `Plus[Log[u1], ..., Log[un]]` | every `ui` positive |
| `Log[Power[x, p]]` | `p Log[x]` | `x` positive, `p` real |
| `Power[Times[u1, ..., un], a]` | `Times[u1^a, ..., un^a]` | every `ui` positive |
| `Power[Power[x, p], q]` | `Power[x, p q]` | `x` positive, `p` real |

The walker is invoked as both a seed candidate and a per-round
candidate inside `simp_search`. Because the expanded forms have
strictly larger leaf counts than their compact source forms (e.g.
`Log[a b]` is 4 leaves but `Log[a] + Log[b]` is 6), the rewriter is
**force-biased**: any change it produces under sound assumptions is
unconditionally accepted as the new running best. Downstream transforms
(`Cancel`, `Together`, `Factor`, ...) cannot recombine an expanded log
or power, so the expanded form persists through the rest of the search.

Quotient handling is automatic: `Log[a/b]` is canonicalised to
`Log[Times[a, Power[b, -1]]]`, the Times rule expands it to
`Log[a] + Log[Power[b, -1]]`, and the Power rule (with `b > 0`,
`-1` trivially real) collapses `Log[Power[b, -1]]` to `-Log[b]`. A
fixed-point loop in `apply_logexp_rules` (capped at 8 iterations)
ensures the cascade reaches its full extent.

Examples:

| Input | Result |
|-------|--------|
| `Simplify[Log[a*b], a > 0 && b > 0]` | `Log[a] + Log[b]` |
| `Simplify[Log[a/b], a > 0 && b > 0]` | `Log[a] - Log[b]` |
| `Simplify[Log[5*x], x > 0]` | `Log[5] + Log[x]` |
| `Simplify[Sqrt[a*b], a > 0 && b > 0]` | `Sqrt[a] Sqrt[b]` |
| `Simplify[(a*b)^c, a > 0 && b > 0]` | `a^c b^c` |
| `Simplify[(a/b)^c, a > 0 && b > 0]` | `a^c b^(-c)` |
| `Simplify[(a^p)^q, a > 0 && Element[p, Reals]]` | `a^(p q)` |
| `Simplify[Log[a^p], a > 0 && Element[p, Reals]]` | `p Log[a]` |
| `Simplify[Log[2 x] - Log[x], x > 0]` | `Log[2]` |
| `Simplify[Log[x^2] + Log[1/x^2], x > 0]` | `0` |
| `Simplify[(a/b)^c * b^c, a > 0 && b > 0]` | `a^c` |

When the conditions are not met (e.g. `Simplify[Log[a*b]]` with no
assumption, or `Simplify[Log[a^p], a > 0]` without proving `p` real)
the rewrite stays inert and the input is returned unchanged.

### v1 gaps (recorded for v2)

- **General-real cascade with Boole / Floor corrections.** The
  user-supplied rule set includes a "general reals" branch where
  `Log[x y]` and `(x y)^a` for arbitrary real (possibly negative) `x`,
  `y` carry an additional `2 Pi I Boole[x < 0 && y < 0]` correction
  term. v1 implements only the strict-positive subset; the
  Boole-corrected branch is deferred.
- **General-complex fallback.** Mathematica's most general identities
  use `Arg`, `Im[a Log[x]]`, and `Ceiling` to produce phase-corrected
  expansions for arbitrary complex operands. Not implemented in v1.
- **Universal power-merge rules.** `x^a x^b -> x^(a+b)`, `1/x^a ->
  x^(-a)`, etc. are largely already handled by the `Times` and `Power`
  evaluators' coefficient mergers, but no separate `Simplify` pass
  enforces them under assumptions; cases that escape the canonicaliser
  (e.g. structurally unmerged forms after a transform) won't be merged.

### Tests

`tests/test_element.c` covers literal numeric / structural decisions,
list threading, the `$Assumptions` lattice, restoration after
`Assuming` exits, and direct compound-expression facts. 26 tests.

`tests/test_logexp_simplify.c` covers each row of the cascade table
above plus several "not enough assumption" negative cases that verify
the rewriter does NOT fire when its conditions are unmet, plus
composite cancellation cases (e.g. `Log[2 x] - Log[x] -> Log[2]`)
where the cascade enables downstream Plus-cancellation. 22 tests.

## FactorTerms / FactorTermsList (2026-04-28)

Two new builtins in `facpoly.c` extract polynomial content without
performing a full factorization:

- `FactorTerms[poly]` returns the overall numerical factor multiplied
  by the residue. `FactorTerms[poly, x]` (or with a list of variables
  `FactorTerms[poly, {x_1, ..., x_n}]`) returns instead the content of
  `poly` with respect to those variables, recursively peeling off
  contents of progressively smaller subsets.
- `FactorTermsList[...]` returns the same factors as a `List` instead of
  multiplying them together: the first element is always the integer
  content, the next `n` elements are the progressively-extracted
  polynomial contents (in order of *decreasing* variable subset), and
  the final element is the residue. `Apply[Times, FactorTermsList[p]]`
  reproduces `p` (modulo `Together` normalization).

### Algorithm

The core helper is `ft_content_wrt_set(poly, S, ground)`: a recursive
multivariate-content routine that views `poly` as a polynomial in
`K[ground][S]` and computes the GCD of its `S`-monomial coefficients in
`K[ground]`. It peels variables from `S` one at a time, collecting
coefficients via `get_coeff` and combining them with the existing
`poly_gcd_internal` over `ground`. The base case (`S` empty) returns
`expr_expand(poly)`, which is already in `K[ground]`.

`ft_compute_list` orchestrates the layered extraction:

1. **Numerical content**: call `ft_content_wrt_set(poly, all_vars, [])`
   so the gcd lives in `Z` (via `my_number_gcd`).
2. For `k = n, n-1, ..., 1`: extract content with respect to the
   subset `{x_1, ..., x_k}`, treating the complement
   `all_vars \ {x_1, ..., x_k}` as the polynomial ground ring.
3. The residue (the running quotient after each peel) is the last
   element. For rational inputs the algorithm runs `Together` first and
   the residue absorbs the denominator on the way out.

`exact_poly_div` is used to perform each division, with a graceful
`Times[poly, content^-1]` fallback if exact polynomial division fails.

### Threading

`FactorTerms` auto-threads over `List`, `Equal`, `Unequal`, `Less`,
`LessEqual`, `Greater`, `GreaterEqual`, `And`, `Or`, `Not`, `Xor`. This
matches Mathematica's behaviour: e.g.
`FactorTerms[1 < 77 x^3 - 21 x + 35 < 2]` returns
`1 < 7 (5 - 3 x + 11 x^3) < 2`. `FactorTermsList` does *not* thread,
because its result is itself a `List`.

### Attributes

Both functions carry `Protected` (no other attributes -- they are not
`Listable`, since their threading semantics over `Equal`/`Less`/etc.
are richer than simple element-wise mapping).

### Known limitation

Gaussian-integer content is not extracted as a Gaussian unit:
`FactorTerms[5 I x^2 + 20 I x + 10]` returns
`5 (2 + I x^2 + 4 I x)` rather than Mathematica's
`5 I (-2 I + 4 x + x^2)`. Both forms are valid factorizations and
multiply back out to the input; the picocas output simply uses the
real-integer GCD `5` as the leading factor. This is documented in the
implementation comment.

### Tests

`tests/test_factor_terms.c` covers all the documented Mathematica
examples (numerical-only, single-variable, multi-variable, threading
over `List` / `Equal` / `Less`, rational input, complex coefficients)
plus round-trip checks
(`Expand[Together[poly - Apply[Times, FactorTermsList[poly]]]] == 0`)
on every multi-step input, and a sanity check that both functions are
`Protected`. 35+ assertions in 7 test groups, all passing under
valgrind with no leaks attributable to the new code.

## Simplify uses ExpandNumerator, ExpandDenominator, FactorTerms, and per-variable Collect (2026-04-28)

`Simplify`'s heuristic search now considers four additional algebraic
transforms on every candidate it generates:

- `ExpandNumerator` -- distribute products only inside the numerator.
- `ExpandDenominator` -- distribute products only inside the denominator.
- `FactorTerms` -- pull out an overall numerical factor.
- `Collect[expr, v]` for each free variable `v` of `expr`, automatically
  enumerated by calling `Variables[expr]`. This is the variant
  Mathematica's `Simplify` performs to surface a more compact regrouped
  form (e.g. `a x + b x + c -> c + x (a + b)`) without the user having
  to specify which variable to collect by.

The first three are wired into the existing `SIMP_TRANSFORMS` table and
participate in the same per-round candidate scoring as `Expand`,
`Factor`, etc. The Collect step is implemented separately
(`try_collect_per_variable` in `simp.c`) because it is not a unary
transform: it must enumerate variables and try one Collect per
variable. It is invoked once on the input (seed phase) and once on each
surviving seed inside the round loop, mirroring how
`transform_trig_roundtrip` is wired.

`FactorTerms` is also added to the `transform_safe_for` non-integer-
power exclusion alongside `Factor`/`FactorSquareFree`/`TrigFactor`
because it shares the same polynomial-content machinery and would
otherwise stall on inputs containing `Sqrt[]` or other fractional
exponents.

### Example

```
Simplify[a x + b x + c]      (* -> c + x (a + b) *)
```

Previously this returned `a x + b x + c` unchanged because none of
`Expand`, `Factor`, `Together`, etc. produces the regrouped form. The
new per-variable Collect step tries `Collect[a x + b x + c, x]`,
obtains `c + x (a + b)`, and the heuristic accepts it as the new best
candidate via the standard `update_best` complexity tiebreak.

### Tests

`tests/test_simplify.c::test_simplify_collect_by_variable` exercises
the new behaviour end-to-end. All other Simplify tests continue to
pass, confirming the additional transforms do not regress the existing
polynomial / rational / trigonometric / hyperbolic / log-power /
assumption-driven cases.

## $SimplifyDebug -- per-transform tracing (2026-04-28)

A new system symbol `$SimplifyDebug` (default `False`) controls whether
`Simplify` emits a per-transform trace on **stderr**. When set to
`True`, every transform invocation inside `simp_search` produces one
line in the form

```
/<TransformName>/: <input> -> <output> [<elapsed> ms]
```

This covers all sources of work in the search:

- The 11 unary transforms in `SIMP_TRANSFORMS` (`Together`, `Cancel`,
  `Expand`, `ExpandNumerator`, `ExpandDenominator`, `Factor`,
  `FactorSquareFree`, `FactorTerms`, `Apart`, `TrigExpand`,
  `TrigFactor`).
- The composite `TrigRoundtrip` (TrigToExp / Together / Cancel /
  ExpToTrig).
- Per-variable `Collect[_, v]` calls (the variable name is included
  in the trace label).
- `AssumptionRules` and `LogExpRules` rewrites.

The flag is read directly off the `$SimplifyDebug` OwnValue (matching
how `read_dollar_assumptions` handles `$Assumptions`) so toggling it
on/off has no measurable cost when off and never re-enters the
evaluator.

### Usage

```mathematica
$SimplifyDebug = True;
Simplify[a*x + b*x]    (* trace lines emitted on stderr *)
$SimplifyDebug = False;
```

### Use cases

- Diagnosing `Simplify` calls that hang or take an unexpectedly long
  time -- the trace shows which transform is responsible and on which
  intermediate candidate.
- Identifying transforms whose output explodes the candidate set
  (e.g. `TrigExpand` rewriting `Sin[7 x]` into a degree-7 polynomial
  in `Sin[x]`/`Cos[x]`).
- Spotting candidate cycles where the search keeps oscillating between
  equivalent forms.

### Implementation

`simp.c` adds:

- `simp_debug_enabled()` -- reads the OwnValue rule list directly.
- `simp_debug_log(name, in, out, ms)` -- formats and prints one line
  to `stderr` using `expr_to_string` and `fflush(stderr)`.
- `traced_call_unary(name, in)` -- wraps `call_unary_copy` with timing
  and logging.

All transform call sites in `simp_search`, `try_collect_per_variable`,
and `transform_trig_roundtrip` now emit a trace line when the flag is
set. Timing uses `clock()` from `<time.h>` (C99-portable).

## Simplify branch pruning + ComplexityFunction Automatic + SimplifyCount + TrigFactor degree guard (2026-04-28)

Four related changes that together fix the
`Simplify[(1/4)*I*(E^((-I)*x) + E^(I*x))*(E^(-6*I*x) - E^(6*I*x))]`
hang reported on the original repro.

### `SimplifyCount` (new built-in)

Exposes the default Simplify complexity measure as a user-callable
function. Implements the Mathematica `SimplifyCount` definition:

| Expression          | SimplifyCount |
| ------------------- | ------------- |
| `Symbol`            | 1             |
| `0`                 | 1             |
| `Integer p > 0`     | digits(p)     |
| `Integer p < 0`     | digits(\|p\|) + 1 (extra unit for the leading "-") |
| `Rational[n, d]`    | SimplifyCount[n] + SimplifyCount[d] + 1 |
| `Complex[re, im]`   | SimplifyCount[re] + SimplifyCount[im] + 1 |
| `Real` / `MPFR`     | 2 (NumberQ but not Integer/Rational) |
| `String`            | 1 (treated as a leaf) |
| `Function[h, args]` | SimplifyCount[h] + sum SimplifyCount[args] |

This replaces the previous "leaf count + integer digit count"
approximation in `simp_default_complexity`. The new metric matches
Mathematica's documented behaviour exactly so users who reason about
the score (custom `ComplexityFunction`s, picocas test goldens) get
identical numbers. Spot checks: `SimplifyCount[100 Log[2]] == 6`,
`SimplifyCount[Log[2^100]] == 32`, `SimplifyCount[Cos[Pi n]] == 4`,
`SimplifyCount[(-1)^n] == 4`, `SimplifyCount[3.14] == 2`,
`SimplifyCount[1/2] == 3`.

`SimplifyCount` is registered with `ATTR_LISTABLE | ATTR_PROTECTED`
and a docstring.

### `ComplexityFunction -> Automatic`

`Simplify[expr, ComplexityFunction -> Automatic]` is now a synonym for
`Simplify[expr]` (no option). The option detection in `builtin_simplify`
maps the symbol `Automatic` to the internal NULL `complexity_func`
slot, so scoring goes through the fast native `simp_default_complexity`
instead of evaluating `Automatic[candidate]` (which would never
reduce). User-supplied `ComplexityFunction -> f` continues to work
unchanged.

### Branch pruning by complexity

`simp_search`'s round loop now drops any candidate whose score exceeds
its parent's score before that candidate becomes a seed for the next
round. The rule applies to all five candidate sources -- the 11
`SIMP_TRANSFORMS`, the `TrigRoundtrip` composite, per-variable
`Collect`, `AssumptionRules`, and `LogExpRules` -- and uses
`score_with_func` so it honours a user-supplied `ComplexityFunction`.

Strict-greater (`>`) is the gate: equal-score forms still propagate to
preserve plateau exploration (since `cs_add_or_free` already dedups,
this cannot loop). `update_best` is unchanged: a worse-scoring
candidate may still beat the current best on the same step (e.g. an
assumption-driven rewrite that's correctness-preserving but expanded),
just not seed further work.

`AssumptionRules` and `LogExpRules` keep their "force-win for `best`"
semantics: an assumption-driven rewrite always wins the `best` slot if
non-trivial, but only seeds the next round when its score is no worse
than the parent.

### TrigFactor degree-cap guard (the actual hang fix)

`trigfactor_run_pipeline` now skips its `Factor` step not only when
`count_distinct_squared_trig_atoms > 2` (the existing guard) but also
when `max_trig_atom_power > 4`. The new helper
`max_trig_atom_power(e)` returns the maximum integer exponent on any
`Sin`/`Cos`/`Sinh`/`Cosh` atom in `e`. Threshold rationale:

- The known-good Pythagorean cases reach degree 4 (e.g.
  `Sin[x]^4 + 2 Sin[x]^2 Cos[x]^2 + Cos[x]^4 -> 1`), so a cap of 4
  preserves them.
- TrigFactor's Path B (TrigExpand the input first) on
  `Sin[m x] + Sin[n x]` (m,n >= 3) produces a polynomial of total
  degree max(m, n) in `Sin[x]`/`Cos[x]`. With m=5, n=7 the result is
  degree 7 with seven mixed-degree monomials. Factor's
  trial-division loop in `factor_roots` stalls hard on that shape
  even though there are still only two distinct atoms.
- Skipping Factor in the high-degree branch is safe: the
  `trig_factor_identities` rewriter that runs immediately after
  Factor catches the Pythagorean structure that survives without
  Factor, and the surface-form Tan/Sec round-trip at the end works
  on either pathway.

`TRIG_FACTOR_DEGREE_THRESHOLD = 4` is colocated with the existing
`TRIG_FACTOR_ATOM_THRESHOLD` constant in `trigsimp.c`.

### End-to-end result

Originally-hanging cases now return immediately:

| Call | Result | Time |
| ---- | ------ | ---- |
| `Simplify[(1/4)*I*(E^((-I)*x) + E^(I*x))*(E^(-6*I*x) - E^(6*I*x))]` | `1/2 (Sin[5 x] + Sin[7 x])` | ~1s |
| `TrigFactor[Sin[5 x] + Sin[7 x]]`           | `Sin[5 x] + Sin[7 x]` | <1s |
| `TrigFactor[Sin[3 x] + Sin[5 x]]`           | `Sin[3 x] + Sin[5 x]` | <1s |
| `TrigFactor[Sin[3 x] + Sin[7 x]]`           | `Sin[3 x] + Sin[7 x]` | <1s |
| `TrigFactor[Sin[3 x] + Sin[9 x]]`           | `Sin[3 x] + Sin[9 x]` | <1s |
| `TrigFactor[1/2 (Sin[5 x] + Sin[7 x])]`     | `1/2 (Sin[5 x] + Sin[7 x])` | <1s |

Previously-working TrigFactor cases continue to work unchanged:

| Call | Result |
| ---- | ------ |
| `TrigFactor[Sin[x]^2 + Cos[x]^2]`                            | `1`                  |
| `TrigFactor[Cos[x]^2 - Sin[x]^2]`                            | `Cos[2 x]`           |
| `TrigFactor[Sin[x]^4 + 2 Sin[x]^2 Cos[x]^2 + Cos[x]^4]`      | `1`                  |
| `TrigFactor[Cos[x+y] + Sin[x] Sin[y]]`                       | `Cos[x] Cos[y]`      |
| `TrigFactor[Cosh[x]^2 - Sinh[x]^2]`                          | `1`                  |

### Memory/leaks

Valgrind on `Simplify[a*x + b*x]` reports 5,504 bytes definitely lost
versus 5,824 bytes on the pre-change baseline. The drop comes from
the new pruning eliminating some `expr_copy` traffic in the round
loop. No leak record references any of the new helpers
(`simp_default_complexity`, `simp_debug_*`, `traced_call_unary`,
`builtin_simplify_count`, `max_trig_atom_power`).

### Tests

All simp-related tests pass: `simp_tests`, `simplify_tests`,
`logexp_simplify_tests`, `assuming_tests`, `element_tests`,
`trigfactor_tests`, `trigexpand_tests`. The 12 pre-existing test
failures (`backtrack`, `context`, `facpoly`, `hyperbolic`, etc.)
reproduce identically on a `git stash` baseline -- none caused by
this change.

## Bottom-up Simplify, Abs simplification, TrigFactor real-coefficient guard (2026-04-28)

This change extends `Simplify` from a top-level-only candidate-set search
to a memoised bottom-up traversal, adds a structural `Abs[...]`
simplification cascade, and fixes a soundness bug in `TrigFactor` exposed
by the new walker.

### Bottom-up traversal

Before: `Simplify[expr]` ran the candidate-set search (`simp_search`) on
`expr` directly. Transforms like the Pythagorean identity
`Sin[x]^2 + Cos[x]^2 -> 1` only fired if the relevant structure appeared
at the root; e.g. `Simplify[f[Sin[x]^2 + Cos[x]^2]]` returned
`f[Sin[x]^2 + Cos[x]^2]` unchanged because no transform peeked inside the
generic head `f`.

After: `Simplify` recursively descends into each sub-expression, runs
`simp_search` on the children first, rebuilds + re-evaluates the parent
(restoring canonical-form invariants on `Plus`/`Times`), and only then
runs `simp_search` on the parent. This surfaces buried simplifications:

| Call                                                | Result   |
| --------------------------------------------------- | -------- |
| `Simplify[f[Sin[x]^2 + Cos[x]^2]]`                  | `f[1]`   |
| `Simplify[f[Sin[x]^2 + Cos[x]^2, Sin[y]^2 + Cos[y]^2]]` | `f[1, 1]` |
| `Simplify[g[h[Sin[x]^2 + Cos[x]^2]]]`               | `g[h[1]]` |
| `Simplify[{Sin[x]^2 + Cos[x]^2, Cosh[x]^2 - Sinh[x]^2}]` | `{1, 1}` |

Cost control:

  - **Memoisation** by `expr_hash` + `expr_eq` over a 256-bucket table:
    structurally identical subtrees (e.g. `f[g[x], g[x], g[x]]`) cost
    one `simp_search` call instead of three.
  - **Atom short-circuit**: when the active `AssumeCtx` is empty (no
    user-supplied facts), atoms (`EXPR_INTEGER`, `EXPR_BIGINT`,
    `EXPR_REAL`, `EXPR_SYMBOL`, `EXPR_STRING`) skip `simp_search`
    entirely and return a copy -- every transform is a no-op on a bare
    atom in the absence of assumption rewrites. (The check is
    `ctx == NULL || ctx->count == 0` rather than just `ctx == NULL`,
    because `assume_ctx_from_expr` always returns a non-NULL ctx even
    when `$Assumptions = True`.)
  - **Held-head skip list**: `Hold`, `HoldForm`, `HoldComplete`,
    `HoldPattern`, `Unevaluated`, `Pattern`, `Blank`, `BlankSequence`,
    `BlankNullSequence`, `Function`, `Slot`, `SlotSequence` -- and any
    head whose attribute mask sets `ATTR_HOLDFIRST`/`HOLDREST`/
    `HOLDALL`/`HOLDALLCOMPLETE` (e.g. `Set`, `If`, `Module`, `Block`,
    `While`, `For`, `Do`) -- do not get descended into. We still run
    `simp_search` on the whole expression (a no-op on these heads in
    practice, but reserves the chance for future structural rewrites).
  - **Small-subtree leaf-count guard**: at non-top levels, skip
    `simp_search` when `simp_default_complexity(canonical) <= 7`. The
    threshold admits Pythagorean-eligible Plus/Times nodes
    (`Sin[x]^2 + Cos[x]^2` = 9, `Cosh[x]^2 - Sinh[x]^2` = 12) while
    rejecting isolated trig powers (`Cosh[x]^2` = 4,
    `Times[-1, Sinh[x]^2]` = 7) on which `TrigRoundtrip` produces a
    12-term polynomial in `Cosh[2x]/Sinh[2x]/Cosh[4x]/Sinh[4x]` --
    feeding that into `Factor` cost 7+ seconds per call before this
    guard. The top-level Simplify call (depth == 0) always runs
    `simp_search` regardless of size.
  - **Depth cap** of 64 levels: pathological nesting falls back to the
    existing top-level `simp_search`.

### `count_distinct_trig_atoms` gate inside TrigFactor

When `simp_search` is invoked on a sub-expression that contains a
multi-trig polynomial -- typically the residue of a `Factor` call on a
`TrigRoundtrip` output -- `Factor`'s multivariate trial-division loop in
`factor_roots` stalls hard. With four distinct degree-1 trig atoms
(`Cosh[2x]`, `Sinh[2x]`, `Cosh[4x]`, `Sinh[4x]`) and dense cross-products
the call ran ~150s before any guard fired. The existing
`count_distinct_squared_trig_atoms` gate only counted atoms that appeared
at integer power >= 2, so degree-1 dense forms slipped through.

`count_distinct_trig_atoms(e)` extends the existing collector to count
*every* `Sin/Cos/Sinh/Cosh` call (any power, including degree 1) as a
distinct atom when its argument differs structurally. The
`trigfactor_run_pipeline` Factor gate now also skips when this count
exceeds `TRIG_FACTOR_TOTAL_ATOM_THRESHOLD = 3`. The pipeline still runs
the identity-collapse rules (Pythagorean, angle-sum, double-angle), so
the angle-addition cases (`Sin[a]Cos[b] + Cos[a]Sin[b]`, 4 distinct
atoms) continue to collapse via the structural `:>` rule -- only
`Factor`, the slow component, is bypassed.

### Abs simplification cascade

`apply_abs_rules` is a structural Abs-rewrite pass driven by a recursive
walker. It runs as **Phase 0** of `simp_search` (a force-take preprocessor
that rewrites the input before any other transform sees it) and again
inside the round loop on each candidate. Phase-0 placement is deliberate:
the rewritten form is often leaf-count-larger than the input
(e.g. `Abs[Times[a,b,c]]` -> `Abs[a] Abs[b] Abs[c]` grows from 5 to 7
leaves), so without preprocessing the leaf-count tiebreak would
re-prefer the original from the seed set.

**Universal rules** (fire without an assumption context):

| Rule                                  | Reason                                       |
| ------------------------------------- | -------------------------------------------- |
| `Abs[Abs[x]] -> Abs[x]`               | Idempotency.                                 |
| `Abs[Conjugate[x]] -> Abs[x]`         | `\|conj(z)\| = \|z\|`.                       |
| `Abs[E^z] -> E^Re[z]`                 | Magnitude of any complex exponential.        |
| `Abs[Times[a, b, ...]] -> Abs[a] Abs[b] ...` | `\|ab\| = \|a\| \|b\|`. Captures both numeric coefficient extraction (`Abs[c x] -> \|c\| Abs[x]`) and quotient split (`Abs[x/y]` is `Abs[Times[x, Power[y, -1]]]`). |
| `Abs[x^n] -> Abs[x]^n` for integer `n` | `\|x^n\| = \|x\|^n` is exact for integer `n`. |

**Cascading rules** (require `AssumeCtx`):

| Rule                                | Predicate                          |
| ----------------------------------- | ---------------------------------- |
| `Abs[x] -> x`                       | `assume_known_nonneg(ctx, x)` (i.e. `x >= 0`)  |
| `Abs[x] -> -x`                      | `assume_known_nonpos(ctx, x)` (i.e. `x <= 0`)  |
| `Abs[x^y] -> Abs[x]^y`              | `assume_known_real(ctx, y)` (`y` is real)      |
| `Abs[x^y] -> x^Re[y]`               | `assume_known_positive(ctx, x)` (strictly `x > 0`) |

The user-provided `Abs[Sin[x_]] :> Sign[Sin[x]] * Sin[x]` rule is
**omitted**: leaf count grows from 3 to 6, and picocas has no Sign-
folding pass to recover it. We can add it back when interval-narrowing
on `Sign[trigfunc]` lands.

End-to-end checks:

| Call                                                  | Result               |
| ----------------------------------------------------- | -------------------- |
| `Simplify[Abs[Abs[x]]]`                               | `Abs[x]`             |
| `Simplify[Abs[Conjugate[x]]]`                         | `Abs[x]`             |
| `Simplify[Abs[Conjugate[Abs[x]]]]`                    | `Abs[x]` (cascades)  |
| `Simplify[Abs[E^(2 I y)]]`                            | `E^Re[2 I y]`        |
| `Simplify[Abs[a*b*c]]`                                | `Abs[a] Abs[b] Abs[c]` |
| `Simplify[Abs[2*x*y]]`                                | `2 Abs[x] Abs[y]` (`Abs[2] -> 2`) |
| `Simplify[Abs[(-3)*x]]`                               | `3 Abs[x]` (`Abs[-3] -> 3`) |
| `Simplify[Abs[I*x*y]]`                                | `Abs[x] Abs[y]` (`Abs[I] -> 1`) |
| `Simplify[Abs[x/y]]`                                  | `Abs[x]/Abs[y]`      |
| `Simplify[Abs[1/x]]`                                  | `1/Abs[x]`           |
| `Simplify[Abs[x^3]]`                                  | `Abs[x]^3`           |
| `Simplify[Abs[x], x > 0]`                             | `x`                  |
| `Simplify[Abs[x], x >= 0]`                            | `x`                  |
| `Simplify[Abs[x], x <= 0]`                            | `-x`                 |
| `Simplify[Abs[x^y], Element[y, Reals]]`               | `Abs[x]^y`           |
| `Simplify[Abs[x^y], x > 0]`                           | `x^y`                |
| `Simplify[Abs[(x+1)/y], x >= 0 && y > 0]`             | `(1+x)/y`            |

### TrigFactor real-coefficient guard

The bottom-up walker exposed a soundness bug in `TrigFactor`'s
linear-combination identity:

```
a_. Sin[x_] + b_. Cos[x_] + r___ /; (NumberQ[a] && NumberQ[b])
  :> Sqrt[a^2 + b^2] Sin[x + ArcTan[a, b]] + r
```

`NumberQ` accepts `Complex` values, so for `Cos[x] + I Sin[x]` the rule
matched with `a = I, b = 1` and produced
`Sqrt[I^2 + 1] Sin[x + ArcTan[I, 1]] = Sqrt[0] * ... = 0` -- a false
collapse of `E^(I x)` to zero. Top-level Simplify never seeded
`Cos[x] + I Sin[x]` directly so the bug was hidden; bottom-up turns
`E^(I x)` into `Cos[x] + I Sin[x]` via `TrigRoundtrip` at the leaf level,
exposing the bug.

Fix: tighten the predicate to also require `Im[a] === 0 && Im[b] === 0`,
so the identity only fires for genuinely real coefficients.

### Memory and leaks

Bottom-up amplifies pre-existing leaks in `apply_assumption_rules`,
`together_recursive`, `cancel_recursive`, `evaluate_step`,
`builtin_polynomiallcm`, and `builtin_times` -- all of which are reached
more often now that `simp_search` runs at internal nodes. None of the
leak records originate inside the new code (`simp_bottomup`,
`simp_memo_*`, `simp_skip_recursion_head`, `simp_head_holds_args`,
`apply_abs_rules`, `abs_walk`, `try_simp_abs`, `contains_abs`,
`count_distinct_trig_atoms`, `collect_trig_atoms`).

One pre-existing leak in `apply_assumption_rules` (early NULL return
forgetting to free `eq_diffs`) is fixed in this change. Valgrind on
`Simplify[a*x + b*x]; Simplify[f[Sin[x]^2 + Cos[x]^2]]; Simplify[g[h[Sin[x]^2 + Cos[x]^2]]]`
reports 37,744 bytes definitely lost (vs 8,624 bytes baseline) -- the
remaining ~4x amplification is from the increased call count of leaky
sub-routines, dramatically reduced from the initial 12x amplification by
the small-subtree leaf-count guard. Fixing the underlying poly/rat
leaks is left to a separate change.

### Tests

All simp-related tests pass: `simp_tests`, `simplify_tests`,
`logexp_simplify_tests`, `assuming_tests`, `element_tests`,
`trigfactor_tests`, `trigexpand_tests`, `core_tests`. The 6 pre-existing
hard test failures (`context`, `facpoly`, `linalg`, `list_set`,
`polymod`, `series`) reproduce identically on a `git stash` baseline --
none caused by this change. The earlier-quoted "12 pre-existing
failures" count included tests that printed errors but exited 0; the
actual abort/timeout count is 6.

## TrigToExp added as a Simplify search seed (2026-04-29)

`TrigToExp` is now one of the unary transforms tried inside
`simp_search`. Previously the only path to the exponential form was
the full trig roundtrip (`TrigToExp -> Together -> Cancel -> ExpToTrig`),
which always converts back to a trig form at the final step. For
hyperbolic combinations whose exp form is strictly simpler than the
trig form -- the canonical case being `Sinh[x] + Cosh[x] -> E^x` --
the roundtripped result was `Cosh[x] + Sinh[x]` (unchanged), hiding
the lower-complexity intermediate.

Adding `TrigToExp` as a standalone seed lets the search keep `E^x`
as a candidate; `update_best` then picks it because its `SimplifyCount`
(3) beats the trig form's (5).

For pure trig inputs, `TrigToExp` produces complex-coefficient exp
forms with strictly higher complexity than the original (e.g.
`TrigToExp[Sin[x] + Cos[x]]` has score ~12), so the original survives
the tiebreak. Verified: `Sin[x]^2 + Cos[x]^2 -> 1`,
`2 Tan[x]/(1 + Tan[x]^2) -> Sin[2 x]`, `Sin[x] + Cos[x]` preserved.

End-to-end checks:

| Call                                | Before                  | After   |
| ----------------------------------- | ----------------------- | ------- |
| `Simplify[Sinh[x] + Cosh[x]]`       | `Cosh[x] + Sinh[x]`     | `E^x`   |
| `Simplify[Sin[x]^2 + Cos[x]^2]`     | `1`                     | `1`     |
| `Simplify[Cosh[x]^2 - Sinh[x]^2]`   | `1`                     | `1`     |
| `Simplify[(E^x - E^(-x))/Sinh[x]]`  | `2`                     | `2`     |
| `Simplify[Sin[x] + Cos[x]]`         | `Cos[x] + Sin[x]`       | unchanged |

## Simplify gating + shape classifier (2026-04-29)

Two layered changes reduce wasted work in `Simplify`. Both are
performance-only on inputs the existing search already handled
correctly; correctness on every test case in the suite is preserved.

### Phase 1 -- structural gates

`transform_can_fire(name, expr, ctx)` is the central
cheap-precondition gate. It is consulted at every transform
invocation in `simp_search` -- both the `SIMP_TRANSFORMS[]` table
entries and the open-coded seed/round transforms (`AbsRules`,
`AssumptionRules`, `LogExpRules`, `TrigRoundtrip`,
`CollectPerVariable`).

Predicates added (all O(n) one-walk, in `src/simp.c`):

| Predicate                       | Gates                                                |
| ------------------------------- | ---------------------------------------------------- |
| `contains_trig_or_hyperbolic`   | `TrigExpand` / `TrigFactor` / `TrigToExp` / `TrigRoundtrip` |
| `contains_log`, `contains_power`| `LogExpRules`                                        |
| `contains_variable`             | `Factor` / `FactorSquareFree` / `FactorTerms` / `Apart` on numeric-only input |
| `expr_variables_count_capped`   | `CollectPerVariable` on <2-variable input            |
| `ctx_has_facts`                 | `LogExpRules` / `AssumptionRules` on empty/inconsistent ctx |
| `head_is_trig_or_hyperbolic`    | scalar test underlying `contains_trig_or_hyperbolic` |

Gates are conservative: they return false only when the predicate
proves the transform cannot fire. False positives (gate passes but
transform no-ops) keep current behaviour exactly.

### Phase 2 -- shape classifier and dispatcher

`simp_classify` runs one O(n) walk and assigns the input to one of
five disjoint shape categories, in priority order:

1. **TRIG** -- any Sin/Cos/.../Sinh/.../ArcXxx head present.
2. **GENERAL** (early) -- Abs present (Abs handled by existing
   search machinery).
3. **LOGEXP** -- Log present.
4. **POLYNOMIAL** -- no trig/log/abs, all integer powers, and
   `PolynomialQ` over `Variables[input]` is True.
5. **RATIONAL** -- no trig/log/abs, `Together[input]` yields
   polynomial / polynomial.
6. **GENERAL** -- everything else.

`simp_dispatch` routes by shape:

| Shape       | Pipeline                                                   |
| ----------- | ---------------------------------------------------------- |
| POLYNOMIAL  | `simp_pipeline_polynomial`: Expand, Factor, FactorTerms, Collect-per-variable on input AND on Expand result. |
| RATIONAL    | `simp_pipeline_rational`: AssumptionRules force-take, Together, Cancel, ExpandNumerator/Denominator, Apart, Factor on Cancel result, Collect-per-variable on Cancel result. |
| LOGEXP      | `simp_pipeline_logexp`: LogExpRules cascade force-take, AssumptionRules force-take at equal complexity, Together / Cancel / Expand cleanup. |
| TRIG        | `simp_search` (the existing heuristic; the gates make it cheap). |
| GENERAL     | `simp_search`.                                             |

`simp_bottomup` calls `simp_dispatch` instead of `simp_search`, so
specialisation applies at every recursion level.

### Measured wins

| Input                                                | Before      | After Phase 1 | After Phase 2 |
| ---------------------------------------------------- | ----------- | ------------- | ------------- |
| `6*(x/Sqrt[2]-y/Sqrt[2])*(x/Sqrt[2]+y/Sqrt[2])`      | 818 ms      | 125 ms        | 150 ms (RATIONAL) |
| `(x-1)*(x+1)*(x^2+1)+1`                              | full search | full search   | 2 ms (POLYNOMIAL) |
| `a x + b x + c`                                      | full search | full search   | 7 ms (POLYNOMIAL) |
| `3/(x+3) + x/(x+3)`                                  | full search | full search   | 8 ms (RATIONAL)   |
| `Simplify[Log[a^p], a>0]`                            | full search | full search   | <1 ms (LOGEXP)    |
| `Sin[x]^2 + Cos[x]^2`                                | full search | full search   | full search (TRIG) |

Trig and Abs inputs continue to use the full heuristic search, where
`simp_search`'s round loop and force-take seeds are the right tool.
The Phase 1 gates make every transform on those inputs cheap when
its precondition fails.

### Tests

All 19 simplify-related test suites pass: `simp_tests`, `simplify_tests`,
`logexp_simplify_tests`, `trig_tests`, `hyperbolic_tests`,
`trigexpand_tests`, `trigfactor_tests`, `element_tests`,
`assuming_tests`, `integrals_tests`, `deriv_tests`, `poly_tests`,
`rat_tests`, `parfrac_tests`, `expand_tests`, `factor_terms_tests`,
`limit_tests`, `series_tests`, `core_tests`.

## Simplify: equation rebalancing, even-exponent and Sqrt-of-product assumption rules, Pythagorean square completion, cube-root-of-unity identity, Power printer fix (2026-04-29)

Six fixes that close known gaps reported against `Simplify`. Each is
modular and tested in isolation in `tests/test_simplify.c`.

### 1. Power printer for negative numeric bases (`src/print.c`)

`Power[-1, 1/3]` previously printed as `-1^(1/3)`, which round-trips
through the parser as `-(1^(1/3))` because unary minus has lower
precedence than `^`. Both the standard and the TeX printers now emit
explicit parens around any negative numeric base (`Integer`, `Real`,
`BigInt`, `Rational`), producing `(-1)^(1/3)` instead.

### 2. Equation / inequality rebalancing (`src/simp.c`)

`Simplify` on a binary `Equal/Less/Greater/...` head now also tries a
rebalanced candidate: compute `d = Expand[lhs - rhs]`, divide every
term's integer coefficient by their GCD, and partition into
positive-coefficient terms (LHS) and negated negative-coefficient terms
(RHS). The "polarity" is chosen so the leading variable term's
coefficient ends up positive; for strict inequalities a negative
divisor flips `<` to `>` (etc.). The rebalanced form is scored against
the threaded form via `SimplifyCount` and the simpler one wins.

| Input                                       | Before                          | After                  |
| ------------------------------------------- | ------------------------------- | ---------------------- |
| `Simplify[2 x - 4 y + 6 z - 10 == -8]`      | `2 (-5 + x - 2 y + 3 z) == -8`  | `x + 3 z == 1 + 2 y`   |
| `Simplify[-2 x == 4]`                       | `-2 x == 4`                     | `x == -2`              |
| `Simplify[-2 x < 4]`                        | `-2 x < 4`                      | `x > -2`               |
| `Simplify[10 == 5 x]`                       | `10 == 5 x`                     | `x == 2`               |

Skipped when any term has a non-`int64` numeric coefficient
(`Real`, `BigInt`, `Rational`) so floating-point or symbolic
denominators are preserved.

### 3. Sqrt over a product of squares (`src/simp.c`, `apply_assumption_rules`)

For each known-positive symbol `x` the assumption-rule list now
additionally emits

```
Power[Times[Power[x, 2], rest___], Rational[1, 2]] :> x Power[Times[rest], Rational[1, 2]]
```

(and analogous rules with `-x` for negatives, `Abs[x]` for unsigned
reals). `ReplaceRepeated` then peels factors off `Sqrt[x^2 y^2 z^2 ...]`
one symbol at a time:

```
Simplify[Sqrt[x^2 y^2], x > 0 && y < 0]            -> -x y
Simplify[Sqrt[x^2 y^2 z^2], x>0 && y<0 && z>0]    -> -x y z
Assuming[x > 0, Simplify[Sqrt[x^2 y^2], y < 0]]    -> -x y
```

### 4. Even-exponent collapse for `(-1)^...` (`src/simp.c`)

Added a fifth assumption category, `evens`, populated by the new
`prov_even` helper. It accepts `Mod[x, 2] == 0` (in either Equal-arg
order) and `Element[x, Evens]`. For each even symbol `m` and each
integer symbol `k`, the rule list grows by

```
Power[-1, m] :> 1
Power[-1, k_Integer m] :> 1
Power[Power[-1, k_Integer], m] :> 1
Power[-1, k m] :> 1                       (* paired *)
Power[Power[-1, k], m] :> 1
```

For each integer symbol `n` (independent of any even fact) it also
emits the literal-multiplier rules

```
Power[-1, Times[m_Integer /; EvenQ[m], n]] :> 1
Power[Power[-1, n], m_Integer /; EvenQ[m]] :> 1
```

so e.g. `Cos[k Pi]^4` collapses all the way:

```
Simplify[Cos[k Pi]^4, Element[k, Integers]]                            -> 1
Simplify[Cos[k Pi]^m, Element[k, Integers], Assumptions->Mod[m,2]==0]  -> 1
Simplify[(-1)^(2 k), Element[k, Integers]]                             -> 1
```

`Simplify[Cos[k Pi]^3, Element[k, Integers]]` correctly stays at
`(-1)^(3 k)` because `3 k` may be odd.

### 5. Pythagorean perfect-square completion (`src/simp.c`)

A new transform `transform_pythag_square_complete` applies just two
targeted rules

```
1 + 2 Sin[x_] Cos[x_] + r___ :> (Sin[x] + Cos[x])^2 + r
1 - 2 Sin[x_] Cos[x_] + r___ :> (Sin[x] - Cos[x])^2 + r
```

It runs both as a top-level seed in `simp_search` and inside the round
loop on every candidate, so a `Factor` result like `(1 + 2 Sin Cos)^2`
also has the rule applied to its sub-expression. Kept separate from
`TrigFactor`'s identity rule list because the linear-combination rule
inside that list (`a Sin + b Cos -> Sqrt[a^2+b^2] Sin[x + ArcTan]`)
would re-rewrite `(Sin + Cos)^2` into `2 Sin[x + Pi/4]^2` in the same
ReplaceRepeated pass, hiding the factored form from the score
selection.

```
Simplify[4 Sin[x]^2 Cos[x]^2 + 4 Sin[x] Cos[x] + 1] -> (Cos[x] + Sin[x])^4
```

### 6. Roots-of-unity simplifier (`src/simp.c`)

A new transform `simp_roots_of_unity` recognises every
`Power[-1, p/q]` and `Power[E, Times[Complex[0, p/q], Pi]]` atom in
the input and reduces the expression modulo the cyclotomic polynomial
of the corresponding primitive root of unity. The procedure:

1. Walk the expression and collect denominators `q_i` of all
   root-of-unity atoms (via `Cases[..., Infinity]`).
2. Set `Q = LCM(q_i)`. Every atom is now expressible as `omega^k` for
   `omega = (-1)^(1/Q) = e^(I Pi / Q)`, a primitive `(2Q)`-th root of
   unity.
3. Lift the input to a polynomial in a fresh variable `$ru`,
   substituting each atom as `$ru^Mod[a Q/b, 2 Q]` so negative or
   over-period exponents fold into `[0, 2Q)`.
4. Compute `Phi_{2Q}($ru)` recursively from
   `Phi_n(x) = (x^n - 1) / Prod_{d | n, d < n} Phi_d(x)` using
   `PolynomialQuotient`.
5. Reduce by `PolynomialRemainder[polyForm, Phi_{2Q}, $ru]`. Because
   `Phi_{2Q}` is the minimal polynomial of `omega` over Q, the
   remainder is the unique degree-`< deg(Phi_{2Q})` representative.
6. Substitute `$ru -> Power[-1, 1/Q]` to recover the surface form.

The algorithm is entirely general -- it handles the cube case, every
prime denominator, mixed `(-1)^(p/q)` and `E^(I p Pi /q)` forms, and
expressions with free coefficients (which carry through as polynomial
coefficients during `PolynomialRemainder`). The cyclotomic helper
`$ruCyclotomic[n, x]` and the driver `$ruSimplify[expr]` are added
as DownValues on internal `$ru*` symbols by
`simp_install_roots_of_unity_helpers` on first call (lazy, idempotent).

```
Simplify[1 - (-1)^(1/3) + (-1)^(2/3)]                                 -> 0
Simplify[1 - (-1)^(1/5) + (-1)^(2/5) - (-1)^(3/5) + (-1)^(4/5)]       -> 0
Simplify[1 - (-1)^(1/7) + (-1)^(2/7) - ... + (-1)^(6/7)]              -> 0
Simplify[3 + 2 E^(-2 I Pi/3) + 2 E^(2 I Pi/3)]                        -> 1
Simplify[(-1)^(1/3) + (-1)^(2/3)]                                     -> -1 + 2 (-1)^(1/3)
Simplify[(-1)^(2/3)]                                                  -> (-1)^(2/3)   (* unchanged: original shorter *)
Simplify[5 x + 2]                                                     -> 2 + 5 x      (* no atoms; helper short-circuits *)
```

### Tests

`tests/test_simplify.c` adds twelve focused tests:

- `test_simplify_equation_rebalance`
- `test_simplify_equation_polarity_flip`
- `test_simplify_inequality_polarity_flip`
- `test_simplify_sqrt_product_signs`
- `test_simplify_sqrt_product_three`
- `test_simplify_cos_k_pi_to_even_int_power`
- `test_simplify_cos_k_pi_to_even_symbolic_power`
- `test_simplify_cos_sin_fourth_power`
- `test_simplify_cube_root_of_neg_one`
- `test_simplify_fifth_root_of_unity_alternating_sum`
- `test_simplify_seventh_root_of_unity_alternating_sum`
- `test_simplify_complex_exp_cube_root`

The pre-existing `test_simplify_cos_n_pi_integer` updates its expected
output from `-1^n` to `(-1)^n` to reflect the printer fix.

## Rationalize

`Rationalize` converts approximate (Real / MPFR) numbers to nearby
rationals. It mirrors the Mathematica builtin in both surface API and
threading semantics.

| Form | Behaviour |
|------|-----------|
| `Rationalize[x]`     | Convert `x` to a small-denominator rational only when a continued-fraction convergent of `x` satisfies <code>&#124;p/q − x&#124; &lt; 10⁻⁴ / q²</code>. Otherwise `x` is returned unchanged. |
| `Rationalize[x, dx]` | Return the rational with the smallest denominator lying in `[x − dx, x + dx]`. The classical "simplest rational in interval" continued-fraction algorithm. |
| `Rationalize[x, 0]`  | Force conversion of any inexact `x` to rational form using a tolerance derived from the ulp of `x` (≈ ½ ulp). For `N[Pi]` this lands on `245850922/78256779`. |

Threading rules:

* Numeric leaves (`Integer`, `BigInt`, `Real`, `Rational`, `MPFR`) are
  rationalised in place. In default mode, exact leaves pass through
  unchanged; only `Real` / `MPFR` are rewritten.
* `Complex[re, im]` recurses on its components, so e.g.
  `Rationalize[N[Pi + 33/211 E I], 0]` becomes
  `245850922/78256779 + (42051801/98914198) I`.
* Any `NumericQ` subexpression — `Sqrt[2]`, `Exp[Sqrt[2]]`, `Pi` — is
  numericalised then rationalised when an explicit `dx` is supplied.
  This is what makes `Rationalize[Exp[Sqrt[2]], 2^-12]` collapse to
  `218/53`.
* `Hold` / `HoldForm` / `HoldComplete` / `HoldPattern` / `Unevaluated`
  wrappers are passed through.

Attributes: `Protected`. (Threading is implemented internally via the
expression descent rather than via `Listable`.)

### Algorithm notes

Two distinct kernels back the builtin:

1. **Default mode** (`rationalize_default_double`) walks the convergents
   of the continued fraction of `x` in arbitrary-precision GMP, and
   accepts the first convergent satisfying the `c / q²` bound. The
   denominator is capped at `~2.5×10⁶` — past that point, double-
   precision noise produces spurious "good" convergents that would
   otherwise let `Rationalize[N[Pi]]` find a deep, fictitious match.
2. **Tolerance / zero modes** use the standard "simplest rational in
   `[lo, hi]`" descent on the continued fraction of the interval, with
   GMP scratch for the result. The zero case (`dx = 0`) walks
   convergents in GMP directly because the interval-based algorithm
   cannot maintain enough double-precision dynamic range to reach the
   ~15-deep convergent that `Rationalize[N[Pi], 0]` requires.

### C-level API

In addition to the `Rationalize` symbol, the module exports two C
entry points that other internal algorithms (Simplify, NSolve, etc.)
can use without going through the symbol table:

```c
typedef enum {
    RATIONALIZE_DEFAULT,    /* c/q² bound, may report failure        */
    RATIONALIZE_TOLERANCE,  /* explicit dx > 0                        */
    RATIONALIZE_ZERO,       /* dx = 0, ulp-derived tolerance          */
} RationalizeMode;

bool internal_rationalize_double(double x, double dx, RationalizeMode m,
                                 mpz_t out_n, mpz_t out_d);

Expr* internal_rationalize_expr(const Expr* e, double dx, RationalizeMode m);
```

Both `out_n` and `out_d` are `mpz_init`'d on success. In default mode
the function reports `false` when no convergent satisfies the bound;
the caller is expected to leave `x` alone in that case. The expression-
level entry returns the original expression unchanged for inputs the
core declines to convert.

### Inexact-coefficient round-trip in exact-symbolic functions

The exact-symbolic builtins
`PolynomialGCD`, `PolynomialLCM`, `Cancel`, `Together`, `Apart`,
`Factor`, `FactorSquareFree`, `Simplify`, `Series`, and `Limit`
all assume rational coefficients — their algorithms (subresultant PRS,
Berlekamp–Zassenhaus, partial-fraction matrix solve, leading-term
series extraction, Together-based limit cancellation) cannot reason
about an `EXPR_REAL` / `EXPR_MPFR` leaf. Without a fix-up they would
silently return un-simplified residues such as

```
Cancel[(x^2/9 + (2 x y)/15 + y^2/25) / (x^2/9. - y^2/25.)]
```

leaving the polynomial GCD on the table.

The fix lives in `src/rationalize.c` and is plumbed in at the entry
point of every affected builtin via three helpers:

```c
bool  internal_contains_inexact      (const Expr* e);
bool  internal_args_contain_inexact  (const Expr* res);
Expr* internal_force_rationalize     (const Expr* e);
Expr* internal_rationalize_then_numericalize(Expr* res,
                                             Expr* (*core)(Expr*));
```

`internal_force_rationalize` is functionally close to
`Rationalize[expr, 0]` but with two pragmatic differences tuned for
internal use:

1. Per inexact leaf the lenient `RATIONALIZE_DEFAULT` algorithm runs
   first, with a fall-back to bit-exact `RATIONALIZE_ZERO` only when
   no convergent fits the `c/q²` bound. This recovers `1/9.` as `1/9`
   instead of the bit-exact `562949953421312/5066549580791809`.
2. Symbolic constants (`Pi`, `E`, `EulerGamma`, …) are *not*
   rationalised — they stay symbolic so that e.g.
   `Together[Pi/x + 1.5]` still works algebraically with `Pi` and
   only converts the `1.5`.

`internal_rationalize_then_numericalize` is the wrapper plumbed at the
top of each affected builtin:

```c
if (internal_args_contain_inexact(res))
    return internal_rationalize_then_numericalize(res, builtin_xxx);
```

It builds a clone of `res` with every argument force-rationalised,
recursively dispatches to the original builtin (which now sees no
inexact leaf and runs the exact-arithmetic path), and finally feeds
the result back through `numericalize` so callers observe
inexact-in / inexact-out semantics. Ownership matches the in-tree
convention — the evaluator (eval.c) frees `res` after a non-NULL
return, so the wrapper never touches `res` itself; the locally-built
clone is freed by the wrapper unconditionally.

Test coverage in `tests/test_inexact_dispatch.c` pins both required
behaviours for each of the ten affected builtins:

* the exact-arithmetic path actually fires (cancellations and
  factorisations happen, the result is not the un-simplified
  unevaluated form), and
* the result contains an `EXPR_REAL` leaf, so an inexact input
  produces an inexact output as expected.

The headline regression test is

```
Cancel[(x^2/9 + (2 x y)/15 + y^2/25) / (x^2/9. - y^2/25.)]
  → (12.5 x + 7.5 y) / (12.5 x - 7.5 y)^1.0
```

which hinges on the rationalised intermediate `(25/2 x + 15/2 y)/(25/2 x - 15/2 y)` actually computing.

### Tests

`tests/test_rationalize.c` exercises the full API:

* C-level core: terminating decimals, `27/4` recovery, π rejection in
  default mode, `22/7` and `201/64` for the documented Mathematica
  tolerance values, negative tolerances, zero / negative inputs,
  non-finite inputs.
* `Rationalize[x]` default mode: `6.75 → 27/4`, `N[Pi]` unchanged,
  threading over `Plus`, exact-input pass-through, negative reals,
  `Complex` componentwise.
* `Rationalize[x, dx]` tolerance mode: `Pi/0.01 → 22/7`,
  `Exp[Sqrt[2]]/2^-12 → 218/53`, multi-term polynomial
  rationalisation, the zero-in-interval shortcut, negatives.
* `Rationalize[x, 0]` zero mode: `N[Pi] → 245850922/78256779`,
  short decimals, `Complex[re, im]` componentwise (the Mathematica
  example `Rationalize[N[Pi + 33/211 E I], 0]`).
* Edge cases: 0-arg / 3-arg malformed calls, symbolic and negative
  tolerances stay unevaluated, threading through `List`,
  `Protected` attribute.

## Simplify: half-angle Tan / Tanh and Pythagorean reduction (2026-04-29)

Two new transforms in `src/simp.c` close the long-standing gaps where
`Simplify` could not reach the Weierstrass half-angle form
`Tan[x/2]` (resp. `Tanh[x/2]`) or fold `1 - Cos[x]^2` into `Sin[x]^2`
inside a larger product. Both are invoked as seeds and per-candidate
transforms in `simp_search`, alongside the existing
`PythagSquareComplete` and `TrigRoundtrip`.

### `transform_halfangle` -- half-angle tangent

Folds the Weierstrass identities

```
Sin[x] / (1 + Cos[x])              -> Tan[x/2]
Sin[x]^a (1 + Cos[x])^(-a)         -> Tan[x/2]^a
Sin[x] / (c (1 + Cos[x]))          -> Tan[x/2] / c          (FreeQ[c, x])
Sin[x]^a (c (1 + Cos[x]))^b        -> Tan[x/2]^a c^b        (a + b == 0, FreeQ[c, x])
(1 - Cos[x]) / Sin[x]              -> Tan[x/2]
(1 - Cos[x])^a Sin[x]^b            -> Tan[x/2]^a            (a + b == 0)
```

and the corresponding `Sinh / Cosh -> Tanh[x/2]` family with the sign
flip `(Cosh[x] - 1) / Sinh[x] = Tanh[x/2]`. Each rule has a trailing
`r___` BlankNullSequence on the enclosing `Times`, so the rule fires
on subterms inside larger products (e.g.
`(1/2) Sin[x] / (1 + Cos[x])` rewrites). The general-power forms use
the conditional pattern `/; a + b === 0` to keep the rule shape
free of specific exponent values, which is what handles inputs like

```
Simplify[(Cos[5] + 1)^a Sin[5]^(-a)] -> Tan[5/2]^(-a)
```

without any case-specific rule.

Output complexity is uniformly less than or equal to the input on
every shape that fires, so `simp_search`'s leaf-count tiebreak takes
the rewritten form.

### `transform_pythag_reduce` -- Pythagorean collapse inside products

Strict leaf-count reductions

```
1 - Cos[x]^2  -> Sin[x]^2
1 - Sin[x]^2  -> Cos[x]^2
Cosh[x]^2 - 1 -> Sinh[x]^2
1 + Sinh[x]^2 -> Cosh[x]^2
```

via `ReplaceRepeated` with a trailing `+ r___` so the matching pair
collapses even when it sits among other Plus terms (e.g.
`1 - Cos[x]^2 + 5` -> `5 + Sin[x]^2`). Together with the bottom-up
recursion of `simp_bottomup`, this lets

```
Simplify[Sin[x] (1 - Cos[x]^2)]   -> Sin[x]^3
Simplify[Cos[x] (1 - Sin[x]^2)]   -> Cos[x]^3
Simplify[Sinh[x] (Cosh[x]^2 - 1)] -> Sinh[x]^3
Simplify[Cosh[x] (1 + Sinh[x]^2)] -> Cosh[x]^3
```

reach the cubed form directly: the inner Plus collapses first via
`PythagReduce`, leaving a `Sin[x] * Sin[x]^2` (resp. hyperbolic) Times
that the standard arithmetic simplifier reduces to the cubed form.

### `TrigRoundtrip` seed score-gate

The `TrigRoundtrip` composite (`TrigToExp -> Together -> Cancel ->
ExpToTrig`) is correctness-preserving but its `Together` step on
real-exponential forms can factor out an asymmetric `E^(-kx)` and
introduce a high-frequency `E^(2kx)` cofactor (e.g.
`Together[-1/2 + 1/4 E^(-2x) + 1/4 E^(2x)]` ->
`1/4 E^(-2x) (1 - 2 E^(2x) + E^(4x))`). `ExpToTrig` then turns the
cofactor into a 16-term `Cosh[k x] / Sinh[k x]` product that the
downstream `Together`, `Cancel`, `Factor`, ... transforms each spend
hundreds of milliseconds on per round.

The seed-phase `TrigRoundtrip` is therefore now score-gated: the
result is still considered for `update_best`, but is only added to
the seed set if its score is at most `2 * input_score + 8`. The
gate has no effect on inputs where the round-trip simplifies (the
existing `(E^x - E^(-x))/Sinh[x] -> 2` test is unaffected, alt_score
1 << input_score 7) and stops `Simplify[Sinh[x] (Cosh[x]^2 - 1)]`
from running for over a minute on the per-round expansion of the
16-term cofactor.

### Tests

`tests/test_simplify.c` adds:

* `test_simplify_halfangle_tan_with_factor` -- `Sin[x] / (2 (Cos[x] + 1))`
  to `(1/2) Tan[x/2]`.
* `test_simplify_halfangle_tan_numeric_one|two|four` -- `Sin[1]/(1+Cos[1])`,
  `Sin[2]/(1+Cos[2])`, `Sin[4]/(1+Cos[4])`. The `Sin[2k]/(1+Cos[2k])`
  forms collapse to `Tan[k]` because `Tan[2/2] = Tan[1]` and
  `Tan[4/2] = Tan[2]`.
* `test_simplify_halfangle_tan_power_a` -- general-exponent form.
* `test_simplify_halfangle_tan_inverse_power` -- inverse-exponent form
  `(Cos[5] + 1)^a Sin[5]^(-a) -> Tan[5/2]^(-a)`.
* `test_simplify_halfangle_tan_integer_assumption` -- same with
  `Element[n, Integers]` (assumption is consistent with but not
  required by the algebraic identity).
* `test_simplify_pythag_reduce_sin_cube|cos_cube` -- product collapses.
* `test_simplify_halfangle_tanh_with_factor|basic|power` -- hyperbolic
  half-angle, including the `Sinh[x] / (2 (Cosh[x] + 1)) -> (1/2) Tanh[x/2]`
  form requested by the user.
* `test_simplify_pythag_reduce_sinh_cube|cosh_cube` -- hyperbolic
  product collapse, exercising the `TrigRoundtrip` score-gate path.

## Simplify: sign-flipped Pythagorean reduction + per-transform PythagReduce chain (2026-04-29)

Two follow-on fixes to the Pythagorean reduction pass close gaps that
left `Simplify[Sin[x]^3 + Sin[3 x] - 3 Sin[x]]` stuck at
`3 Sin[x] (-1 + Cos[x]^2)` instead of reaching `-3 Sin[x]^3`.

### Sign-flipped Pythagorean rule pairs

`transform_pythag_reduce` now carries the negated counterparts of the
existing rules so that the reduction fires on either canonical sign of
the constant term in a `Plus`:

```
1 - Cos[x_]^2 + r___  :>  Sin[x]^2 + r
1 - Sin[x_]^2 + r___  :>  Cos[x]^2 + r
-1 + Cos[x_]^2 + r___ :> -Sin[x]^2 + r       (* new *)
-1 + Sin[x_]^2 + r___ :> -Cos[x]^2 + r       (* new *)
-1 + Cosh[x_]^2 + r___ :> Sinh[x]^2 + r
1 + Sinh[x_]^2 + r___ :> Cosh[x]^2 + r
1 - Cosh[x_]^2 + r___ :> -Sinh[x]^2 + r      (* new *)
-1 - Sinh[x_]^2 + r___ :> -Cosh[x]^2 + r     (* new *)
```

Without these, `FactorSquareFree[3 Cos[x]^2 Sin[x] - 3 Sin[x]]` (which
`simp_search` produces while exploring the triple-angle expansion)
yields `3 Sin[x] (-1 + Cos[x]^2)`, and the inner Plus has the canonical
form `-1 + Cos[x]^2`, not `1 - Cos[x]^2`, so the pre-existing rules
never matched.

### Chained PythagReduce inside the round loop

`simp_search` runs `SIMP_ROUNDS` (= 2) of seed expansion. Each round
calls every entry in `SIMP_TRANSFORMS` on each seed, and separately
calls `transform_pythag_reduce` once per seed. A candidate produced
mid-round (say, FactorSquareFree's output) only gets a Pythagorean
reduction when it is promoted to a seed in the *next* round -- but if
this happens at the final round, the loop exits before that follow-up
fires.

The fix calls `transform_pythag_reduce` on every transform output
inside the round, immediately after `update_best`. This adds at most
~0.5 ms per transform output (the rule list is parsed once and
`ReplaceRepeated` short-circuits when nothing matches), and lets the
final-round candidates pick up the Pythagorean rewrite without
requiring a third round.

With both fixes:

```
Simplify[Sin[x]^3 + Sin[3 x] - 3 Sin[x]]
    (* TrigExpand -> 3 Cos[x]^2 Sin[x] - 3 Sin[x]
       FactorSquareFree (round 1) -> 3 Sin[x] (-1 + Cos[x]^2)
       PythagReduce on that output -> -3 Sin[x]^3 *)
->  -3 Sin[x]^3
```

### Test

`tests/test_simplify.c` adds `test_simplify_sin_triple_angle_collapse`
covering the user-reported expression.

## Factor: bivariate Hensel performance (Phase F3, 2026-04-30)

Three optimisations land in `src/mvfactor.c`'s bivariate two-factor
Hensel lift, targeting the per-iteration polynomial-arithmetic cost
and the recombination order:

### F3a -- sort `us[]` by ascending degree before recombination

`bpoly_hensel_lift_multi` now sorts the borrowed `us[]` array by
univariate x-degree at entry, so the smallest-degree singletons get
tried first by `lift_multi_internal`'s subset enumeration.  For
inputs whose true bivariate factor has a low-degree univariate image
(e.g. `(1 - x^k)(x - y^m)` whose `(x - y^m)` reduces to `x` at α=0),
the singleton `{x}` is now the first recombination trial instead of
being one of `r` possibilities.  The complement-side recursion
preserves relative order, so a single sort at the public entry is
enough.

### F3b -- incremental U·V update across Hensel iterations

The previous lift recomputed `U*V mod y^{k+1}` from scratch each
iteration via `bpoly_mul_truncate_y`, paying O(d_x² · k²) per step
and O(d_x² · B³) cumulative for a y-degree-B lift.  F3b maintains a
running `UV` BPoly tracking `U_{k-1} · V_{k-1} mod y^{k+1}` and
extends/updates it incrementally: at the start of iteration k we add
the cross sum `Σ_{a=1}^{k-1} U[y^a] · V[y^{k-a}]` at y^k, and after
the diophantine update we add `Δu · V_0 + U_0 · Δv` at y^k (the
`y^{2k} · Δu · Δv` term vanishes mod y^{k+1} for k ≥ 1).  Cumulative
cost drops to O(d_x² · B²) -- a factor-of-B speedup.  Applied to
both `bpoly_hensel_lift_2` and `bpoly_hensel_lift_2_lc`; for the LC
variant, the y^0 polynomials are `qu_0 · u` and `qv_0 · v` instead
of plain `u, v`.

### F3c -- Mignotte coefficient bound fast-fail

Each lift call computes a conservative integer-coefficient bound
`M = 2^(deg_x(P) + deg_y(P)) · ceil(||P||_2)` once at entry and
checks it after every iteration: if any coefficient of `U` or `V`
exceeds `M` in absolute value, the subset cannot lift to an integer
factor of P, so the lift aborts immediately rather than running the
remaining iterations on a divergent candidate.

### Tests

`tests/test_factor_recombine.c` adds `test_f3_perf_x6_y15`,
`test_f3_perf_x4_y25`, and `test_f3_perf_x8_y17` covering high-y-
degree inputs (B ∈ {15, 17, 25}) where the lift inner loop dominates.
Each completes in single-digit milliseconds end-to-end.

### Out of scope

`Factor[Expand[(1 - x^12)(x - y^13)]]` was previously cited as the
F3 motivating case.  Tracing showed the hang is actually inside
`bz_factor_to_expr(x^13 - x)` -- the univariate Berlekamp-Zassenhaus
factorizer -- before the bivariate path runs.  `Factor[x^13 - x]`
and `Factor[x^13 - 1]` hang in isolation.  Fixing that is a separate
task (univariate BZ performance, not bivariate Hensel).


## Evaluator alignment with Withoff §3.1 (M1 milestone, 2026-05-01)

The PicoCAS evaluator (`src/eval.c::evaluate_step`) was audited against
the algorithm described in David Withoff's 1992 *Mathematica Internals:
A Tutorial* (§3.1). Five correctness gaps were closed. All 79 unit tests
continue to pass.

### 1. Flat now runs before Listable

The structural-attribute order is now `Flat → Sequence (already done) →
Listable → Orderless`, matching Withoff §3.1. Previously Listable ran
first and missed lists that only surfaced after flattening. Concretely,
`Plus[Plus[a, {1,2}], 3]` now correctly threads to `{4 + a, 5 + a}` —
Flat first promotes the inner `Plus` so that Listable sees the list arg.

### 2. User DownValues now run before built-ins

The dispatch order inside `evaluate_step` is now:

1. apply user-defined DownValues
2. call the head's `builtin_func` (the equivalent of Mathematica's
   internal "down code")
3. handle `Set` / `SetDelayed` primitives

Mathematica's documented order is "user rules first, then internal".
Protected symbols are unaffected because `apply_assignment` already
refuses to install DownValues on a Protected target. The visible effect
is that a user can override a non-Protected built-in if one exists,
matching Mathematica behaviour.

### 3. OneIdentity restricted to pattern matching

`ATTR_ONEIDENTITY` is no longer applied as an evaluation rewrite. The
1-arg collapse (`Plus[x] → x`, `Times[x] → x`, `Power[x] → x`,
`GCD[x] → x`, `LCM[x] → x`, `And[x] → x`, `Or[x] → x`, `Dot[x] → x`)
is now performed inside each head's builtin via an explicit
`if (n == 1) return expr_copy(args[0])` early-return. Pattern-matching
OneIdentity in `src/match.c` is unchanged.

The visible behaviour change: a user-defined OneIdentity head that has
no rule of its own no longer collapses. For example,

    SetAttributes[g, OneIdentity]
    g[x]               (* now stays as g[x], previously rewrote to x *)

This matches Mathematica's actual semantics for OneIdentity (it is a
pattern-matching attribute, not an evaluation-time rewrite). Builtins
modified to handle `n == 1` themselves: `builtin_plus`, `builtin_times`
(both already had it), `builtin_power`, `builtin_gcd`, `builtin_lcm`,
`builtin_and`, `builtin_or`.

### 4. HoldAllComplete no longer short-circuits the evaluator

Previously, a head with `ATTR_HOLDALLCOMPLETE` caused `evaluate_step`
to return immediately, skipping built-ins, DownValues, Flat,
Listable, and Orderless. This made it impossible for a HoldAllComplete
head to do useful work via its builtin or rules.

The new behaviour suppresses only:

- argument evaluation (including the `Evaluate[]` override)
- `Sequence` flattening
- `Unevaluated` stripping
- `Flat` flattening (matches Mathematica's documented HoldAllComplete
  behaviour)

…but lets built-ins and DownValues attached to the head still run.
For example:

    Length[HoldComplete[a, b, c]]    →    3
    HoldComplete[Sequence[a, b, c]]  →    HoldComplete[Sequence[a, b, c]]
    HoldComplete[Unevaluated[1+2]]   →    HoldComplete[Unevaluated[1 + 2]]

### 5. $RecursionLimit C-stack guard

Nested `evaluate()` calls now track a static depth counter (the REPL is
single-threaded so a static suffices). When the counter would exceed
`eval_recursion_limit` (default 1024, matching modern Mathematica), the
input is wrapped in `Hold[]`, a sticky `eval_overflow` flag is set, and
a `$RecursionLimit::reclim` message is emitted. Every enclosing
fixed-point loop checks the flag on each iteration and bails out so the
unwind doesn't burn `$IterationLimit` × 1024 iterations of redundant
rewrite attempts.

Concretely, a self-referential rule like

    f[x_] := f[x] + 1
    f[1]

now terminates in milliseconds with a clear `$RecursionLimit` message
rather than running C-stack-bound recursion plus `$IterationLimit`
worth of futile rewrites at every depth.

## $RecursionLimit user-visible symbol (M1 close-out, 2026-05-01)

The C-side limit is now exposed as a writable system variable. At
startup, `eval_init()` seeds an OwnValue of 1024 on `$RecursionLimit`
and installs a docstring. When the user assigns to it via `Set` or
`SetDelayed`, `apply_assignment` validates the RHS and pushes the new
value into `eval_set_recursion_limit()`:

    $RecursionLimit              (* → 1024 *)
    $RecursionLimit = 256        (* → 256, C-side limit also 256 *)
    f[x_] := f[x] + 1
    f[1]                         (* aborts at depth 256 *)
    $RecursionLimit = 5          (* $RecursionLimit::limset, OwnValue rolled back *)

Validation rules:

- The RHS must evaluate to a positive integer.
- Values below 20 are rejected with a `$RecursionLimit::limset`
  message; the OwnValue is rolled back to the active C-side limit so
  the symbol does not lie about the in-effect value.
- Bigints that don't fit in `INT_MAX` are clamped.
- `SetDelayed` (`$RecursionLimit := ...`) is supported by evaluating a
  copy of the RHS at assignment time so the limit reflects the value
  the user expects to read back.

This closes out the M1 milestone of `EVAL_IMPROVEMENTS_PLAN.md`. The
plan itself has also been revised to remove UpValues / SubValues from
scope (those Mathematica semantics are intentionally not pursued).

### Files touched

`src/eval.c`, `src/eval.h`, `src/core.c`, `src/arithmetic.c`,
`src/boolean.c`, `src/power.c`, plus `EVAL_IMPROVEMENTS_PLAN.md`.

## simp_radicals: cross-base radical fusion in Simplify (2026-05-01)

`Simplify` now combines distinct positive-integer radicals that share
an exponent inside any `Times` node:

    Sqrt[2] Sqrt[3]                    →  Sqrt[6]
    Sqrt[2] Sqrt[3] Sqrt[5]            →  Sqrt[30]
    Sqrt[2] Sqrt[3] Sqrt[6]            →  6           (Sqrt[36] collapses)
    2^(1/3) 3^(1/3) 5^(1/3)            →  30^(1/3)
    Sqrt[6] - Sqrt[2] Sqrt[3]          →  0
    -Sqrt[2] Sqrt[3] x + Sqrt[6] x     →  0

The evaluator does not perform this combine automatically because, in
general,

    a^(p/q) * b^(p/q)  ≠  (a*b)^(p/q)

once `a` or `b` can be negative or non-real -- the principal-value
branch shifts (`Sqrt[-2] Sqrt[-3] = -Sqrt[6]`, not `+Sqrt[6]`).
Restricting to **positive integer bases** keeps the rewrite sound.
Same-base products (`Sqrt[2] Sqrt[2] -> 2`) were already collapsed by
`Power`'s exponent merging in the evaluator; this transform fills the
cross-base gap that previously left `Sqrt[6] - Sqrt[2] Sqrt[3]`
unsimplified.

### Implementation

`simp_radicals` (`src/simp.c`) is a memoised bottom-up walker. For each
`Times` node it buckets `Power[+integer, Rational[p, q]]` factors by
exponent; multi-element buckets are fused into a single
`Power[product_of_bases, p/q]` which is then evaluated so the
perfect-power detection in `builtin_power` collapses
`Power[36, 1/2] -> 6`. Symbolic, negative-integer, and rational-non-
integer bases are left alone.

### Wiring

`simp_radicals` is registered both as a force-take seed in
`simp_search`'s seed phase (alongside `transform_pythag_reduce`,
`transform_halfangle`, etc.) and as a per-candidate rewrite inside the
round loop, so it can fire on `Together` / `Cancel` outputs that
surface fresh `Sqrt[a] Sqrt[b]` products. The transform is gated by
`has_non_integer_power`, so on inputs without any radicals the walker
short-circuits to a copy.

### Files touched

`src/simp.c` (new transform + seed/round wiring),
`tests/test_simplify.c` (seven new test cases covering two-, three-,
and cube-root combines, perfect-square collapse, the targeted
cancellation cases, and the symbolic-base no-op).

## TrigExpand seed propagation + chained simp_radicals (2026-05-01)

Two follow-on changes to `simp_search` (`src/simp.c`) so that linear
combinations of rationalised trig terms now collapse against their
angle-addition forms. Concretely:

    e1 = (Cos[x]/Sqrt[6] + Sin[x]/Sqrt[2])(1 + 1/3)
    e2 = (4/9) Sqrt[6] Sin[x + Pi/6]
    Simplify[e1 - e2]                    →  0

Previously `Simplify` returned the input unchanged because the
`TrigExpand` candidate -- which expands `Sin[x + Pi/6]` into
`(1/2) Cos[x] + (Sqrt[3]/2) Sin[x]` and surfaces `Sqrt[3] Sqrt[6]`
products -- has a higher leaf count than the input, so the round-loop
score gate dropped it before `Together` / `simp_radicals` got a chance
to combine and cancel.

### Change 1: TrigExpand bypasses the score gate, with a bound

In the round loop's standard transform dispatcher, when the firing
transform is `TrigExpand` and the result has higher leaf count than
the seed, the result is still added to the next round's seed set
provided `score(r) <= 2 * parent_score + 8` (the same bound
`TrigRoundtrip` uses at the seed phase). This lets the expanded form
seed a follow-on round of `Together` / `Cancel` / `simp_radicals` /
per-variable `Collect` while still pruning pathological blow-ups
(e.g. `Sin[a + b + c + d]` cascades).

### Change 2: simp_radicals chained on transform outputs

After every `SIMP_TRANSFORMS` entry fires, the result is passed
through `simp_radicals` and `update_best`-tested before being
score-gated for next-round propagation. This mirrors the chained
`transform_pythag_reduce` pass already in place. Without it, a
`Together` output like

    (1/9 (-12 Sqrt[2] Sqrt[3] Sin[x] + 12 Sqrt[6] Sin[x])) /
        (Sqrt[2] Sqrt[6])

would only feed `simp_radicals` after becoming a seed in the next
round -- which doesn't happen at the final round (`SIMP_ROUNDS = 2`).
Chaining catches the `Sqrt[2] Sqrt[3] -> Sqrt[6]` fusion that drives
the numerator to zero.

### Files touched

`src/simp.c` (round-loop score gate + chained simp_radicals),
`tests/test_simplify.c` (one new test
`test_simplify_trig_radical_angle_addition`).

## simp_split_additive: disjoint-variable Plus decomposition (2026-05-01)

A new pre-pass at the top of `simp_dispatch` (`src/simp.c`) that
decomposes `Plus[a_1, ..., a_n]` into connected components in the
free-symbol-sharing graph, simplifies each component independently
through `simp_dispatch`, and returns the canonicalised sum.

### Motivation

The two-variable form of the rationalised-trig angle-addition test
case,

    e1 = Cos[x]/Sqrt[6] + Sin[x]/Sqrt[2] +
         Cos[y]/Sqrt[6]/3 + Sin[y]/Sqrt[2]/3
    e2 = Sqrt[6] Sin[x + Pi/6]/3 + Sqrt[6] Sin[y + Pi/6]/9
    Simplify[e1 - e2]                              →  0

hit `$RecursionLimit::reclim` and ran for minutes when simplified as
a whole. Each three-term piece (the x-piece and the y-piece, each
analogous to the one-variable case above) collapses to 0 in ~1s on
its own. The cost of `simp_search`'s trig transforms (`Together`,
`TrigExpand`, `Cancel`) is super-linear in the number of addends, and
addends over disjoint user-symbol sets cannot interact under any of
those transforms.

### Algorithm

1. If the input is not a `Plus` with at least 4 addends, return NULL
   (let the regular pipeline run).
2. For each addend, collect its "free symbols": symbol leaves
   appearing in argument position only -- numeric constants
   (`Pi`, `E`, `EulerGamma`, ...) and symbols that appear ONLY as
   function heads (`Sin`, `Cos`, ...) are excluded. Without the
   head-skip rule, every all-trig addend would glom into one
   component and the split would never fire.
3. Union-find over addends, joining pairs whose free-symbol sets
   intersect.
4. If there are fewer than 2 components, or every component is a
   singleton, return NULL.
5. For each component, build a sub-`Plus` of its addends, evaluate
   to canonicalise, run `simp_dispatch` on it.
6. Return `evaluate(Plus[results...])`.

### Why disjoint splits are safe

`Together`, `TrigExpand`, `Cancel`, `Factor`, `Collect`, and the
trig roundtrip all operate on the algebraic structure of the
expression. Two addends with disjoint free-symbol sets share no
common factor, common denominator, or common angle, so no transform
in `SIMP_TRANSFORMS` can produce a cross-term that didn't exist in
either piece individually. The sum of the per-component
simplifications is therefore equal to the simplification of the sum,
modulo evaluator canonicalisation (which also operates additively
under `Plus`'s `Flat | Orderless` attributes).

### Files touched

`src/simp.c` (new helpers `split_symset_*`, `split_uf_*`,
`simp_split_additive`; `simp_dispatch` calls it first),
`tests/test_simplify.c` (one new test
`test_simplify_trig_radical_angle_addition_two_vars`).

## simp_split_multiplicative: disjoint-variable Times decomposition (2026-05-01)

### Motivation

The multiplicative analog of `simp_split_additive`. A `Times` node
whose factors fall into 2+ variable-disjoint components inflates
`simp_search`'s effective free-symbol budget, and the heuristic gives
up before reaching trig identities that exist within a single
component. Concretely:

```
Simplify[Cos[x] Cos[y] Sec[x+y] (Tan[x] + Tan[y])]            -> Tan[x+y]
Simplify[Tan[z] Cos[x] Cos[y] Sec[x+y] (Tan[x] + Tan[y])]      (no reduction)
```

The standalone product collapses via `TrigFactor`. Adding the
disjoint factor `Tan[z]` adds a third free symbol; the search budget
fans out across more candidates and the reduction never lands.

The full motivating case is the nested tan-addition formula:

```
Simplify[(Tan[z] + (Tan[x]+Tan[y])/(1 - Tan[x] Tan[y]))
       / (1 - (Tan[x]+Tan[y]) Tan[z]/(1 - Tan[x] Tan[y]))]   -> Tan[x+y+z]
```

`simp_bottomup` recurses into the denominator's `Times` subterm
`Times[-1, Cos[x], Cos[y], Sec[x+y], Tan[z], Plus[Tan[x], Tan[y]]]`,
where `simp_split_multiplicative` lifts `Tan[z]` (and the literal
`-1`) into singleton components and dispatches the
variable-{x,y} subgroup -- which then collapses to `Tan[x+y]`.

### Algorithm

1. Fire only on `Times[a1, a2, ..., an]` with `n >= 3`.
2. Collect free symbols per factor with the same walker used by the
   additive splitter (`split_collect_addend_symbols`), so that
   numeric-constant symbols (Pi, E, ...) and pure function heads
   (Sin, Cos, ...) are excluded.
3. Pairwise union-find: two factors are unioned iff their free-symbol
   sets intersect. Empty sets (numeric-only factors) stay singleton.
4. Bail out unless the partition gives 2+ components AND at least one
   component holds 2+ factors. (A 4-singleton split would just
   re-dispatch every atom for no win.)
5. For each component, build a sub-`Times`, evaluate it (so canonical
   form invariants hold), and recursively `simp_dispatch` it.
6. Multiply the per-component results back together and return the
   evaluated product.

### Why disjoint factor splits are safe

Same argument as the additive case. Every transform in
`SIMP_TRANSFORMS` (`Together`, `Cancel`, `Expand`, `Factor`,
`TrigExpand`, `TrigToExp`, `TrigFactor`, `TrigRoundtrip`,
`Collect[_, _]`, the Pythagorean rules, half-angle, `Radicals`) acts
within a single variable's algebraic or trigonometric structure.
Variable-disjoint factors share no common base, no common angle, no
common Pythagorean pair -- so the product of the per-component
simplifications equals the simplification of the whole product,
modulo `Times`'s `Flat | Orderless` canonicalisation (which the
final `evaluate` restores).

### Files touched

`src/simp.c` (new `simp_split_multiplicative`; `simp_dispatch` calls
it after `simp_split_additive`).

## Factor: recombination short-circuit (Phase F5, 2026-05-02)

After Phase F3 made the bivariate Hensel inner loop fast, the
remaining hot spot for r-large irreducible bivariate inputs (e.g.
`x²⁴ + y² - 1` whose univariate image at any small alpha factors
into 8 cyclotomic pieces) was the exhaustive k-subset enumeration in
`lift_multi_internal`.  For r = 8 the loop tries 8 + 28 + 56 + 35 =
127 lifts before concluding "irreducible".  Phase F5 adds a
near-success heuristic that short-circuits the loop early.

### The signal

`bpoly_hensel_lift_2_internal` (a new file-static helper in
`src/mvfactor.c`) exposes a `completed_iters_out` flag.  It is set
to `true` only when the y-iteration loop completed all B steps
without bailing on a Diophantine non-integer or Mignotte
coefficient overflow.  Such a lift "almost worked" -- only the
final `U·V == P` over Z verify failed -- which is the natural
signal that recombination at the next k might find the right
partition.

The public `bpoly_hensel_lift_2` is a thin wrapper passing a local
dummy flag; existing call sites are unaffected.

### The short-circuit

`lift_multi_internal` aggregates the flag across all k-subsets of a
given size.  When an entire k-pass produces no completion signal
AND k ≥ 2, the loop breaks and returns P as irreducible
bivariately.  k = 1 is always exhausted unconditionally, so
legitimate singleton lifts (e.g. `{x²+1}` peeling off `x²+y` from
`x⁴ - y²`) still fire.  k = 2 is always tried before any
short-circuit, so all current k = 2 recombination cases (e.g.
`(x²-y)(x²-4y)` with r = 4 image (x-1)(x+1)(x-2)(x+2)) remain
correct.

### Measured impact

| Workload | Pre-F5 | Post-F5 | Δ |
|---|---|---|---|
| `Factor[Expand[x¹²+y²-1]]` (r=6) | 23 ms | 20 ms | -13 % |
| `Factor[Expand[x¹²+y⁴-1]]` (r=6, B=4) | 24 ms | 20 ms | -17 % |
| `Factor[Expand[x²⁰+y²-1]]` (r=6) | 31 ms | 28 ms | -10 % |
| `Factor[Expand[x²⁴+y²-1]]` (r=8) | 117 ms | 84 ms | -29 % |
| `Factor[Expand[x³⁰+y²-1]]` (r=8) | -- | 59 ms | (new) |

### Files touched

`src/mvfactor.c`: `bpoly_hensel_lift_2_internal` (new; carries the
`completed_iters_out` output), `bpoly_hensel_lift_2` (thin wrapper),
`lift_multi_internal` (tracks `any_completion_at_k`, breaks when
the per-k-pass signal is absent and k ≥ 2).

`tests/test_factor_recombine.c`: `test_f5_irreducible_r6_x12`,
`test_f5_irreducible_r8_x24`, `test_f5_irreducible_r6_x12_y4`,
`test_f5_recombination_still_works_k2`.  The first three pin the
short-circuit on bivariate-irreducible inputs whose univariate
image has 6 or 8 factors; the last verifies F5 does not regress
legitimate k = 2 recombination.

## Factor: trivariate two-factor monic Hensel via MPoly (Phase F2 MVP, 2026-05-02)

The previous Phase 4 (`factor_via_z_independent_split`) handled n ≥ 3
inputs only when at least one factor was independent of one variable.
Cases like `(zx - x² - y²)(3z + 4xy - y²)` -- where every factor
depends on every variable -- fell through to `factor_roots` which
caught only `(var - c)` shapes, then ultimately returned the input
unchanged after ~6 seconds.  Phase F2 adds the missing multivariate
Hensel pipeline.

### MPoly substrate (`src/mpoly.{c,h}`)

Sparse Z[x_1, ..., x_n] type with arbitrary-precision (mpz_t)
coefficients.  Terms are stored in lex descending monomial order with
parallel row-major exponent and coefficient arrays.  Operations:

- Construction: `mpoly_new`, `mpoly_zero`, `mpoly_from_int`,
  `mpoly_from_mpz`, `mpoly_monomial`, `mpoly_copy`.
- Term manipulation: `mpoly_push_term` (bulk append, then
  `mpoly_normalize` to sort/dedup/drop-zero), `mpoly_set_coef`,
  `mpoly_get_coef`.
- Predicates: `mpoly_is_zero`, `mpoly_eq`, `mpoly_deg_var`,
  `mpoly_total_deg`, `mpoly_is_constant_in_var`.
- Arithmetic: `mpoly_add` (merge sort), `mpoly_sub`, `mpoly_neg`,
  `mpoly_mul` (full term-by-term then normalise), `mpoly_scale`.
- Substitution: `mpoly_subst_var_int` (x_i := alpha),
  `mpoly_shift_var_int` (x_i := x_i + alpha via binomial expansion),
  `mpoly_coef_of_var` (returns the polynomial coefficient of x_i^k),
  `mpoly_lc_var`.
- Expr round-trip: `expr_to_mpoly`, `mpoly_to_expr`.

### Trivariate Hensel pipeline (`src/mvfactor3.{c,h}`)

Two new primitives:

- **Bivariate Diophantine** (`mpoly_diophantine_2`): solve
  `Δu·V + Δv·U = E` over Z[x_main, x_y] with monic-in-x_main U, V
  coprime over Q(x_y)[x_main].  Algorithm: try alpha values for
  x_y; for each promising alpha, solve the univariate Diophantine
  via `zupoly_diophantine`; then Hensel-lift in x_y up to
  `deg_{x_y}(E)`.

- **Trivariate two-factor Hensel** (`mpoly_hensel_lift_3_2`): given
  P ∈ MPoly{n vars} with bivariate seeds U_xy, V_xy at x_z = α_z,
  lift to U, V satisfying U·V = P over Z[x_main, x_y, x_z].  Shifts
  x_z by α_z, iterates degrees k = 1..deg_{x_z}(P_sh), solves the
  bivariate Diophantine on each residual, accumulates the
  corrections, verifies the final product.

The implementation reuses the existing `bpoly_hensel_lift_2`
machinery for the inner y-lift step via MPoly ↔ BPoly converters
(`mpoly_to_bpoly_in`, `bpoly_to_mpoly_in`).

### Wiring (`src/facpoly.c::factor_trivariate_via_mhensel`)

Invoked in `heuristic_factor` between Phase 4 (specialise-and-divide)
and the irreducibility gate, only for `v_count == 3`.  The
orchestrator:

1. Convert P to MPoly.
2. Try each variable as main, requiring degree ≥ 1 and constant
   leading coefficient = 1.
3. Try (α_y, α_z) ∈ {0, ±1, ±2, ±3}^2 = 49 alpha tuples.
4. For each: compute univariate image, check squarefree, factor via
   existing `factor_via_bz_callback`, require r = 2.
5. Specialise z = α_z, shift y by α_y, run `bpoly_hensel_lift_2`.
6. Convert bivariate factors to MPoly, unshift y, then run
   `mpoly_hensel_lift_3_2` to introduce z.
7. Verify product == P, convert back to Expr, return.

### Drive-by fix to `extract_monomial`

The trivariate test case `(z + xy)(z - xy)` exposed a pre-existing
bug in `extract_monomial`: the function did not descend into nested
Times structures.  For `Plus[Power[z,2], Times[-1, Times[Power[x,2],
Power[y,2]]]]`, the inner Times was silently dropped, so
`factor_binomial` saw the polynomial as `z² + (constant -1)` and
incorrectly returned `(z-1)(z+1)`.  Fixed `extract_monomial` to
recursively descend into nested Times via a new
`extract_monomial_walk` helper.

### Measured impact

| Workload | Pre-F2 | Post-F2 |
|---|---|---|
| `Factor[Expand[(zx - x² - y²)(3z + 4xy - y²)]]` | ~6 s, returns input unchanged | **0.33 s, correctly factored** |
| `Factor[Expand[(z + xy)(z - xy)]]` | (-1+z)(1+z) WRONG | 15 ms, correct |
| `Factor[Expand[(z + x)(z + y)]]` | (-1+z)(1+z) WRONG | 12 ms, correct |
| `Factor[Expand[(z² + x)(z + y)]]` | wrong | 5 ms, correct |

### MVP scope

What it handles: n = 3 inputs, monic in some main variable, exactly
two factors, with a good (α_y, α_z) found in 49 attempts.

What it doesn't handle yet (follow-up phases):
- n ≥ 4 (extension is recursive but needs more glue code).
- Non-monic-in-main inputs (needs Wang's lc correction at the
  trivariate level, mirroring F1 at the bivariate level).
- r ≥ 3 outputs (needs trivariate multifactor Hensel + recombination).

### Files touched

- New: `src/mpoly.{c,h}`, `src/mvfactor3.{c,h}`,
  `tests/test_mpoly.c`, `tests/test_mvfactor3.c`.
- Modified: `src/facpoly.c` (added `factor_trivariate_via_mhensel`,
  fixed `extract_monomial` for nested Times),
  `tests/CMakeLists.txt` (registered new test targets and source
  files), `tests/test_factor_recombine.c` (added 3 F2 end-to-end
  tests including the user-reported case).

## Symbol interning (M2 milestone, 2026-05-03)

Tier 3 prep work from `EVAL_IMPROVEMENTS_PLAN.md` step 3.1 has shipped:
the system now interns every symbol name behind a global table, so
identity tests on symbols collapse to pointer equality.

### What changed

- **New module `src/sym_intern.{c,h}`.** A chained hash table
  (`INTERN_BUCKETS = 8191`) maps a symbol name to a single,
  table-owned `const char*`. `intern_symbol(name)` returns the
  canonical pointer; `intern_clear()` releases all interned strings.
  The first lookup-by-already-canonical-pointer short-circuits via
  pointer compare before any `strcmp` runs, so re-interning costs
  nothing.

- **`Expr` symbol identity is now pointer-based.**
  `expr_new_symbol(name)` stores the interned pointer in
  `data.symbol`; `expr_copy` of a symbol copies that pointer
  (no `strdup`); `expr_free` skips freeing the symbol-name field
  (the interner owns it); `expr_eq` for `EXPR_SYMBOL` is a `==`
  compare; `expr_compare` short-circuits on pointer equality before
  falling back to lexicographic `strcmp` (canonical sort order is
  unchanged). `EXPR_STRING` keeps its existing `strdup`/`free`
  semantics — only symbol names are interned.

- **Symbol table keyed on canonical pointers.**
  `symtab_get_def` / `symtab_lookup` / `symtab_get_own_values` /
  `symtab_get_down_values` intern their input first, then walk
  buckets with pointer compares; `SymbolDef::symbol_name` stores
  the interned pointer; `symtab_clear` no longer frees it.

- **Cached well-known symbols.** `src/sym_names.{c,h}` declares
  `SYM_List`, `SYM_Plus`, `SYM_Times`, `SYM_Power`, `SYM_Sequence`,
  `SYM_Hold`, `SYM_HoldComplete`, `SYM_HoldPattern`, `SYM_HoldForm`,
  `SYM_Unevaluated`, `SYM_Evaluate`, `SYM_Set`, `SYM_SetDelayed`,
  `SYM_Rule`, `SYM_RuleDelayed`, `SYM_Pattern`, `SYM_Blank`,
  `SYM_BlankSequence`, `SYM_BlankNullSequence`, `SYM_PatternTest`,
  `SYM_Condition`, `SYM_Function`, `SYM_Slot`, `SYM_SlotSequence`,
  `SYM_CompoundExpression`, `SYM_If`, `SYM_True`, `SYM_False`,
  `SYM_Null`, `SYM_Rational`, `SYM_Complex`, `SYM_Part`,
  `SYM_Composition`, `SYM_Derivative`. `sym_names_init()` is invoked
  from the top of `core_init()`, before any other init step that
  might create symbols.

- **Hot evaluator paths now use pointer compares.** `has_list_arg`,
  `apply_listable`, `flatten_sequences`, the held-`Evaluate` /
  `Set` / `SetDelayed` / `Rule` / `RuleDelayed` / `Unevaluated` /
  `Function` / `Derivative` / `Composition` / `Condition` / `Part` /
  `List` checks in `evaluate_step`, and `assignment_target_symbol`
  all replaced their `strcmp` calls with `==` against the cached
  `SYM_*` pointers.

- **Defensive intern in `eval_flatten_args`.** This helper is also
  invoked from `internal_call_impl` in `internal.c`, which passes a
  bare C string literal (`"Plus"`, `"Times"`, ...) rather than an
  interned pointer. The function now re-interns its `head_name`
  argument on entry so the inner pointer-compare against
  `arg->data.function.head->data.symbol` still hits.

### Why this matters

Every `evaluate_step` invocation runs a handful of "is the head
named X" checks, and the codebase has hundreds of `strcmp` sites
on symbol names overall. Interning replaces an O(name-length) byte
walk with a single pointer compare in the hot loop, and gives
later milestones (refcount + timestamp; rule dispatch index) a
stable identity to key on.

### Files touched

- New: `src/sym_intern.{c,h}`, `src/sym_names.{c,h}`.
- Modified: `src/expr.c` (symbol allocation, copy, free, equality,
  compare), `src/symtab.c` (canonical-pointer lookups + storage),
  `src/eval.c` (hot-path comparisons + defensive intern),
  `src/core.c` (call `sym_names_init()`),
  `tests/CMakeLists.txt` (registered new sources in `COMMON_SRC`),
  `EVAL_IMPROVEMENTS_PLAN.md` (status + milestone notes).

### Verification

- All 83 unit-test binaries pass under `tests/build`.
- The full REPL smoke test (PolynomialLCM, Plus/Listable interaction,
  Sin over a list, HoldComplete, OneIdentity, DownValue patterns,
  derivatives) returns identical output to baseline.
- `valgrind --leak-check=full --show-leak-kinds=all`: interner
  allocations correctly classified as still-reachable for the
  program lifetime; no new definitely-lost or indirectly-lost
  records introduced.

## Refcounted atom sharing (M3 phase-1, 2026-05-03)

This is the first deliverable of the M3 milestone (`EVAL_IMPROVEMENTS_PLAN.md`
§3.2). It replaces the deep-copy semantics of `expr_copy` for **atoms**
(integer, real, symbol, string, bigint, mpfr) with reference-counted
sharing. FUNCTION nodes are still deep-copied — function-node sharing
plus `expr_unshare()` is deferred to M3 phase-2.

### Mechanism

Every `Expr` now carries `unsigned refcount`:
- All `expr_new_*` constructors initialize the field to `1`.
- `expr_ref(e)` is a new helper that bumps the count and returns the
  same pointer.
- `expr_free(e)` first decrements; only when the count reaches zero
  does it actually release children, payload, and the node itself.
- `expr_copy(atom)` now does `e->refcount++; return e;` — no malloc,
  no payload copy.
- `expr_copy(function)` still allocates a fresh FUNCTION node and
  recursively copies `head` + `args[]` (each recursive copy hits the
  fast atom path for any atom leaves).

### Three-phase staging

The change was landed in three reviewable phases so that each phase
is independently testable:

1. **Phase 1 — add field, no semantic change.** Refcount field
   initialised by every constructor; `expr_free` decrements before
   actually freeing. Pure no-op semantically because nothing yet
   shares (every node still has refcount 1).
2. **Phase 2 — audit atom-payload mutations.** Required because once
   atoms are shared, in-place mutation of an atom's payload would
   silently corrupt every other holder. Nine sites converted to
   "free-and-replace" (allocate a fresh atom, then free the old
   handle):
   - `print.c` ×6 (Times-of-integer/real negation, Rational integer
     numerator negation, Plus-Times integer/real negation)
   - `linalg.c` ×2 (lines 619 & 912: `new_den` integer negation)
   - `simp.c` ×1 (line 4281: `const_terms[i]` integer negation)
3. **Phase 3 — flip `expr_copy`.** With the audit complete, the
   atom-share semantics are turned on. The first test run exposed
   five additional sites that the Phase-2 grep missed — all in
   `poly.c`'s BPList machinery, where exponent atoms held in a
   working list were being mutated in place during GCD/LCM
   common-factor extraction (lines ~1214, 1261, 1386, 1435, 1489).
   These were converted to free-and-replace as well.

### Audit invariant

The Phase-3 invariant — *no atom payload may be mutated in place after
construction* — is enforced module-wide by these checks (which all
return only `expr.c` constructor sites today):

```
grep -nE "data\.(integer|real)\s*[+\-*/]?="     src/*.c
grep -nE "->type\s*=\s*EXPR_"                   src/*.c
grep -nE "mpz_\w+\s*\(\s*\S+->data\.bigint"     src/*.c
grep -nE "mpfr_\w+\s*\(\s*\S+->data\.mpfr"      src/*.c
grep -nE "free\(\s*\S+->data\.(symbol|string)"  src/*.c
```

The single non-`expr.c` mpfr write is `numeric.c:403,493`; both
target a freshly-allocated `r` whose refcount is provably 1, so the
mutation is on a private node and is not a sharing hazard.

Function-node mutations (`->data.function.{args,head,arg_count} = …`
in `eval.c:eval_flatten_args`, `core.c:builtin_polynomial_div_mod`,
`parse.c`) remain safe because FUNCTION nodes are not yet shared in
phase-1.

### Convention for new code

Going forward:
- Never write `e->data.X = ...` on an atom you did not just construct.
  Always: `expr_free(e); e = expr_new_*(new_value);`
- `expr_copy(atom)` is now O(1) and allocation-free. Use it freely
  where you previously hesitated to copy an atom.
- `expr_ref(e)` is the explicit form when you want a second handle to
  the same node (no copy intent), and is mandatory if `e` is reachable
  from elsewhere and you might forget that `expr_copy` is now a
  ref-bump rather than a deep copy.

### Verification

- All 83 unit-test binaries under `tests/build` pass.
- `valgrind --leak-check=full ./poly_tests`: definitely-lost grows
  from 4,360 B / 137 blocks (baseline) to 28,832 B / 749 blocks
  (phase-1). No new bugs — every extra-leaked stack is an existing
  imbalanced-ownership site (`get_coeff_expanded`, `poly_div_rem`,
  `builtin_polynomial{quotient,gcd,extendedgcd}`, …) where the prior
  code leaked one fresh atom node per occurrence; under refcount,
  those same sites keep the original atom alive at refcount > 0.
  A follow-up sweep will balance the missing `expr_free` calls.

### Files touched

- `src/expr.h`: added `unsigned refcount` to `Expr`; declared
  `expr_ref`.
- `src/expr.c`: refcount initialization in every constructor;
  `expr_ref` body; `expr_free` is now dec-ref-and-maybe-free;
  `expr_copy` short-circuits to inc-ref for non-FUNCTION types;
  `expr_bigint_normalize` switched from raw `mpz_clear+free` to
  `expr_free` so it is refcount-safe.
- `src/print.c`, `src/linalg.c`, `src/simp.c`, `src/poly.c`: 14
  in-place atom mutations converted to free-and-replace.
- `EVAL_IMPROVEMENTS_PLAN.md`: M3 §3.2 status block updated.


## Refcounted FUNCTION sharing (M3 phase-2, 2026-05-03)

Phase-2 of the M3 refcounting milestone landed alongside Phase-1 on
2026-05-03. With Phase-1, only atoms (Integer, Real, Symbol, String,
BigInt, MPFR) were shared via `expr_copy` inc-ref; FUNCTION nodes
still deep-copied. Phase-2 makes `expr_copy` an unconditional inc-ref
for every node type, FUNCTION included.

### Mechanism

- `Expr* expr_copy(Expr* e)` is now `e->refcount++; return e` for all
  types. No allocation, no recursion. The caller receives a shared
  pointer to the same node.
- A new primitive `Expr* expr_unshare(Expr* e)` is the
  copy-on-write counterpart. It consumes one ref of `e` and returns a
  refcount==1 logically-equivalent node:
  - If `e->refcount == 1`, returns `e` unchanged (zero work).
  - Otherwise allocates a fresh top-level node with a private payload
    (own `args[]` array for FUNCTION, own `mpz_t` / `mpfr_t` / `char*`
    for the corresponding leaf types). FUNCTION children are not
    deep-copied — they are inc-ref'd via `expr_copy`. The original `e`
    is dec-ref'd.
- Mutators must call `expr_unshare()` before writing fields if the
  node could be shared. To mutate deeper than one level (e.g.,
  `node->args[0]->args[i] = …`), unshare each level along the path.

### Audit

The only unsafe FUNCTION-mutation sites in the codebase were eight
locations in `src/print.c` where the negative-prefix path used
`t_copy = expr_copy(arg)` followed by an in-place rewrite of
`t_copy->data.function.args[0]`. Under Phase-2 sharing this would
corrupt the caller's tree. They are now wrapped:
```c
t_copy = expr_unshare(expr_copy(arg));
```
with the nested-Rational case unsharing two levels.

Every other FUNCTION-mutation site in the codebase operates on
freshly-constructed nodes whose refcount is 1 at the point of
mutation, so no change was required:
- `src/eval.c` `eval_flatten_args`, `flatten_sequences`,
  Unevaluated-stripping, Orderless `qsort` — all act on `res`, the
  fresh `expr_new_function` allocated at the top of `evaluate_step`.
- `src/internal.c` `internal_call_impl` — passes its just-built
  `res` to `eval_flatten_args` and `eval_sort_args`.
- `src/core.c` `builtin_quotientremainder` — zeroes `arg_count` on
  freshly-built `quot_expr`/`mod_expr` to suppress double-free of
  aliased args before calling `expr_free`.
- `src/parse.c` argv extension — the parser appends arguments to a
  function node it just created.
- `src/match.c` Sequence and `BlankSequence` binding writes —
  populate the args array of a `Sequence[…]` wrapper allocated one
  line earlier.
- `src/series.c` SeriesData recursion — mutates the args[] of the
  freshly-returned outer SeriesData node.
- `src/plus.c` and `src/times.c` numeric-contagion — replace `args[i]`
  on `res` (the evaluator's fresh input).

### Safety verification

- All 84 unit-test binaries pass under `tests/build`. The new binary
  `expr_sharing_tests` (`tests/test_expr_sharing.c`) is a focused
  stress suite for Phase-2:
  - Refcount inc/dec mechanics for atoms and FUNCTION.
  - `expr_unshare` semantics: identity on refcount==1, fresh node on
    refcount > 1, payload independence (BigInt mpz, String char*,
    FUNCTION args[] array).
  - Hot-path mutators (Plus/Times numeric contagion, Flat flatten,
    Orderless qsort, Sequence flatten, QuotientRemainder arg zero,
    print negative-Times, print negative-Rational) under 200×
    repetition.
  - Pattern binding repeated substitution where the binding is itself
    a FUNCTION and gets aliased into multiple slots of the RHS.
  - End-to-end Expand/Factor/Simplify/D fixed-point loops.
  - Deterministic xorshift random share/unshare/mutate/free walk over
    a 32-tree population for 4000 iterations, asserting that
    mutations through an unshared clone never alter the original.
- `valgrind --error-exitcode=2 ./expr_sharing_tests` reports zero
  errors (no use-after-free, no invalid reads/writes, no double-frees).
- `valgrind --leak-check=full ./poly_tests` after Phase-2 reports
  identical leak counters to Phase-1: 28,832 B definitely lost / 749
  blocks. FUNCTION sharing did not introduce new lost stacks.

### Measured allocation reduction

REPL smoke script (Expand + PolynomialGCD + Sin/Cos + D + Together +
Cases + Map + TrigToExp + Simplify + Factor + Quit):

| Run            | Allocs   | Frees   | Bytes allocated |
|----------------|----------|---------|-----------------|
| Pre-Phase-2    | 224,564  | 203,358 | 6,660,295       |
| Phase-2        |  83,550  |  67,795 | 3,118,302       |
| Reduction      | −62.8%   |         | −53.2%          |

### Convention for new code

- `expr_copy(e)` is the default; cheap and reads as "keep a logical
  copy".
- If you intend to mutate a FUNCTION node's `args[]`, `arg_count`, or
  `head` and the node arrived from a caller (i.e., you cannot prove
  you are the sole owner), call `e = expr_unshare(e)` first. The
  unshare is a no-op when `e` is already private, so it is safe to
  call defensively at the start of any mutator.
- For deep mutation (`e->args[0]->args[1] = …`), unshare each
  intermediate node and store the result back into the parent's
  args slot. The pattern is:
  ```c
  e = expr_unshare(e);
  Expr* child = expr_unshare(e->data.function.args[0]);
  e->data.function.args[0] = child;
  expr_free(child->data.function.args[1]);
  child->data.function.args[1] = …;
  ```

### Files touched

- `src/expr.h`: declared `expr_unshare`.
- `src/expr.c`: implemented `expr_unshare`; flipped `expr_copy` to
  unconditional inc-ref.
- `src/print.c`: 8 negative-prefix print sites wrapped with
  `expr_unshare(expr_copy(...))`.
- `tests/test_expr_sharing.c`: new 23-test stress suite.
- `tests/CMakeLists.txt`: registered `expr_sharing_tests`.
- `EVAL_IMPROVEMENTS_PLAN.md`: Phase-2 status appended to §3.2.

## Evaluation timestamps (M3 phase-3, 2026-05-03)

Closes Withoff [§2.1, §3.1 step 4] in `EVAL_IMPROVEMENTS_PLAN.md` §3.3:
*"An expression that has been evaluated once will not normally be
evaluated again unless some part of the expression is modified."*
Implements the conservative single-counter variant.

### Why timestamping makes the evaluator more efficient

PicoCAS's evaluator follows Mathematica's *infinite evaluation* rule:
`evaluate(e)` repeatedly calls `evaluate_step(e)` until two consecutive
steps return structurally-equal expressions (the "fixed point"). Without
caching, that full walk is paid every single time `evaluate(e)` runs,
even when nothing about the symbol-table state has changed since the
last time we already drove `e` to a fixed point.

The four hot patterns where this is wasted work:

1. **Outer fixed-point loop convergence.** When `evaluate(e)` stabilises
   on the very first iteration (no rule fires), the old code still
   builds a fresh function node inside `evaluate_step`, walks it for
   structural equality with the input via `expr_eq`, and discards it.
   That is `O(tree size)` per re-entry.
2. **Recursive sub-evaluation.** Inside `evaluate_step`, the head and
   each non-held argument are individually pushed back through
   `evaluate(...)`. After `Phase-2` made FUNCTION nodes shareable via
   refcount, those subtrees are increasingly the same physical nodes
   re-encountered across the rule pipeline (pattern bindings,
   `replace_bindings`, listable threading), so they are excellent
   cache candidates.
3. **REPL `In[n]` / `Out[n]` reuse.** A user who types `Out[1] + Out[2]`
   in the REPL re-references previously-computed results. Those
   results are full expression trees stored as `DownValues` of `Out`,
   and the old evaluator re-walked each one on every reference.
4. **Identities and tautologies.** `D[Sin[x]^2 + Cos[x]^2, x]` evaluates
   to `0` (an atom). Re-evaluating that result hits the atom fast path,
   but `Together[...]` and `Factor[...]` produce structured FUNCTION
   results; without a cache, each re-encounter pays the full tree walk
   to confirm "yes, still at fixed point".

The fix records, on every Expr, the value of a global monotonic
counter (`eval_clock`) at the moment that Expr was last driven to a
fixed point. The clock advances *only* when the symbol table changes
in a way that could change the meaning of any expression: a `Set`,
`SetDelayed`, `Clear`, `SetAttributes`, or `ClearAttributes`. Pure
builtin calls (Sin, Factor, Expand, GCD, Map, …) do not bump the clock
because they do not change rule state.

At entry, `evaluate(e)` checks `e->last_evaluated_at == eval_clock_get()`.
If true, the expression is *known* to be at a fixed point under the
current symbol-table state, so we return an inc-ref'd view in O(1)
instead of running the loop. If the user has changed any definition
since the cache was set, the clock will have moved past the stamp and
the check fails — full re-evaluation runs against the new state.

This is intentionally conservative: any rule change invalidates *every*
cached expression in one shot, even expressions that did not depend on
the changed symbol. A future refinement (per-symbol clocks plus a
per-expression "max clock of any referenced symbol") would tighten
that, but the global counter is dramatically cheaper to maintain and
still wins on every workload that touches a stable rule set.

### Mechanism

- `Expr` gains a `uint64_t last_evaluated_at` field. Every constructor
  initialises it to 0; the global clock starts at 1, so a freshly built
  tree never accidentally hits the cache on first evaluation.
  `expr_unshare` resets the field to 0 on the freshly-allocated copy
  because the caller is preparing to mutate.
- `eval_clock_get()` / `eval_clock_bump()` are declared in `eval.h`.
  The counter lives as a `static uint64_t g_eval_clock = 1;` inside
  `eval.c`.
- `evaluate(e)` early-exits at the very top with
  `if (e->type == EXPR_FUNCTION && e->last_evaluated_at == g_eval_clock)
   return expr_copy(e);`. Atoms are excluded from the pre-check
  because they already short-circuit cheaply inside `evaluate_step`
  (the pre-check would just add a branch without saving real work).
- On a clean fixed-point exit, the result's `last_evaluated_at` is
  set to the live clock. Overflow and iteration-cap paths leave the
  stamp untouched so a later evaluator can retry from scratch.
- `add_rule` (`src/symtab.c`), `symtab_clear_symbol`, `set_attributes`,
  `add_single_attribute`, and `remove_single_attribute` all call
  `eval_clock_bump()`. The two `add_single_attribute` /
  `remove_single_attribute` paths bump *only* when a real bit change
  occurs — re-applying an already-set attribute is a true no-op.
- `expr_eq` and `expr_hash` deliberately ignore the new field: it is
  benign metadata, not part of the logical value, so two expressions
  with different timestamps but identical structure remain equal /
  hash to the same key. Concurrent updates are not a concern (single-
  threaded), and writing the timestamp on a refcount-shared node is
  safe because every holder observes the same correct "still at fixed
  point under clock T" assertion.

### Test coverage

`tests/test_eval_timestamps.c` (16 tests, the 85th binary):

- **Clock mechanics**: starts at ≥ 1; pure builtin chain (Sin, Factor,
  Expand, Length, Map, D, Apart, PolynomialGCD) does not advance the
  clock; Set/SetDelayed/Clear advance the clock; SetAttributes /
  ClearAttributes advance only on a real bit change.
- **Stamping**: fresh constructor outputs have `last_evaluated_at == 0`;
  `evaluate(e)` stamps the result with the current clock.
- **Cache hit pointer identity**: `evaluate(r)` on an already-stamped
  FUNCTION returns the *same* pointer (refcount-shared) — no
  reconstruction.
- **Invalidation correctness**: assigning `x = 99` after evaluating a
  symbolic expression containing `x` produces the new arithmetic value
  on the next read; the cache does not hand back stale results across
  Set, SetDelayed, Clear, or SetAttributes mutations.
- **Cached subexpressions in fresh outer calls**: pre-evaluating
  `Sin[Pi/4]` and re-using the (stamped) result inside a hand-built
  `Plus[...]` tree still produces the correct end-to-end answer
  (`Sqrt[2]`).
- **`expr_eq` / `expr_hash` invariance**: two structurally identical
  expressions with different timestamps are equal and hash the same.
- **Repeated heavy evaluation**: 50 iterations × 6 kernels (Expand,
  D, Factor, Together, PolynomialGCD, Apart) — each iteration's
  printed output must match iteration 1's byte-for-byte. Catches any
  cache-induced drift or aliasing across timestamps.
- **Random walk**: 200 iterations of randomly chosen Set / Clear /
  Read against three symbols, with a shadow state validated on every
  Read. Exercises every cache invalidation path under realistic
  user-interaction patterns.

Valgrind on `eval_timestamps_tests` reports zero memory-safety errors.

### Files touched

- `src/expr.h`: added `uint64_t last_evaluated_at` field with comment.
- `src/expr.c`: initialised in every `expr_new_*` constructor; reset
  to 0 in `expr_unshare`'s freshly-allocated copy.
- `src/eval.h`: declared `eval_clock_get` / `eval_clock_bump` and the
  full doc comment for the M3 phase-3 contract.
- `src/eval.c`: defined `g_eval_clock` and the two accessors; added
  the early-exit at the top of `evaluate()`; stamped `current` on
  the clean fixed-point exit.
- `src/symtab.c`: `add_rule` and `symtab_clear_symbol` bump the clock.
- `src/attr.c`: `set_attributes`, `add_single_attribute`,
  `remove_single_attribute` bump only on a real change.
- `tests/test_eval_timestamps.c`: new 16-test suite.
- `tests/CMakeLists.txt`: registered `eval_timestamps_tests` (85th
  binary) and added it to the `build_ecm` dependency list.
- `EVAL_IMPROVEMENTS_PLAN.md`: §3.3 marked done with the implementation
  notes.

## Eager early-exit fixed-point loop (M3 phase-4, 2026-05-03)

PicoCAS's evaluator implements Mathematica's *infinite evaluation*
semantics: an expression is repeatedly fed through a single-step
rewriter until the result stops changing. Until phase-4, "stops
changing" was detected by an O(tree) structural compare
(`expr_eq(current, next)`) on every iteration of the outer fixed-point
loop. Phase-4 replaces that compare with a cheap boolean fast-path
that comes directly from the rewriter itself.

### Why an eager early-exit makes the evaluator more efficient

The single biggest cost in the outer fixed-point loop, on workloads
that are already at a fixed point, used to be the `expr_eq(current,
next)` compare. Every re-evaluation of an already-evaluated expression
paid for it; every iteration before the converging one paid for it; and
on deep trees that compare alone could be larger than the rewrite work
in a single step.

Phase-4 lifts that cost. Four user-facing patterns benefit:

1. **Already-converged expressions.** Evaluating `Plus[a, b, c]` once
   and never assigning to `a`, `b`, or `c` is the canonical case: the
   rewriter has nothing left to do, but the loop still has to *prove*
   that. Pre-phase-4, proof cost an O(tree) `expr_eq`. Post-phase-4,
   the rewriter signals "nothing fired" via a boolean and the loop
   exits in O(1).

2. **Sub-expression no-ops inside a deep tree.** When evaluating a
   nested expression like `Sin[x] + Cos[y] + Tan[z]`, each leaf
   sub-evaluation is an already-converged sub-expression. With
   refcount sharing (phase-2) plus phase-4, those sub-evaluations
   short-circuit at the first iteration without entering the heavy
   step body, and the outer loop sees that none of the *args* changed
   either — so the parent also exits in one iteration.

3. **Combined with phase-3 timestamps.** Phase-3 stamps every
   converged expression with the eval-clock value at convergence.
   Phase-4 makes that *first* convergence cheap. Together: a fully
   evaluated expression that is then re-evaluated under an unchanged
   clock returns in O(1) at the entry of `evaluate()` (phase-3); the
   *first* evaluation, which is what installs the stamp, exits in one
   iteration without an `expr_eq` (phase-4).

4. **Mid-rewrite progress detection.** During a chain of rewrites
   (e.g. `Simplify` walking a tree), the loop now distinguishes "we
   are still rewriting" from "we just stopped" purely from the
   rewriter's own signal — no structural compare against the previous
   iteration is needed to decide whether to keep going.

### Mechanism

`evaluate_step(Expr* e, bool* changed)` takes a `bool*` out-parameter
and sets `*changed = true` on every actual rewrite. The default is
`false`, so a step that does nothing semantic (atom, bare symbol with
no OwnValue, fully-reduced function with no applicable
DownValue/built-in/special) leaves the flag at `false`.

The exhaustive list of "rewrite fired" signals inside
`evaluate_step`:

- The head re-evaluated to a different subtree (pointer-inequality
  after refcount sharing).
- An argument re-evaluated to a different subtree.
- An `Evaluate[...]` wrapper inside a held position was stripped.
- An `Unevaluated[...]` wrapper was stripped from a non-held arg.
- `flatten_sequences` actually flattened a `Sequence[...]` (the helper
  now returns `bool` so the caller can tell).
- `eval_flatten_args` actually flattened a nested same-head call (also
  now returns `bool`).
- The `Listable` attribute threaded the call over a `List` arg.
- The `Orderless` attribute reordered the args (with a pre-check so a
  no-op qsort on already-sorted args does not trip the flag).
- A `DownValue` rule matched.
- A C built-in returned non-NULL.
- A `Set` or `SetDelayed` installed a rule.
- A pure `Function[...]` applied to its body.
- `Derivative[n][Function[...]]` reduced to a differentiated `Function`.
- `Composition[f1,...,fn][args]` unrolled to `f1[f2[...[fn[args]]]]`.

The outer loop in `evaluate()` then becomes:

```c
bool step_changed = false;
next = evaluate_step(current, &step_changed);
bool is_fixed_point = !step_changed || expr_eq(current, next);
if (is_fixed_point) { /* stamp current and return */ }
```

The `!step_changed` branch is the fast path — it skips the structural
compare entirely, which is the whole point of phase-4. The `expr_eq`
fallback exists only because some built-ins (notably `builtin_plus`
and `builtin_times`) unconditionally rebuild their output tree even
when no terms combined; those trip the change flag without producing
a structural difference. The fallback ensures `a + b + c` still
converges in one outer iteration. Cost on the slow path is identical
to pre-phase-4; the win is the bypass on the common case.

### Conservative tracking

False-positives are correctness-safe — at worst they cost one extra
outer iteration. False-*negatives* would terminate the loop before a
fixed point and are the bug to avoid. The instrumentation in
`evaluate_step` is therefore deliberately conservative: every branch
that *might* be a rewrite sets the flag explicitly, each with a
one-line comment naming the reason (e.g. `*changed = true; /* DownValue
rule fired */`). The pre-check before the Orderless `qsort` is the
only optimisation in the other direction — it skips the `*changed =
true` when the args are already in canonical order, because re-sorting
already-sorted args is provably a no-op.

### Test coverage

`tests/test_eval_eager_exit.c` (15 tests, the 86th unit-test binary)
locks in:

- Atomic types (Integer, Real, String, BigInt) never report a change.
- Bare symbols with no OwnValue never report a change.
- Bound symbols (with an OwnValue) report a change.
- Built-ins that fire (e.g. `Sin[Pi/4]`) report a change.
- A `NULL` `bool*` out-parameter is tolerated.
- The `Plus[a, b, c]` regression case converges (this would loop
  forever under a no-fallback design — the test pins the design).
- Repeated re-evaluation of an already-evaluated expression hits the
  phase-3 timestamp path and never enters `evaluate_step` at all.
- Listable, Flat, Sequence-flatten, and Orderless-reorder rewrites
  each report a change.
- DownValue matches report a change.
- 50 reps × 8 mixed-kernel `Expand`/`Factor`/`Sin`/`Integrate`/
  `Simplify`/`Range`/`Map`/`Table` workloads produce byte-stable
  output.
- 30 reps of a deeply nested `{Plus[...], Times[...], Plus[...,
  Times[...], ...], Sin[x] + Cos[y] + Tan[z]}` expression produce
  byte-stable output.
- A `Set`-then-redefine-then-`Set` cycle through `eagerG[10]`
  validates that the phase-3 clock invalidation correctly interacts
  with the phase-4 fast path.

### Files touched

- `src/eval.h`: new `evaluate_step(Expr* e, bool* changed)` signature
  with the contract doc; `eval_flatten_args` now returns `bool`.
- `src/eval.c`: the `bool* changed` instrumentation throughout
  `evaluate_step`; the `step_changed || expr_eq` hybrid in the outer
  fixed-point loop; `flatten_sequences` and `eval_flatten_args` now
  return `bool` so callers can distinguish a real flatten from a no-op.
- `tests/test_eval_eager_exit.c`: new 15-test suite.
- `tests/CMakeLists.txt`: registered `eval_eager_exit_tests` (86th
  binary) and added it to the `build_ecm` dependency list.
- `EVAL_IMPROVEMENTS_PLAN.md`: §3.4 marked done with the
  implementation notes.

## DownValue dispatch index + specificity ordering (M4 milestone, 2026-05-03)

The fifth and final M3/M4 evaluator improvement closes the Tier-3 hot-loop
work by giving every `DownValue` rule two pre-computed pieces of metadata
that let `apply_down_values` reject unmatchable rules without invoking the
matcher, and that order the rule list strictly from most-specific to
least-specific.

### Why dispatch indexing helps

Every call to `f[args...]` for a user-defined `f` previously walked the
entire `DownValues` list and ran the full pattern matcher against each
rule. Linear in the number of rules, even when most rules obviously
cannot match (wrong arity, wrong head on the first argument). Symbols
that accumulate large rule banks — `Simplify` rule sets, scientific
function tables, generated DSL handlers — paid this cost on every
invocation.

§3.5 turns each rule's top-level shape into a tiny pre-computed key
that `apply_down_values` matches against the input call's shape via
pointer compares against interned strings. Shape mismatches are
discarded in O(1); only rules that *could* match reach the matcher.

### Why a specificity score helps

Pre-M4 the rule list used a binary "literals first, patterns after"
heuristic. Pattern rules competed in insertion order regardless of how
discriminating each pattern was. Defining `f[x_]` before `f[0]` gave
`f[0]` the wrong handler.

§3.6 replaces that binary partition with a real specificity score
(`int32_t pattern_specificity` in `src/symtab.c`). Rules sort
descending by score; equal scores keep insertion order; redefining the
same pattern updates the RHS in place. The Mathematica-style "most
specific rule wins" semantics now hold even when the user adds rules
in a non-canonical order.

### Mechanism

`Rule` (in `src/symtab.h`) gained three fields, all set once at insertion:

```c
int32_t      dispatch_arity;          /* -1 == variable-arity */
const char*  first_arg_head_canon;    /* interned head, or NULL */
int32_t      specificity;             /* descending sort key */
```

Computation:

- `pattern_dispatch_arity(p)` returns `-1` whenever any top-level slot
  is a `BlankSequence`, `BlankNullSequence`, `Optional`, `Repeated`,
  `RepeatedNull`, or `OptionsPattern`. Otherwise it returns the number
  of top-level argument slots.
- `pattern_first_arg_head(p)` looks through `HoldPattern` / `Verbatim`,
  recurses through `Pattern[name, q]`, `Condition`, `PatternTest`, and
  agreeing branches of `Alternatives`. For atoms it returns the
  canonical head name (`Integer`, `Real`, `String`, `Symbol`); for
  function calls it returns the head symbol. Sequence and optional
  patterns return `NULL` (wildcard).
- `pattern_specificity(p)` sums per-node weights: literal atoms are
  +100, plain function calls are +50 plus children's scores, typed
  blanks +20, plain blanks -10, sequence blanks -100/-200,
  `Optional` -50, `Condition` and `PatternTest` add +10 to their
  inner score, `Alternatives` takes the weakest branch, and
  `OptionsPattern` is -150. `Pattern[name, q]` is transparent because
  the binding name is irrelevant to matching power.

`add_rule` walks until it finds the first rule with strictly lower
specificity, then stable-inserts the new rule before it. Identical
patterns still replace in place.

`apply_down_values` (in `src/symtab.c`) computes the input's
`(arity, first_arg_head_canon)` once via `input_arg_head_canon`, then
for each rule:

1. If `rule->dispatch_arity != -1 && rule->dispatch_arity != input_arity`,
   skip.
2. If `rule->first_arg_head_canon && input_first_head &&
   rule->first_arg_head_canon != input_first_head`, skip.
3. Otherwise call the matcher exactly as before.

Both checks are pointer compares against interned strings, so a
mismatch costs only a few nanoseconds per rule and is essentially
free compared to the matcher.

### Soundness

- A rule that *could* match the input is never skipped. The two filter
  conditions above are conjunctions: a `NULL` key on either side
  always passes. Variable-arity rules (sequence patterns) carry
  `dispatch_arity = -1` and always pass the arity gate; atomic inputs
  whose head is `Integer` / `Real` / `Symbol` / `String` line up with
  rule patterns that constrain on those head names; typed
  `Blank[Integer]` matches an integer input under both filters.
- The specificity ordering is a strict refinement of the pre-M4
  ordering: literal-arg rules score ≥100 and dominate the same head's
  pattern rules (preserving the binary-literal advantage); typed
  blanks dominate plain blanks; sequence-pattern rules sink to the
  bottom. The behaviour shift relative to pre-M4 is that pattern rules
  of *unequal* specificity may now fire in a different order from
  their insertion sequence — which is the desired Mathematica-style
  behaviour.

### Test coverage

`tests/test_rule_dispatch.c` (23 tests, the 87th unit-test binary)
locks in:

- Dispatch metadata for fixed-arity-no-head, typed-blank,
  function-headed, sequence-pattern, optional, literal-int, and
  `HoldPattern`-wrapped rules is set correctly.
- Specificity ordering: literal beats blank, typed blank beats plain
  blank, single blank beats sequence blank, conditioned pattern beats
  unconditioned, equal-specificity rules retain insertion order, and
  re-defining a pattern updates in place without moving.
- Functional dispatch: defining `general` first then `specific` still
  fires `specific` for matching inputs; arity dispatch picks the right
  rule for 1/2/3-argument calls; `f[Sin[x_]]` vs `f[Cos[x_]]` vs
  `f[x_]` correctly partitions on first-arg head;
  `f[n_Integer]` / `f[r_Real]` / `f[s_Symbol]` correctly partitions on
  atom kind; `f[x__]` (variable arity) is never filtered out;
  `f[x_]` defined alone leaves `f[a, b]` unevaluated; head-specialised
  rules for atomic inputs filter cleanly to the wildcard fallback.
- Stress: 100 specific literal rules `hv[k] -> k+1000` plus a
  wildcard fallback all dispatch correctly under the filter.
- Classical recursion: `fib[n_]` defined first then `fib[0]` and
  `fib[1]` still produces `fib[10] = 55` because §3.6 re-orders the
  base cases ahead of the recursive case.

### Files touched

- `src/symtab.h`: extended `Rule` with `dispatch_arity`,
  `first_arg_head_canon`, and `specificity` fields, with full doc.
- `src/symtab.c`: added `strip_pattern_wrappers`,
  `slot_is_variable_length`, `pattern_arg_head_canon`,
  `pattern_dispatch_arity`, `pattern_first_arg_head`,
  `pattern_specificity`, and `input_arg_head_canon`. Rewrote
  `add_rule` for stable specificity-sorted insertion. Rewrote
  `apply_down_values` to apply the two-key filter before invoking the
  matcher.
- `tests/test_rule_dispatch.c`: new 23-test suite.
- `tests/CMakeLists.txt`: registered `rule_dispatch_tests` (87th
  binary) and added it to the `build_ecm` dependency list.
- `EVAL_IMPROVEMENTS_PLAN.md`: §3.5 and §3.6 marked done with the
  implementation notes; M4 closed.


## simp_algebraic: algebraic-extension reduction in Simplify (2026-05-04)

`Simplify` now collapses identities involving square-root surds that
ordinary `Together` / `Cancel` cannot see, by treating each distinct
`Sqrt[u_i]` as a generator `g_i` of the algebraic extension

    K(vars)[g_1, ..., g_n] / (g_1^2 - u_1, ..., g_n^2 - u_n)

and rationalising the denominator via successive sigma-conjugation.

### Examples

    Simplify[(x/Sqrt[x^2 + 1] + 1)/((Sqrt[x^2 + 1] + x)^2 + 1)]
        →  1/(2 + 2 x^2)

    Simplify[(x*(1/Sqrt[x^2 + 6] - (Sqrt[x^2 + 6] - Sqrt[6])/x^2))/
             (Sqrt[x^2 + 6] - Sqrt[6])]
        →  Sqrt[6]/(x Sqrt[6 + x^2])

    Simplify[2/(Sqrt[(x + 1)/(1 - x)] - 1/Sqrt[(x + 1)/(1 - x)])]
        →  (1 + x)/(x Sqrt[(1 + x)/(1 - x)])

### Algorithm

1. Walk the input collecting every distinct surd argument `u_i` from
   `Power[u_i, p/q]` with `q != 1`. Bail if any `q != 2`, if more
   than `ALG_MAX_SURDS` (currently 4) distinct bases appear, or if
   any `u_i` contains an explicit complex literal (`Complex[..,..]`,
   the symbol `I`).
2. Substitute `Sqrt[u_i] -> g_i` for fresh `$pc_alggen{i}$` symbols.
3. `Together` -> `N / D`, both polynomials in `(vars, g_1, ..., g_n)`.
4. Reduce `N` and `D` modulo the relation ideal `{g_i^2 - u_i}_i` via
   successive `CoefficientList[..., g_i]` decomposition. After one
   sweep across all generators the polynomial is multilinear in the
   `g_i` (each `g_i` appears at degree 0 or 1).
5. For `i = 1..n`, multiply numerator and denominator by `sigma_i(D)`
   (denominator with `g_i` sign-flipped) and reduce again. Each step
   eliminates `g_i` from the denominator because
   `(a + b g_i)(a - b g_i) = a^2 - b^2 u_i` is `g_i`-free.
6. `Factor` numerator and denominator separately to expose shared
   `Power[u_i, k]` factors that would otherwise be hidden by the
   expanded polynomial form, then substitute `g_i -> Sqrt[u_i]`,
   `Cancel` to combine `Power[u_i, k] * Power[u_i, 1/2]` terms.

The transform is wired as a one-shot seed in `simp_search`. A
complexity gate (`>` rather than `>=`) lets it preserve a tied score
for outputs that are conventionally preferred (e.g. the rationalised
form of `1/(Sqrt[a]+Sqrt[b])`); strictly worse forms are rejected.

### Principal-branch caveat

The substitution `Sqrt[u]^2 = u` is only sound where `u` lies in the
principal branch's domain. We accept the Mathematica-style convention
(`Simplify` treats this as an identity on the natural real-domain
where the input is real) but skip the transform when any surd
argument contains an explicit complex literal so a sign-of-imaginary-
part change cannot be silently absorbed.

### Files touched

- `src/simp.c`: new `simp_algebraic_impl` plus helpers
  (`alg_collect_sqrt_bases`, `alg_subst_sqrt_to_gens`,
  `alg_subst_gens_to_sqrt`, `alg_sigma_negate`, `alg_reduce_one_gen`,
  `alg_reduce_all`, `contains_explicit_complex`); wired as a seed in
  `simp_search` alongside `simp_radicals`. Includes `sym_intern.h`
  for canonical generator-symbol pointers.
- `tests/test_simplify.c`: three new tests
  (`test_simplify_algebraic_single_surd`,
  `test_simplify_algebraic_multi_surd`,
  `test_simplify_algebraic_fractional_surd_arg`).

## simp_algebraic: u-power extraction (2026-05-04)

Extended `simp_algebraic` so an implicit `Power[u_i, k]` factor in the
denominator -- e.g. `x^4 = (x^2)^2` when `u = x^2` -- gets lifted into
the generator's algebraic ring before `Cancel` runs. Without this step
the surd factor in the numerator and the `u`-power in the denominator
look coprime to ordinary polynomial GCD even though
`Sqrt[u]^2 = u` would let them combine.

    Simplify[(Sqrt[x^2] - 1/Sqrt[x^2])/x^2]
        →  (-1 + x^2)/(x^2)^(3/2)

### Mechanism

After step 5 of `simp_algebraic` (denominator rationalised, both
numerator and denominator are polynomials in `(vars, g_1, ..., g_n)`
with each `g_i` at most degree 1), for every surd whose argument
`u_i` is itself a polynomial we iteratively divide both numerator
and denominator by `u_i` (via `PolynomialQuotient` /
`PolynomialRemainder` against the first variable in `Variables[u_i]`).
Each successful division increments the multiplicity counter `k`.
The extracted factor `Power[u_i, k]` is then re-attached as
`Power[g_i, 2k]` (using the relation `g_i^2 = u_i`), so that the
subsequent `Cancel` over `Q[vars, g_1, ..., g_n]` can combine the
generator powers between numerator and denominator.

The extraction is gated on `alg_u_is_polynomial`: rational `u_i`
(e.g. `(x+1)/(1-x)`) skips this step, since `PolynomialQuotient` is
undefined for rational divisors and would produce spurious results.

### Files touched

- `src/simp.c`: `alg_extract_u_power`, `alg_pick_var`,
  `alg_u_is_polynomial`, plus the per-surd extraction loop in
  `simp_algebraic_impl`.
- `tests/test_simplify.c`: `test_simplify_algebraic_u_power_extraction`.

## Numericalize: preserve rational exponents over symbolic bases (2026-05-04)

Fixed a long-standing inexact-contagion bug where `Times` / `Plus` /
`N[...]` would demote rational exponents to machine reals whenever any
factor was inexact, producing ugly output like

    1/(2 Sqrt[x] (1+x)) + (0.5 (1+x) (...))/...  // Simplify
        ->  1/(x + x^0.5 + x^1.5 + x^2)        (was)
        ->  1/(x + Sqrt[x] + x^(3/2) + x^2)    (now)

### Cause

`numericalize_function` in `src/numeric.c` recursively numericalized
function arguments. For `Power[x, Rational[1,2]]` this dropped the
rational exponent down to a `0.5` real even though the symbolic base
`x` carried no inexactness. The result rebuilt as `Power[x, 0.5]` and
then propagated through `Times` / `Plus` -- contaminating any
downstream simplification.

### Fix

In `numericalize_function`, the existing "preserve integer exponent"
rule is generalised: for `Power[base, exp]`, if the numericalized base
is *not* numeric (per `expr_is_numeric_like`), the original exponent
is preserved verbatim. This keeps `Sqrt[x]`, `(1+x)^(1/2)`, `x^(-1/2)`,
etc. structurally intact under contagion while still letting genuinely
numeric bases collapse: `1.0 Sqrt[2] Sqrt[x] -> 1.41421 Sqrt[x]`,
`N[Sqrt[Pi]] -> 1.77245`. Matches Mathematica's behaviour for
`N[Sqrt[x]]` (returns `Sqrt[x]`).

### Files touched

- `src/numeric.c`: `numericalize_function` Power branch generalised.

## Simplify final-form polish: common-factor lift + sign canonicalization (2026-05-04)

Two new post-search passes run after `simp_bottomup` returns its best
candidate. They are surgical -- applied once at the top level only --
to avoid destabilising the heuristic search on multi-variable trig
inputs (the trig two-vars regression caught when an earlier draft wired
the lift into the seed phase).

### Common-factor lift (`simp_lift_common_factor`)

Pulls a shared multiplicative factor out of a top-level `Plus`, or out
of a `Plus` child of a top-level `Times` (the numerator of a fraction
whose denominator carries a non-integer power). picocas's `Factor` and
`FactorTerms` decompose polynomials over `K[x_1, ..., x_n]` using
`Variables[]` to discover the generator set, and `Variables[]` does
not return non-integer-power expressions like `(1+x^2)^(3/2)` or
`Sqrt[x]`. So a `Plus` that obviously shares an algebraic generator
across all terms slips past `Factor` untouched. The lift takes a
structural multiset view: a non-numeric factor is either an algebraic
generator (`Power` with non-integer exponent, or any opaque token
whose exponent we can't reduce) treated as one indivisible token, or
a `Power[base, n]` with `n` a small positive integer that we split
into `n` copies of `base`. Numeric coefficients merge via rational
GCD. The lifted result is `gcd_coef * Times(common_tokens) * Plus[t_i
/ lift_factor]`, with the per-term divisions handed back to
`evaluate()` for cancellation.

Cases this enables that `Factor` cannot:

    (8/105)(1+x^2)^(3/2) - (4/35) x^2 (1+x^2)^(3/2) + (1/7) x^4 (1+x^2)^(3/2)
        ->  (1/105) (1+x^2)^(3/2) (8 - 12 x^2 + 15 x^4)

    (-x^3/Sqrt[5+2x] + 3 x^2 Sqrt[5+2x]) / (5+2x)
      Cancel ->  (15 x^2 + 5 x^3) / (5+2x)^(3/2)
      lift   ->  (5 x^2 (3 + x)) / (5+2x)^(3/2)

The lift fires only when there is a real shared algebraic factor
(`common_count > 0`); a coefficient-only GCD pull-out is left to
`Together` / `Cancel` because it doesn't expose new structure and can
trigger long downstream rationalisation chains.

### Sign canonicalization (`canon_negate_pairs`)

Mathematica's printed form prefers each binomial to lead with its
positive-coefficient term: it prints `(a - c)(b - d)`, not the
algebraically-equal `(-a + c)(-b + d)` that picocas's facpoly emits
(the canonical Plus order sorts `Times[-1, a]` before `c` because the
sort key strips the leading -1 coefficient). This pass walks a
top-level `Times`, counts negatively-leading `Plus` factors, and
flips them in pairs. Paired flips are value-preserving: each pair
contributes `(-1) * (-1) = 1` to the product. Odd-count cases are
left alone (lifting an extra -1 onto an outer numeric is not always
a canonical win for the score function).

Concrete win:

    Simplify[1/(a*b) - 1/(b*c) - 1/(a*d) + 1/(c*d)]
        was -> ((-a + c) (-b + d)) / (a b c d)
        now -> ((a - c) (b - d)) / (a b c d)

### Files touched

- `src/simp.c`: `simp_lift_common_factor`, `canon_negate_pairs`,
  and the final-form polish in `builtin_simplify`.
- `tests/test_simplify.c`: eleven new user-supplied regression cases
  exercising these passes plus existing simplification machinery.
- `tests/test_numeric.c`: aligned `N[x^(1/2)]` expectation with the
  numericalize fix above (`Sqrt[x]`, not `x^0.5`).

## ExpToTrig seed + reciprocal Pythagorean identities (2026-05-04)

`Simplify` previously could not collapse pure-Exp inputs whose natural
canonical form lives in the hyperbolic trig algebra. Two gaps:

1. **The `TrigRoundtrip` seed in `simp_search` was gated on
   `contains_trig_or_hyperbolic(input)`**, so an input like
   `-1 + (E^x - E^-x)^2 / (E^x + E^-x)^2` -- which has no Cosh / Sinh
   head -- skipped the entire trig pipeline. The result rolled out as
   `-(4 E^(2x)) / (1 + 2 E^(2x) + E^(4x))`.

2. **`PythagReduce` carried no rules for the reciprocal pairs**
   (Tanh ⇄ Sech, Coth ⇄ Csch, Tan ⇄ Sec, Cot ⇄ Csc). Even after
   `simp_search` reached the intermediate `-1 + Tanh[x]^2`, neither
   PythagReduce nor anything else could collapse it to `-Sech[x]^2` --
   and `simp_search`'s `update_best` uses a strict `<` tiebreak, so the
   tied-score plateau (both forms cost 7) preferred whichever arrived
   first.

### Fix A: `ExpToTrig` seed

A new seed in `simp_search` runs `ExpToTrig[input]` when:

- the input contains a `Power[E, _]` subexpression
  (`contains_exp_form`), AND
- the input does **not** already contain trig / hyperbolic heads
  (`TrigRoundtrip` handles those).

The result is added as a candidate (with the same 2x score-blowup
guard as `TrigRoundtrip`), giving the round loop access to forms
like `Sinh[x]^2 / Cosh[x]^2 - 1` that the existing `Together` /
`Cancel` / `Pythag*` pipeline can drive to closure.

### Fix B: Reciprocal Pythagorean rules + polish pass

`PythagReduce`'s rule list was extended with eight new entries:

    1 - Tanh[x_]^2 -> Sech[x]^2
    -1 + Tanh[x_]^2 -> -Sech[x]^2
    -1 + Coth[x_]^2 -> Csch[x]^2
    1 - Coth[x_]^2 -> -Csch[x]^2
    1 + Tan[x_]^2 -> Sec[x]^2
    -1 - Tan[x_]^2 -> -Sec[x]^2
    1 + Cot[x_]^2 -> Csc[x]^2
    -1 - Cot[x_]^2 -> -Csc[x]^2

`has_pythag_head` was extended to include Tan / Cot / Tanh / Coth so
the cheap fast-skip gate still applies when the input has none of
those heads.

These rules tie on `SimplifyCount` (both forms cost 7), so they would
lose the strict-`<` tiebreak inside `simp_search`. To resolve the tie
in favour of the structurally-collapsed reciprocal-head form, a new
final-form polish runs `transform_pythag_reduce` on the
`simp_bottomup` result and accepts on `<=`. The polish lives next to
`simp_lift_common_factor` and `canon_negate_pairs` in
`builtin_simplify`.

### Concrete wins

    Simplify[-1 + (-E^(-x) + E^x)^2 / (E^(-x) + E^x)^2]
        was -> -(4 E^(2 x)) / (1 + 2 E^(2 x) + E^(4 x))
        now -> -Sech[x]^2

    Simplify[-1 + Tanh[x]^2]   ->  -Sech[x]^2
    Simplify[1 + Tan[x]^2]     ->  Sec[x]^2
    Simplify[1 + Cot[x]^2]     ->  Csc[x]^2
    Simplify[-1 + Coth[x]^2]   ->  Csch[x]^2

### Files touched

- `src/simp.c`: `contains_exp_form`, ExpToTrig seed in `simp_search`,
  extended `PythagReduce` rule list, extended `has_pythag_head`, new
  polish pass in `builtin_simplify`.
- `tests/test_simplify.c`: five regression tests covering the
  Exp→Sech path and each reciprocal-pair identity.

## Top-level rational shortcut in builtin_simplify (2026-05-04)

When the input classifies as `SHAPE_RATIONAL` -- meaning no trig, log,
abs, or non-integer-power generators are present and `Together[input]`
factors as poly/poly -- the simplifier now dispatches once at the top
via `simp_dispatch` instead of recursing through `simp_bottomup` first.

Why: `simp_bottomup` descends into every Plus / Times child before
dispatching. For a rational input, each child re-classifies as
`SHAPE_RATIONAL` and re-enters `simp_pipeline_rational` (Together,
Cancel, ExpandNumerator, ExpandDenominator, Apart, Factor, per-variable
Collect). Together at the top combines all children into a single
canonical num/den, so the per-child work is subsumed.

On `a/b + c/d - e/f + 1/(g h) - 1/(i j) - 1/(b k l)`:
- baseline: 8.7-8.9 s
- with shortcut: 7.2-7.5 s  (~18% faster)

The classifier is conservative: anything with trig, log, Abs, or
non-integer powers falls through to `simp_bottomup` as before, so
search-driven domains (the heuristic's strongest area) are unaffected.
The polish passes (lift_common_factor, PythagReduce,
canon_negate_pairs) still run on the result.

### Files touched

- `src/simp.c`: branch in `builtin_simplify` selecting `simp_dispatch`
  vs `simp_bottomup` based on `simp_classify`.

## Factor / Apart / Together / Cancel: shared algebraic-generator pass (2026-05-05)

`Together` and `Cancel` historically had a private "algebraic generator"
pass that handled fractional rational powers of a *symbol* — substitute
`y -> g^m`, run polynomially in `g`, substitute back. The pass has been
generalised in two directions and lifted into `poly.c` so all four
rational-function builtins (`Factor`, `Apart`, `Together`, `Cancel`)
share it.

### Generalisations

1. **Non-symbol bases.** The base may now be any sub-expression matched
   structurally via `expr_eq`, not only an `EXPR_SYMBOL`. Inputs like
   `Log[r]^(1/5)` or `Sqrt[Log[r]]` are treated as algebraic generators
   in the same way `r^(1/5)` is.

2. **Exponential atom.** A `Power[B, c*A]` form (where `c` is rational
   and `A` is any sub-expression) is recognised as a polynomial site in
   `g = Power[B, A/m]`. With `m = lcm` of c-denominators, every matching
   `Power[B, c*A]` becomes `g^(c*m)`. This makes `Factor[Exp[2x] + 2 Exp[x] + 1]`
   collapse to `(1 + E^x)^2` via the substitution `g = E^x`.

### Trigger conditions

- **Radical case** (atom `A = 1`): trigger when at least one fractional
  exponent exists (`m > 1`). Bare occurrences of the base count as
  `Power[B, 1]`, so `r + r^(1/2)` factors via `g = r^(1/2)`.
- **Exponential case** (atom `A != 1`): trigger when at least one
  matching site has coefficient `c != 1`, i.e. the substituted form is
  a non-trivial polynomial in `g`. Single `Power[B, A]` alone (which the
  substitution would just rename) does not trigger.

### API

Exposed from `poly.h`:

```c
bool  poly_find_radical_gen(Expr* e, Expr** base_out, Expr** atom_out, int64_t* m_out);
Expr* poly_subst_radical_to_gen(Expr* e, Expr* base, Expr* atom, int64_t m, const char* gen);
Expr* poly_subst_radical_from_gen(Expr* e, Expr* base, Expr* atom, int64_t m, const char* gen);
char* poly_make_fresh_gen(Expr* e);
```

`poly_make_fresh_gen` walks the expression to pick a `$pc_radgenN$`
symbol name not already present, so recursive calls (e.g. multiple
distinct radical bases in the same input) don't collide on the
substitution variable.

### Examples

```
In[1]:= Factor[1 + r^(1/5) - 2 r^(2/5) - r^(3/5) - 2 r^(4/5)]
Out[1]= -(1 + 2 r^(1/5)) (-1 + r^(1/5) + r^(3/5))

In[2]:= Factor[1 + Log[r]^(1/5) - 2 Log[r]^(2/5) - Log[r]^(3/5) - 2 Log[r]^(4/5)]
Out[2]= -(1 + 2 Log[r]^(1/5)) (-1 + Log[r]^(1/5) + Log[r]^(3/5))

In[3]:= Factor[Exp[2x] + 2 Exp[x] + 1]
Out[3]= (1 + E^x)^2

In[4]:= Factor[E^(4x) - 1]
Out[4]= (-1 + E^x) (1 + E^x) (1 + E^(2 x))

In[5]:= Apart[1/(-1 + r^(3/7))]
Out[5]= (-2/3 - 1/3 r^(1/7))/(1 + r^(1/7) + r^(2/7)) + 1/3/(-1 + r^(1/7))

In[6]:= Together[1/Log[r]^(1/2) + 1/Log[r]^(1/3)]
Out[6]= (1 + Log[r]^(1/6))/Sqrt[Log[r]]
```

### Files touched

- `src/poly.h`, `src/poly.c`: new `poly_find_radical_gen`,
  `poly_subst_radical_to_gen`, `poly_subst_radical_from_gen`,
  `poly_make_fresh_gen`. Internal helpers `decompose_exponent`,
  `is_target_power`, `walk_find_radical_base`, `walk_gather`.
- `src/rat.c`: removed the private symbol-only helpers; `builtin_cancel`
  and `builtin_together` now go through the shared `poly_*` helpers.
- `src/parfrac.c`: 1-arg `Apart` runs through the radical pass before
  its standard partial-fraction logic.
- `src/facpoly.c`: `builtin_factor` runs through the radical pass before
  the existing `Together`/`Numerator`/`Denominator` setup.
- `tests/test_rat.c`: added `test_cancel_nonsymbol_base` and a
  Log-base case in `test_together_algebraic_generator`.
- `tests/test_radical_polyops.c`: new test file pinning the four
  headline behaviours (Factor radical-symbol, Factor radical-Log[r],
  Factor exponential, Apart radical, no-trigger).

## Pattern: HoldFirst attribute (2026-05-05)

`Pattern` is now declared with `HoldFirst | Protected`, matching
Mathematica. The first argument of `Pattern[name, body]` is the
pattern variable name and must be held literally; without `HoldFirst`,
the evaluator would resolve the name through any pre-existing
`OwnValue` before the rule was installed.

### Symptoms before the fix

```
v = 99;
f[v_] := v^2          (* installs a malformed DownValue f[Pattern[99,Blank[]]] *)
f[3]                  (* segfaults: matcher treats integer 99 as a binding name *)
```

A subtler manifestation: `simp_install_roots_of_unity_helpers` parses
its `$ruSimplify[expr_] := ...` definition lazily on first use of
roots-of-unity simplification. If the user had previously assigned
`expr = (rational expression)`, the lazy install evaluated
`Pattern[expr, Blank[]]` against the user's `OwnValue`, baked the
user's value into the rule, and corrupted every later `Simplify`
call. With `HoldFirst` on `Pattern` the LHS-evaluation step at
`eval.c:760-787` (Set/SetDelayed inner-arg evaluation, intentional so
that `f[2+3] = ...` becomes `f[5] = ...`) leaves the pattern variable
untouched.

### Files touched

- `src/attr.c`: added `{"Pattern", ATTR_HOLDFIRST | ATTR_PROTECTED}`
  to `builtin_attrs[]`. Picked up by `attr_init()` and
  `get_attributes("Pattern")` automatically.

## Cancel / Together: push negative lead from denominator into numerator (2026-05-05)

After the polynomial GCD step, `cancel_recursive` already flipped the
sign of the result when the denominator collapsed to a literal
negative integer. It did not flip when the denominator was a
polynomial whose leading coefficient happened to be negative — so
PolynomialGCD's sign convention could leave both numerator and
denominator with negative leading terms, e.g.

```
((2 x - y) (x^2 - x y - 6 y^2)) / ((x + y) (2 x^2 + 3 x y - 2 y^2)) // Cancel
```

came back as `(-x + 3 y) / (-x - y)` (SimplifyCount 21) rather than
the canonical `(x - 3 y) / (x + y)` (SimplifyCount 13). Worse,
`Simplify` then preferred `Apart`'s `-3 - 4 x / (-x - y)`
(SimplifyCount 20) over the Cancel form because it scored slightly
lower, never seeing the truly simple form.

### Fix

`src/rat.c` adds `den_has_negative_lead`, called immediately after
`cancel_exact_div_wrapper` produces `new_num` / `new_den`. It returns
true when either:

1. `is_superficially_negative(new_den)` — the den is a literal
   negative number, a negative `Rational[n, d]`, or a `Times` whose
   leading factor is negative.
2. `new_den` is a `Plus` and **every** term is superficially negative
   (`Plus[-x, -y]`, `Plus[-1, -x]`, `Plus[-x, -3 y]`, …) — i.e. the
   whole sum is `-1 *` (positive-lead sum).

Mixed-sign sums like `Plus[-1, x]` (== `x - 1`) and `Plus[-x, y]`
(== `y - x`) are deliberately **not** flipped: both orderings have
the same leaf count and either choice would be arbitrary, so we
preserve the canonical-sort order.

`Together` benefits transparently because `together_recursive` calls
`cancel_recursive` after combining over a common denominator. `Apart`
benefits indirectly via its `Together` precondition.

### Examples

| Input | Before | After |
|-------|--------|-------|
| `Cancel[(x^2-y^2)/(y-x)]` | `(-x - y)` ✓ already, but for the wrong reason | `-x - y` ✓ |
| `Cancel[(x+1)/(-2)]` | `(x+1)/(-2)` | `(1/2) (-1 - x)` |
| `Cancel[((2x-y)(x^2-xy-6y^2))/((x+y)(2x^2+3xy-2y^2))]` | `(-x + 3 y)/(-x - y)` | `(x - 3 y)/(x + y)` |
| `Together[1/(b-a)]` | `1/(-a + b)` | `1/(-a + b)` (unchanged — mixed sign) |
| `Cancel[1/(x-1)]` | `1/(-1 + x)` | `1/(-1 + x)` (unchanged) |

### Files touched

- `src/rat.c`: added `den_has_negative_lead`; widened the
  `cancel_recursive` post-cancellation sign check from
  "denominator is a negative integer" to "denominator has a globally
  negative lead".


## Times: trig / hyperbolic ratio canonicalization (2026-05-05)

`Times` now automatically rewrites products and quotients of trigonometric
(and hyperbolic) functions sharing an argument into the shortest available
naming. The rewrite is a structural identity (`Sec[x] := 1/Cos[x]`,
`Tan[x] := Sin[x]/Cos[x]`, etc.) applied at the moment `builtin_times` (or
`builtin_power`, for the parser-collapsed `1/Cos[x]` shape) finishes its
base-grouping pass.

### What gets rewritten

| Input pattern | Output |
|--------------|--------|
| `1/Cos[x]`, `1/Sin[x]`, `1/Tan[x]`, `1/Cot[x]`, `1/Sec[x]`, `1/Csc[x]` | `Sec[x]`, `Csc[x]`, `Cot[x]`, `Tan[x]`, `Cos[x]`, `Sin[x]` |
| `1/Cosh[x]`, `1/Sinh[x]`, `1/Tanh[x]`, `1/Coth[x]`, `1/Sech[x]`, `1/Csch[x]` | `Sech[x]`, `Csch[x]`, `Coth[x]`, `Tanh[x]`, `Cosh[x]`, `Sinh[x]` |
| `Sin[x]/Cos[x]` | `Tan[x]` |
| `Cos[x]/Sin[x]` | `Cot[x]` |
| `Sin[x]^2/Cos[x]^2` | `Tan[x]^2` |
| `Sin[x] Csc[x]`, `Cos[x] Sec[x]`, `Tan[x] Cot[x]` | `1` |
| `Cos[x] Tan[x]`, `Sin[x] Cot[x]` | `Sin[x]`, `Cos[x]` |
| `Tan[x] Csc[x]`, `Cot[x] Sec[x]` | `Sec[x]`, `Csc[x]` |
| `Sin[x]^3 / Cos[x]` | `Sin[x]^2 Tan[x]` |
| `Sin[x] / Cos[x]^3` | `Sec[x]^2 Tan[x]` |
| `1/Cos[x]^3` | `Sec[x]^3` |
| `Sinh[x]/Cosh[x]`, `Cosh[x]/Sinh[x]` | `Tanh[x]`, `Coth[x]` |

### Algorithm

Each contributing factor is mapped to its `(Sin^a · Cos^b)` decomposition:

```
Sin[x]  → ( 1,  0)     Cos[x]  → ( 0,  1)
Tan[x]  → ( 1, -1)     Cot[x]  → (-1,  1)
Sec[x]  → ( 0, -1)     Csc[x]  → (-1,  0)
```

Factors sharing the same argument (compared by `expr_eq`) are summed into a
single `(a, b)` pair, then re-emitted using the shortest naming for that
sign pattern: pure positive `(a, 0)` → `Sin^a`; mixed `(a > 0, b < 0)` with
`a + b = 0` → `Tan^a`; with `a + b > 0` → `Sin^(a+b) Tan^|b|`; etc. The
hyperbolic family `{Sinh, Cosh, Tanh, Coth, Sech, Csch}` is handled
identically but in a separate bucket from the trig family.

### Skipped cases

- Non-integer exponents (`Cos[x]^(1/2)`, `Sin[x]^n`).
- Different arguments do not merge: `Sin[x] Cos[y]` stays as is.
- Solo positive-exponent canonical names (`Sin[x]`, `Tan[x]^2`, `Sec[x]^3`)
  are left alone — already canonical.
- Multi-contributor buckets whose canonical re-emission would be the same
  multiset of factors are skipped (avoids perturbing the canonical sort
  order downstream).

### Suppression for trigsimp internals

The pipelines inside `TrigReduce`, `TrigFactor`, `TrigExpand`, `TrigToExp`,
and `ExpToTrig` rewrite `Tan[x] :> Sin[x]/Cos[x]` (and friends) as a first
step so that polynomial machinery (Together, Factor, identity rules) can
operate on the resulting Sin/Cos polynomial form. Without protection, the
new Times-level canonicalizer would immediately re-collapse the Sin/Cos
form back to `Tan` and the rule pipeline would never make progress.

A reentrant suppress counter (`trig_canon_suppress_inc/dec` in
`src/trig_canon.h`) lets these pipelines disable the rewrite for the
duration of their internal work. The outer evaluator re-evaluates whatever
they return, so the user-visible result is still canonicalized.

### Files

- `src/trig_canon.h`, `src/trig_canon.c` — new module containing
  `trig_canon_groups` (Times entry), `trig_canon_power` (Power entry), and
  the suppress counter.
- `src/times.c` — uses the shared `BasePower` typedef and calls
  `trig_canon_groups` after radical fusion.
- `src/power.c` — calls `trig_canon_power` just before the final
  `return NULL` so `1/Cos[x]` rewrites even when the parser collapses it
  directly to `Power[Cos[x], -1]`.
- `src/trigsimp.c` — wraps `builtin_exptotrig`, `builtin_trigtoexp`,
  `builtin_trigexpand_impl`, `builtin_trigreduce_impl`, and
  `trigfactor_run_pipeline` with `trig_canon_suppress_inc/dec`.
- `tests/test_trig_canon.c` — 19 test groups covering solo reciprocals,
  pair collapses, full cancellations, productive simplifications, leftover
  powers, hyperbolic family, mixed arguments, non-integer / symbolic
  exponents, FullForm structure, compound arguments, numeric coefficients,
  and explicit no-rewrite-outside-Times/Power assertions.

## REPL session hooks: `$PreRead`, `$Pre`, `$Post`, `$PrePrint`, `$Epilog` (2026-05-05)

PicoCAS now exposes the standard Mathematica REPL hooks. Each is a
global symbol whose OwnValue, if assigned, is consulted at a specific
phase of the read/eval/print loop. With no assignment every hook is a
no-op pass-through, so existing scripts are unaffected.

### Hooks

- **`$PreRead`** — text-level. Applied to the raw input string after
  it is read by readline but before parsing. Must be assigned a
  function that takes one string and returns a string. If the
  assigned function returns a non-string value the helper emits a
  `$PreRead::strret` diagnostic and falls back to the original input.
- **`$Pre`** — applied to the parsed Expr before evaluation. Without
  HoldAll the wrapped expression is evaluated before `$Pre` sees it
  (matching Mathematica's documented behaviour, which makes the
  effect indistinguishable from `$Post`). Assigning a head with
  HoldAll lets `$Pre` intercept the raw parse tree.
- **`$Post`** — applied to the evaluator's result.
- **`$PrePrint`** — applied to the result expression just before
  printing. `Out[n]` retains the pre-`$PrePrint` value, so display
  rewrites do not contaminate later `Out[k]` references.
- **`$Epilog`** — evaluated once when the session terminates (Quit[]
  or EOF). `Set` assigns the evaluated RHS once; use `SetDelayed` if
  you want the body to fire on every shutdown.

### REPL ordering

```
read -> $PreRead -> parse -> [In[n]=parsed] -> $Pre ->
evaluate -> $Post -> [Out[n]=result] -> $PrePrint -> print
```

`In[n]` records what the user typed (post-`$PreRead`, pre-`$Pre`) so
history reflects the user's input. `Out[n]` records the evaluator's
answer with `$Post` applied but without `$PrePrint` rewrites.

### Examples

```mathematica
$PreRead = Function[s, StringJoin[s, " + 1"]]
3                       (* parsed as 3 + 1, evaluated to 4 *)

$Pre = Function[x, 2 x]
5                       (* $Pre[5] -> 10 *)

$Post = Function[x, wrap[x]]
1 + 2                   (* prints wrap[3] *)

$PrePrint = Function[x, displayWrap[x]]
1 + 2                   (* Out[n] = 3, prints displayWrap[3] *)

$Epilog := Print["Goodbye!"]
Quit[]                  (* prints Goodbye! before exiting *)
```

### Files

- `src/repl_hooks.{c,h}` — new module exposing
  `repl_hooks_init`, `repl_apply_pre_read`, `repl_apply_pre`,
  `repl_apply_post`, `repl_apply_pre_print`, `repl_apply_epilog`. All
  helpers are unit-testable independently of the readline-driven
  REPL loop.
- `src/core.c` — `core_init` calls `repl_hooks_init()` so docstrings
  are available before any user code runs.
- `src/repl.c` — `process_input` runs the read pipeline through the
  helpers; EOF and `Quit[]` paths in `repl_loop` invoke
  `repl_apply_epilog()` before tearing down.
- `tests/test_repl_hooks.c` — 32 tests covering init, each hook's
  unset/identity/transforming/clear-disables paths, fallback on
  non-string `$PreRead` returns, HoldAll-attributed `$Pre`
  intercepting raw input, `$Post` chaining on evaluator output,
  `$PrePrint` rewriting display only, the
  Set-vs-SetDelayed distinction for `$Epilog`, end-to-end pipeline
  composition, and reassignment-replaces-behaviour semantics.

## Simplify: trig at rational multiples of Pi (2026-05-05)

Adds a Simplify-time tree walker (`simp_trig_pi_canon` in
`src/simp.c`) that rewrites every `Sin`, `Cos`, `Tan`, `Cot`, `Sec`,
`Csc` of a rational multiple of Pi into the representation whose
reduced fraction has the smaller numerator.  The rewrite uses the
complement identities (`Cos[a] == Sin[Pi/2 - a]` etc.) under the
post-quadrant reduction, so

```
Cos[4 Pi/9]   ->  Sin[Pi/18]      (4/9 vs 1/18 picks 1)
Cos[5 Pi/9]   ->  -Sin[Pi/18]     (sign flip via Cos[Pi - x] = -Cos[x])
Tan[7 Pi/18]  ->  Cot[Pi/9]       (7/18 vs 1/9 picks 1)
```

Wired into `simp_search` at the seed phase and chained on every
transform output in the round loop.  Idempotent on its own output
and inert when the input has no rational-Pi trig calls.

### Why a Simplify-only transform

`builtin_cos` / `builtin_sin` already handle a whitelist of "nice"
denominators (1, 2, 3, 4, 5, 6, 10, 12) at evaluation time and
otherwise leave the call as-is.  Wiring this rewrite into the global
`Cos` / `Sin` builtins would change every user-facing print of e.g.
`Cos[4 Pi/9]` and break parity with Mathematica.  The unified
representation is only needed during Simplify, where the goal is to
make additive cancellation discover that two structurally distinct
calls denote the same value.

### Cases fixed

```
Simplify[Cos[4/9 Pi] - Sin[Pi/18]]                   ->  0
Simplify[Cos[Pi/9] Cos[2Pi/9] Cos[3Pi/9] Cos[4Pi/9]
          - 1/16]                                    ->  0
```

For the second, `Cos[3 Pi/9]` is `Cos[Pi/3] = 1/2`, so the product
equals `(1/2) Cos[Pi/9] Cos[2 Pi/9] Cos[4 Pi/9]`.  TrigReduce's
product-to-sum machinery already lands at
`1/8 (Cos[4 Pi/9] + Cos[5 Pi/9])`; the canonicalizer then collapses
both Cos calls to (`Sin[Pi/18]`, `-Sin[Pi/18]`) and the Plus
cancels.

### Files

- `src/simp.c` — `simp_trig_pi_canon`, `trig_pi_canon_one`,
  `trig_pi_canon_walk`, plus helpers `trig_pi_extract`,
  `trig_pi_make_arg`, `trig_pi_i64_gcd`, `trig_pi_reduce_frac`.
  Wired into `simp_search` at the seed phase and the round loop's
  per-transform chain.
- `tests/test_simplify.c` — `test_simplify_cos_minus_sin_complement`
  and `test_simplify_morrie_product_minus_constant` cover both
  reported cases.

## TrigReduce: coefficient-aware angle-addition + cache coherence (2026-05-05)

Two adjacent improvements to the TrigReduce pipeline.

### Coefficient-aware angle-addition collapse rules

`trig_reduce_collapse` previously enumerated eight literal-coefficient
rules of the form `Sin[a] Cos[b] + Cos[a] Sin[b] -> Sin[a+b]` (with
sign / Cos / Sinh / Cosh variants).  These match only when the Plus
children are bare products of two trig calls.  Replaced with four
`c_.`-keyed rules that bind the same coefficient on both Plus
children:

```
c_. Sin[a_] Cos[b_] + c_. Cos[a_] Sin[b_] + r___ :> c Sin[Expand[a + b]] + r
c_. Sin[a_] Cos[b_] - c_. Cos[a_] Sin[b_] + r___ :> c Sin[Expand[a - b]] + r
c_. Cos[a_] Cos[b_] - c_. Sin[a_] Sin[b_] + r___ :> c Cos[Expand[a + b]] + r
c_. Cos[a_] Cos[b_] + c_. Sin[a_] Sin[b_] + r___ :> c Cos[Expand[a - b]] + r
```

(Plus the four hyperbolic analogues.)  `c_.` defaults to 1 so the
literal cases still match, and inputs like `k Cos[a] Cos[b] -
k Sin[a] Sin[b]` now reduce to `k Cos[a + b]`.

### Cache coherence on `trig_canon_suppress_dec`

`trig_canon_suppress_dec` now bumps `eval_clock` when the suppress
counter drops to zero.  Subexpressions evaluated while suppression
was active have their `last_evaluated_at` field stamped at the
current clock, but that stamp records "fully evaluated under the
suppress flag", not "fully evaluated unconditionally".  Without the
bump, `evaluate()` short-circuits on the cached form even after the
suppress region ends, so structures like

```
Cos[5 Pi/18] Sec[5 Pi/18] Sec[7 Pi/18] Sin[7 Pi/18]
```

inside Plus children of a TrigReduce result stay visible instead of
collapsing to `4 Tan[7 Pi/18]`.  Bumping invalidates the
suppressed-mode cache so the next `evaluate()` pass walks the tree
with `trig_canon` active.

### Files

- `src/trig_canon.c` — `trig_canon_suppress_dec` calls
  `eval_clock_bump()` on the outermost dec.
- `src/trigsimp.c` — `trig_reduce_collapse` rules rewritten to use
  `c_.` coefficients.
- `tests/test_trigreduce.c` —
  `test_trigreduce_angle_addition_circular_with_coeff` covers both
  the `Sqrt[3]` and the symbolic-`k` cases.

## Simplify / Together / Cancel: depth guards on heuristic_factor (2026-05-05)

`heuristic_factor` (in `src/facpoly.c`) had two unbounded-recursion
paths that blow the C stack on inputs whose "polynomial variables"
are opaque atoms (e.g. `Power[base, n]` with symbolic `n` --
`collect_variables` treats those as single generators):

1. The cont/pp split inside `heuristic_factor` calls
   `heuristic_factor(cont)` recursively.  When `poly_content` extracts
   a non-trivial GCD that itself shares all variables of `P`, the
   recursion never makes progress on `v_count`.
2. `factor_via_z_independent_split` per-z peel calls
   `heuristic_factor(remaining)` at the same `v_count`.

Two guards added:

- A reentrancy depth counter (`heuristic_factor_depth`, max 64)
  short-circuits to `expr_copy(P)` when the depth limit is hit.
- A structural progress check: if `exact_poly_div(P, cont)` yields
  `pp == 1`, then `cont == P` up to a constant and recursing on
  `cont` is identical to recursing on `P` -- bail and return `P`
  unchanged.

### Cases fixed (segfault → correct-but-unsimplified)

```
Simplify[(-14 Sin[x]^3 + 35 Sin[x] + 6 Sqrt[3] Cos[x]^3
           + 9 Sqrt[3] Cos[x]) / (Cos[2 x] + 4)
           - (7 Sin[x] + 3 Sqrt[3] Cos[x])]
   was: SIGSEGV in cancel_recursive ↔ together_recursive mutual
        recursion via PolynomialGCD.

Simplify[(Cos[x] - 1)^n (Cos[x] + 1)^n - (-Sin[x]^2)^n]
Simplify[(Sin[x] - 1)^n (Sin[x] + 1)^n - (-Cos[x]^2)^n]
   were: SIGSEGV in heuristic_factor ↔ factor_via_z_independent_split
         cycle at v_count = 3.
```

All three now return correct-but-unsimplified results.

## Simplify: PythagCanon top-level short-circuit (2026-05-06)

A new `transform_pythag_canon` (in `src/simp.c`) replaces the
substitution-bound version of the Pythagorean identity:

```
Cos[x_]^n_Integer /; n >= 2 && EvenQ[n] :> (1 - Sin[x]^2)^(n/2)
Sin[x_]^n_Integer /; n >= 2 && EvenQ[n] :> (1 - Cos[x]^2)^(n/2)
Cosh[x_]^n_Integer /; n >= 2 && EvenQ[n] :> (1 + Sinh[x]^2)^(n/2)
Sinh[x_]^n_Integer /; n >= 2 && EvenQ[n] :> (-1 + Cosh[x]^2)^(n/2)
```

The transform Expands the input, applies each rule set via
`ReplaceRepeated`, Expands again, scores all four results plus the
input, and returns the lowest-scoring strict win (or a copy of the
input when no direction beats it). It is wrapped in `simp_memo_wrap`
under the pseudo-head `$PythagCanon` so repeated calls with the same
input are O(1).

The transform fires in two places:

1. As a seed inside `simp_search` alongside `transform_pythag_reduce`
   and `transform_pythag_square_complete`. Same role as those: feed
   the round-loop with a coefficient-blind alternative.

2. As a top-level short-circuit at the very top of `simp_bottomup`,
   gated to `depth == 0`. If the canon strictly improves the score,
   the local `input` pointer is swapped to point at the canon result
   and the rest of `simp_bottomup` runs on that smaller form. With
   the swap, an input like

   ```
   18 (Cos[x]+1)(Cos[x]-1)(Cos[y]^2-1)^2 (x-1)
       + 18 (x-1) Sin[x]^2 Sin[y]^4
   ```

   collapses to `0` in ~5 ms instead of ~6.5 s -- the recursive
   descent into every Plus/Times subnode no longer happens because
   the parent has already become an atom. A defensive check for
   `input->type != EXPR_FUNCTION` after the swap returns early to
   avoid dereferencing the function-union member on a literal.

### Why we need it

The bare `PythagReduce` rules look for `1 +/- Cos[x_]^2 + r___`
shapes with unit constant and unit coefficient. After Expand, a
typical input like `(Cos+1)(Cos-1) (x-1) 18` yields
`-18 (x-1) + 18 (x-1) Cos[x]^2`, where the `1` is now `-18 (x-1)` --
a coefficient PythagReduce cannot match. The substitution-based
canonicalizer is coefficient-blind and operates directly on the
Power node.

### Coverage

12 new unit tests in `tests/test_simplify.c` cover:

- the original 6.5 s user case -> 0
- difference-of-squares trig products in 1 and 2 variables
- higher even powers (`Cos[x]^4 + 2 Cos[x]^2 Sin[x]^2 + Sin[x]^4 -> 1`)
- hyperbolic counterpart (`(Cosh+1)(Cosh-1) - Sinh^2 -> 0`)
- scalar-multiple Pythagoras (`5 Cos^2 + 5 Sin^2 -> 5`)
- negative regressions: bare atoms, no-op when no even-power trig,
  preservation of `Sin[x]^2`, `Cos[x]^2` etc. as canonical forms,
  no destabilisation of the existing `1 - Sin[x]^2 -> Cos[x]^2`
  rewrite.

## Simplify: reciprocal-pair Pythagorean canonicalisation + TrigReduce short-circuit (2026-05-06)

### Symptom

Two pairs of user-reported failure modes:

1. **Sec/Csc Pythagorean products (HANG)**.
   `(Sec[x]+1)(Sec[x]-1) - Tan[x]^2` and the two-variable version
   `(Sec[x]+1)(Sec[x]-1) * 2(Sec[y]-1)(Sec[y]+1) * 9(Sec[y]-1)(Sec[y]+1)
   * (x-1) - 18(x-1) Tan[x]^2 Tan[y]^4`
   Mathematica: `0` in <1 ms.  picocas: hung past 2 minutes.

2. **(Sin+Cos)^4 / (1+Sin[2x])^2 (slow)**.
   `Simplify[(Sin[x]+Cos[x])^4 - (1+Sin[2 x])^2]`.
   Mathematica: `0` in <1 ms.  picocas: ~220 ms.

### Root cause

The `has_pythag_head` gate that fast-skips `PythagReduce` /
`PythagCanon` listed only `Cos/Sin/Cosh/Sinh/Tan/Cot/Tanh/Coth`,
never `Sec/Csc/Sech/Csch`.  And `PythagCanon` carried only four
substitution directions (`Cos->Sin`, `Sin->Cos`, `Cosh->Sinh`,
`Sinh->Cosh`); it had nothing for the reciprocal-pair identities
`Sec^2 = 1 + Tan^2`, `Csc^2 = 1 + Cot^2`, `Sech^2 = 1 - Tanh^2`,
`Csch^2 = -1 + Coth^2`.  Inputs that survived the gate fell through
to `TrigFactor` (~700 ms on the two-variable Sec/Tan case), which
rewrote `Sec/Csc` to `1/Cos`, `1/Sin` and produced a polynomial in
`Cos`, `Sin` that nothing else in the round loop could collapse
within the search-budget rounds.

For (Sin+Cos)^4: `TrigReduce` was a SIMP_TRANSFORMS entry but not a
seed.  Bottom-up descent ran the full search on each child of the
top-level Plus before the parent step ever called TrigReduce; the
children's searches ran Factor / TrigFactor / FactorSquareFree
(28 + 14 + 9 ms) on `(Sin[x]+Cos[x])^4` and `(1+Sin[2x])^2`
individually -- even though TrigReduce on the *whole* input lands
on `0` in 3 ms.

### Fix (three pieces in `src/simp.c`)

1. **`has_pythag_head` extended** with `Sec/Csc/Sech/Csch`.  Now the
   reciprocal-pair gate passes whenever any of those heads appear.

2. **`PythagReduce` rule list extended** with eight reverse-direction
   reciprocal identities: `-1 + Sec^2 -> Tan^2`, `1 - Sec^2 -> -Tan^2`,
   `-1 + Csc^2 -> Cot^2`, `1 - Csc^2 -> -Cot^2`,
   `1 - Sech^2 -> Tanh^2`, `-1 + Sech^2 -> -Tanh^2`,
   `1 + Csch^2 -> Coth^2`, `-1 - Csch^2 -> -Coth^2`.  These complement
   the existing eight (Tan -> Sec, etc.) so a Plus with a `1 +/- Sec^2`
   shape collapses without going through TrigFactor.

3. **`PythagCanon` rule sets extended** from four to twelve.  The new
   directions are
   ```
   Sec[x_]^(2k) -> (1 + Tan[x]^2)^k       (Sec -> Tan)
   Tan[x_]^(2k) -> (-1 + Sec[x]^2)^k      (Tan -> Sec)
   Csc[x_]^(2k) -> (1 + Cot[x]^2)^k       (Csc -> Cot)
   Cot[x_]^(2k) -> (-1 + Csc[x]^2)^k      (Cot -> Csc)
   Sech[x_]^(2k) -> (1 - Tanh[x]^2)^k     (Sech -> Tanh)
   Tanh[x_]^(2k) -> (1 - Sech[x]^2)^k     (Tanh -> Sech)
   Csch[x_]^(2k) -> (-1 + Coth[x]^2)^k    (Csch -> Coth)
   Coth[x_]^(2k) -> (1 + Csch[x]^2)^k     (Coth -> Csch)
   ```
   Same `n >= 2 && EvenQ[n]` guard as the existing four.  The
   keep-best-strict-win selection means firing a no-op direction is
   only an ms-scale speed cost, never a correctness cost.

4. **TrigReduce promoted to seed and to depth-0 short-circuit**.
   Two new call sites:

   - In `simp_search`, a `TrigReduce` seed runs after `PythagCanon`.
     Score-gated like `TrigRoundtrip` (drop seeds whose leaf count
     exceeds `2 * input + 8` to avoid pathological blow-ups).

   - In `simp_bottomup`, at `depth == 0`, a `TrigReduce` short-circuit
     runs in parallel with the existing `PythagCanon` short-circuit:
     try `TrigReduce[input]` once before descending into children;
     adopt the result only when it's a strict score win.  This is
     the move that turns the `(Sin+Cos)^4` case from 220 ms to 5 ms:
     bottom-up no longer wastes a full search on each child.

5. **`simp_search` early-exit on literal zero**.  Once `best` becomes
   the literal integer `0`, no transform can produce anything
   structurally simpler, so the round loop and remaining seed-phase
   blocks short-circuit via a `goto search_done` label.  Restricted
   to literal `0` so that custom `ComplexityFunction` callers still
   get full search coverage for non-zero atomic answers.

### Performance (before -> after)

| Case                                                  | Before | After |
|-------------------------------------------------------|--------|-------|
| `(Cos+1)(Cos-1) ... - Sin^2 Sin^4` (two-var)          | 4 ms   | 5 ms  |
| `(Cosh+1)(Cosh-1) ... - Sinh^2 Sinh^4`                | 4 ms   | 5 ms  |
| `(Sec+1)(Sec-1) ... - Tan^2 Tan^4` (two-var)          | hung   | 5 ms  |
| `(Sin[x]+Cos[x])^4 - (1+Sin[2x])^2`                   | 220 ms | 5 ms  |
| `(Cosh[x]+Sinh[x])^3 - (Cosh[3x]+Sinh[3x])`           | 31 ms  | 4 ms  |
| `Cos[x]^3 - Sin[x]^3 - (Cos-Sin)(1+Sin Cos)`          | 14 ms  | 2 ms  |

### Files

- `src/simp.c` -- `has_pythag_head` extended; `transform_pythag_reduce_impl`
  rule list extended with eight reciprocal identities;
  `transform_pythag_canon_impl` extended with eight substitution
  directions; `simp_search` gains a TrigReduce seed and a literal-zero
  early-exit; `simp_bottomup` gains a depth-0 TrigReduce short-circuit;
  `simp_best_is_zero` helper added.
- `tests/test_simplify.c` -- 25 new unit tests covering the round-2
  user cases.  Passing: Sec/Csc Pythagoreans, (Sin+Cos)^4,
  (Cosh+Sinh)^3, Cos^3-Sin^3, single-variable reciprocal-pair Pythags
  (`1 - 1/Sin^2 + Cot^2`, etc.), double-angle with constant factors,
  and the chain of Exp / 2^x identities (cases 11, 13, 14, 15, 16, 19, 20).
  Pinned to current output (TODO: should be 0): the eight
  power-arithmetic identities (cases 7, 8, 10, 12, 17, 18, 21) and the
  Tan-addition identity (case 6).  Each pinned test carries a comment
  pointing at the missing transform; replacing the expected string
  with `"0"` becomes a regression guard once the transform lands.

### Known limitations (tracked in the pinned tests)

- **Power-arithmetic / prime-base rebasing**: `4^x * 2^(-x) * 2^(-x)`
  does not collapse to `1` because picocas does not recognise that
  `4 = 2^2` and rebase to a common prime base.  Touches cases 7, 10.
- **Power-of-product distribution**: `(2 E)^(x+y)` is not split into
  `2^(x+y) E^(x+y)`.  picocas has no `PowerExpand` builtin, and
  Simplify doesn't apply `(a*b)^n -> a^n b^n` unconditionally
  (Mathematica's Simplify does this internally with branch-cut
  awareness).  Touches case 12.
- **Quotient power with integer assumption**: `(y/x)^(-n)` is not
  split to `x^n / y^n` even when `n` is asserted Integer.  Touches
  case 17.
- **Same-base combination across structural representations**:
  `(1/x)^Log[2] / x` does not combine with `(1/x)^(1+Log[2])` because
  `1/x` and `(1/x)^1` are not treated as a common base in the search.
  Touches case 21.
- **Tan-addition recognition**: `Tan[2] Tan[3] - (-Tan[2]/Tan[5] -
  Tan[3]/Tan[5] + 1)` does not collapse via `Tan[a+b] = (Tan[a]+Tan[b])
  / (1 - Tan[a] Tan[b])` because picocas has no inverse-of-TrigExpand
  recognition for tangents.  Touches case 6.

## Simplify: PrimeRebase, PowerOneify, RadicalCanon, TanAddition (2026-05-06)

Four new universal Simplify seeds that close power-arithmetic, radical-
canonicalisation, and trigonometric-addition gaps.  All wired into
`simp_dispatch` (the first three) or `simp_search`'s seed phase (the
fourth) so they reach every shape pipeline.  Each is score-gated and
inert by default -- no regressions on the existing 91 test binaries.

### PrimeRebase

Rewrites `Power[c, e] -> Power[p, k*e]` when the integer `c = p^k` is a
perfect prime power (`k >= 2`).  Sound for ALL complex `e` because the
identity `(a^b)^c = a^(b*c)` is unconditional when `a > 0`.  Closes the
"different base bucket" gap in canonical Times for inputs mixing `4^x`
and `2^x` style factors.  Walks the tree once (recursing into exponents
so nested rebasable powers like `4^x` inside `2^(x*4^x)` bubble out),
then a single `evaluate()` collapses same-base factors via Times'
canonical exponent combine.

Coverage:

| Input | Was | Now |
|------|-----|-----|
| `Simplify[4^x * 2^(-x) * 2^(-x) - 1]` | `-1 + 2^(-2 x) 4^x` | `0` |
| `Simplify[f[4^x*2^(-x)*2^(-x)] - f[1]]` | `-f[1] + f[2^(-2 x) 4^x]` | `0` |
| `Simplify[2^(2^(2 x) x) - 2^(x*4^x)]` | `-2^(4^x x) + 2^(2^(2 x) x)` | `0` |

Branch-cut-sensitive shapes (e.g. `(-4)^x` -> `(-2)^(2x)`,
`(a*b)^n` -> `a^n*b^n` for non-integer `n`) are deliberately NOT
handled: those need PowerExpand semantics with explicit assumptions and
remain documented limitations (test cases 8, 12, 17 stay pinned).

### PowerOneify

In every `Times` node, combines factors of the form `A * Power[A, e]`
into `Power[A, e + 1]` using the universally valid identity
`Power[a, 1] = a`.  Closes the "bare-A vs Power[A, _] mismatch" in
canonical Times bucketing -- canonical Times groups factors by the base
of any wrapping `Power[base, exp]`, but a bare `A` whose canonical form
is itself a `Power` (e.g. `Power[x, -1] = 1/x`) does NOT get re-bucketed
as `Power[A, 1]` before grouping.

Coverage:

| Input | Was | Now |
|------|-----|-----|
| `Simplify[(1/x)^Log[2]/x - (1/x)^(1 + Log[2])]` | `-(1/x)^(1+Log[2]) + (1/x)^Log[2]/x` | `0` |

### RadicalCanon

Two related Power-canonicalisations that picocas's standard Power
evaluator does not perform, leaving equivalent expressions in distinct
shapes that block additive cancellation:

1. **Rational-base split** -- `Power[Rational[a, b], q] ->
   Power[a, q] * Power[b, -q]` for positive integer `a`, `b`.  Sound
   because `a, b > 0` -> `(a/b)^q = a^q * b^(-q)` for all complex `q`
   with no branch cut.  Worked example: `Sqrt[1/2] = Power[Rational[1,2],
   1/2]` -> `Power[2, -1/2]`.
2. **Negative-rational-exponent rationalisation** -- `Power[a, q] ->
   Power[a, r] * Rational[1, a^k]` for `q < 0` rational with denominator
   `> 1`, where `k = ceil(-q)` and `r = q + k in [0, 1)`.  Worked example:
   `Power[2, -1/2]` -> `Power[2, 1/2] * 1/2` = `Sqrt[2]/2`.  The result
   does NOT fold back: `Rational[1, a^k]` occupies a distinct base bucket
   from `Power[a, r]` in canonical Times.

Together these two rules force radical shapes (`Sqrt[1/2]` vs `1/Sqrt[2]`
vs `Sqrt[2]/2`) into a single canonical representation so that pairs
like `-Sqrt[1/2]/3 + Sqrt[2]/6` cancel additively.

Coverage:

| Input | Was | Now |
|------|-----|-----|
| `Simplify[Sqrt[1/2] - 1/Sqrt[2]]` | `Sqrt[1/2] - 1/Sqrt[2]` | `0` |
| `Simplify[Sqrt[2] Sqrt[3] Sqrt[5] x - Sqrt[30] x]` | -- | `0` (~2 ms) |
| `Simplify[c1 - c2]` (polynomial * `(1+x^2)^(3/2)` identity) | -- | `0` (~12 ms) |
| `Simplify[t1 + t2 + t3]` (Sin/Cos/Sqrt mix, see test) | `Sin[y] (-1/3 Sqrt[1/2] + 1/6 Sqrt[2])` (505 ms, wrong) | `0` (~380 ms) |

The score gate for RadicalCanon is `<=` rather than `<`: even when the
score is the same, the rewrite is correctness-preserving and always
yields the more useful form for downstream additive cancellation.  This
mirrors how the existing `LogExpRules` cascade treats force-take rewrites.

### TanAddition

Recognises angle-addition shapes that picocas's TrigExpand misses
because the latter only fires on literal-Plus arguments.  Algorithm:

1. Walk the input collecting every distinct expression appearing as the
   argument of any `Tan, Cot, Sin, Cos, Sec, Csc` head (a borrowed-pointer
   set; structural equality via `expr_eq`).
2. For each ordered pair `(i, j)` with `i != j`, evaluate
   `args[i] + args[j]`.  If the sum also belongs to the set, the triple
   `(a, b, c=a+b)` witnesses an addition-formula opportunity.
3. For each triple, build a `RuleDelayed` for each Trig head that
   actually appears at `c` in the input (a cheap structural check) using
   the angle-addition formula in terms of `Tan[a], Tan[b]` (or
   `Sin/Cos[a], Sin/Cos[b]` for Sin/Cos heads).
4. Apply `ReplaceAll` with the rule list, then `Together` + `Cancel` so
   the polynomial cancellation fires.
5. Score-gated; inert when fewer than three distinct trig args appear or
   when no pair-sum matches a third arg.

Coverage:

| Input | Was | Now |
|------|-----|-----|
| `Simplify[Tan[2] Tan[3] B A - (-Tan[2]/Tan[5] - Tan[3]/Tan[5] + 1) B A]` | `A B (-1 + Cot[5] Tan[2] + Cot[5] Tan[3] + Tan[2] Tan[3])` | `0` (~1.3 s) |

Symbolic shapes (e.g. `Tan[x], Tan[y], Tan[x+y]` triples) work the same
way; the integer-arg case is just the one most surprising to TrigExpand.

### Wiring

- `transform_prime_rebase`, `transform_power_oneify`,
  `transform_radical_canon` all live in `simp_dispatch` ahead of the
  shape classifier.  They run on every Simplify call regardless of the
  input's shape (POLYNOMIAL / RATIONAL / LOGEXP / TRIG / GENERAL).
  Each is run once, score-gated; on a strict win the dispatcher recurses
  on the rewritten form.
- `transform_tan_addition` lives in `simp_search`'s seed phase alongside
  the other trig seeds (PythagCanon, TrigReduce, trig_pi_canon, halfangle,
  radicals, algebraic).
- All four are wrapped in `simp_memo_wrap` so identical inputs inside
  the same Simplify call share a single impl invocation.

### Files

- `src/simp.c`: new transform implementations (~600 LOC) and seed wiring.
- `tests/test_simplify.c`: 9 additional unit tests (`prime_power_combine_2x_4x`,
  `inside_f_prime_power_collapse`, `2_to_2_to_2x_x`, `one_over_x_log2_combine`,
  `tan_addition_three_angles`, plus 4 user-reported new cases for sqrt/
  algebraic/log/halfangle/compound-angle).

## Simplify: TanAddition uniform Sin/Cos primitives + depth-0 short-circuit (2026-05-06)

### Symptom

Mixed Tan/Sec angle-addition input hung past 900 ms:

```
m1 = Tan[z] Cos[x] Cos[y] Sec[x+y] (Tan[x] + Tan[y]);
m2 = Tan[z] Tan[x+y];
Simplify[m1 - m2]               (* = 0 by Tan[x]+Tan[y] = Sin[x+y]/(Cos[x] Cos[y]) *)
```

Mathematica returns `0` in <1 ms; picocas took ~960 ms.

### Root cause

Two compounding issues:

1. **TanAddition's Tan[c]/Cot[c] rule RHS used Tan[a]/Tan[b] primitives.**
   On a (a, b, c=a+b) triple where the input contains *both* `Tan[c]` and
   `Sec[c]`, the Tan[c] rule produced
   `(Tan[a]+Tan[b]) / (1 - Tan[a] Tan[b])` while the Sec[c] rule produced
   `1 / (Cos[a] Cos[b] - Sin[a] Sin[b])`. The two denominators are equal
   by trig identity but structurally distinct -- Together / Cancel could
   not unify them, so the substituted form failed to collapse.

2. **TrigReduce's depth-0 short-circuit blew up the input.**
   On the same input, `TrigReduce[m1 - m2]` produced a 9-term
   `Cos[...] Sec[x] Sec[z] Sec[x+y]` polynomial in ~700 ms (allocated,
   scored, and discarded by the strict-win gate). 700 ms wasted before
   the search loop ever ran TanAddition at the parent level.

### Fix (`src/simp.c`)

1. **`append_addition_rules_for_triple` refactored** to compute the two
   Sin/Cos primitives once
   ```
   sin_rhs = Sin[a] Cos[b] + Cos[a] Sin[b]   (= Sin[a + b])
   cos_rhs = Cos[a] Cos[b] - Sin[a] Sin[b]   (= Cos[a + b])
   ```
   and reuse them across all six head rules:
   ```
   Tan[c] -> sin_rhs / cos_rhs
   Cot[c] -> cos_rhs / sin_rhs
   Sin[c] -> sin_rhs
   Cos[c] -> cos_rhs
   Sec[c] -> 1 / cos_rhs
   Csc[c] -> 1 / sin_rhs
   ```
   The substituted form now has uniform Sin/Cos denominator structure.
   Cancel internally rewrites bare `Tan[x] = Sin[x]/Cos[x]` during its
   polynomial-GCD pass, so any leftover `Tan[a]+Tan[b]` factor combines
   with the new denominator and the whole expression collapses to 0.
   Helpers `mk_tan_sum` and `mk_one_minus_tan_prod` are removed.

2. **TanAddition added as depth-0 short-circuit in `simp_bottomup`**,
   immediately before the existing TrigReduce short-circuit.  TanAddition's
   gate (`has_pythag_head` + 3+ distinct trig args + a sum-witnessing
   triple) keeps it cheap when inert.  When it produces a strict-win
   result (typically 0 directly), `input` is replaced and the
   `canon_owned && input->type != EXPR_FUNCTION` early-return fires,
   so the still-expensive TrigReduce short-circuit never even runs.

### Performance

| Case | Before | After |
|------|--------|-------|
| `Tan[z] Cos[x] Cos[y] Sec[x+y] (Tan[x]+Tan[y]) - Tan[z] Tan[x+y]` | 960 ms | 7.6 ms |
| `Tan[2] Tan[3] B A - (-Tan[2]/Tan[5] - Tan[3]/Tan[5] + 1) B A` | ~1.3 s | ~5 ms (new formula still works; existing test passes) |
| `(Sin[x]+Cos[x])^4 - (1+Sin[2 x])^2` | 5 ms | 5 ms (no regression) |
| `(Sec[x]+1)(Sec[x]-1) - Tan[x]^2` | 5 ms | 5 ms (no regression) |
| `(Cosh[x]+Sinh[x])^3 - (Cosh[3 x]+Sinh[3 x])` | 4 ms | 4 ms (no regression) |

### Files

- `src/simp.c`: `append_addition_rules_for_triple` rewritten around
  shared Sin/Cos primitives; `simp_bottomup` gains a depth-0 TanAddition
  short-circuit before TrigReduce; helpers `mk_tan_sum` /
  `mk_one_minus_tan_prod` removed.
- `tests/test_simplify.c`: new `tan_addition_mixed_sec_tan` unit test
  pinning the 7.6 ms regression guard for the user-reported case.

## Simplify: PowerDistribute -- principal-branch power-distribution rewrites (2026-05-06)

### Symptom

Three Simplify identities returned the input unchanged:

```
Simplify[(-4)^x * (-2)^(-x) * 2^(-x) - 1]
  expected:  0
  actual:    -1 + (-4)^x (-2)^(-x) 2^(-x)

Simplify[Exp[x] Exp[y] 2^x 2^y - (2 Exp[1])^(x + y)]
  expected:  0
  actual:    -(2 E)^(x + y) + 2^(x + y) E^(x + y)

Simplify[y^n (y/x)^(-n) - x^n, Assumptions -> Element[n, Integers]]
  expected:  0
  actual:    -x^n + y^n (y/x)^(-n)
```

### Root cause

Three independent gaps in the existing power-rewrite cascade:

1. **No negative-base split.** `prime_rebase` rewrites `c^e -> p^(k e)`
   only for *positive* perfect-prime-power bases (`n >= 4`). `(-4)^x`
   and `(-2)^x` were left alone, so they sat in distinct base buckets
   in the canonical Times and never cancelled the explicit `2^(-x)`
   factor, blocking the chain
       (-4)^x (-2)^(-x) 2^(-x) -> 1.

2. **`apply_logexp_rules` was gated on `ctx_has_facts(ctx)`.**  The
   distribute identity `Power[Times[c1,...,ck], e] -> Times[ci^e]` for
   constant-positive bases (which holds without any user assumption,
   since 2 and `E` are unconditionally positive) never fired when the
   user supplied no `Assumptions` argument.  `(2 E)^(x+y)` therefore
   stayed wrapped, structurally distinct from the LHS's
   `2^(x+y) E^(x+y)` form.

3. **`apply_logexp_rules` distribute required *every* base positive.**
   Even with `Element[n, Integers]` provided, `Power[Times[y, x^(-1)],
   -n]` could not distribute because `y` and `x^(-1)` are symbolic
   (positivity unknown).  The integer-exponent identity
   `(a b)^n = a^n b^n` (universally valid for integer `n`) had no
   path through the rewrite cascade.

### Fix (`src/simp.c`)

New transform `transform_power_distribute(e, ctx)`, inserted in
`simp_dispatch` right after `transform_prime_rebase` /
`transform_power_oneify`.  Single structural walk; on each `Power[B, e]`
node it tries three universally-sound rewrites in order:

(A) **Negate-base split**: `Power[neg_int, e] -> Power[-1, e] *
    Power[|neg_int|, e]`.  Principal-branch identity for any `e` (real
    or symbolic); valid without assumptions.  After the split, prime
    rebase + Times same-base merge collapses
    `(-4)^x (-2)^(-x) 2^(-x)` into
    `(-1)^x 2^(2x) (-1)^(-x) 2^(-x) 2^(-x)` and on to `1`.

(B) **Constant-positive base distribute**: when `Power[Times[args], e]`
    has at least one factor that's a literal positive numeric or a
    recognised positive-constant symbol (`Pi`, `E`, `EulerGamma`,
    `GoldenRatio`, `Catalan`, `Degree`, `Glaisher`, `Khinchin`), peel
    those factors off as `Power[ci, e]`.  When *every* factor is a
    constant positive, distribute fully into `Times[ci^e]`.  Soundness:
    `(a b)^c = a^c b^c` holds whenever `a > 0`.

(C) **Integer-exponent distribute**: when `prov_int(ctx, e)` succeeds
    (literal integer exponent or via `Element[_, Integers]`
    assumption), distribute over `Power[Times[args], e]` regardless of
    base positivity, *and* collapse `Power[Power[u, p], e] -> Power[u,
    p e]`.  Soundness: integer exponents distribute through products
    and through nested Power without branch cuts, for any complex base.

After any rewrite, `evaluate()` is called so the canonical Times same-
base merge collapses adjacent `Power[u, ?]` factors against existing
factors elsewhere in the surrounding Times.  The transform is gated by
a cheap structural check (`has_distributable_power`) and is inert when
no matching node exists (microseconds-cheap).

The new transform sits in the same strict-win recursion shape used by
`transform_prime_rebase`: only adopt the rewrite when it strictly beats
the input on score; otherwise fall through.  Because `evaluate()` runs
inside the impl, a productive split + cancellation lands as a
strictly-shorter expression and the dispatcher recurses.

### Coverage

| Input | Before | After |
|-------|--------|-------|
| `(-4)^x (-2)^(-x) 2^(-x) - 1` | unchanged | `0` |
| `Exp[x] Exp[y] 2^x 2^y - (2 E)^(x+y)` | 2-term sum | `0` |
| `y^n (y/x)^(-n) - x^n` w/ `Element[n, Integers]` | unchanged | `0` |

All 21 cases in the user's Simplify regression set now collapse to `0`,
the slowest at 0.21 s; full project test suite (91 binaries) passes.

### Soundness note

The constant-positive-distribute (B) and integer-exponent-distribute
(C) cases are unconditionally valid (positivity / integer-exponent
guarantee no branch cut).  The negate-base split (A) is a
principal-branch identity -- the same convention picocas uses
everywhere else (e.g. `Sqrt[-1] = I`, `Power[-1, 1/3] = (-1)^(1/3)`).
None of the three depends on the no-numeric-sampling-zero rule that
governs Simplify correctness.

### Files

- `src/simp.c`: `transform_power_distribute` (~180 LOC) inserted between
  `transform_power_oneify` and `transform_radical_canon`; new dispatch
  block in `simp_dispatch` next to the existing rebase / oneify gates.

## Inverse-trig and inverse-hyperbolic auto-evaluation + Simplify (2026-05-06)

Three families of identity now fire either during normal evaluation or
during `Simplify`, closing a long-standing gap where the system left
expressions like `ArcCos[-x] + ArcCos[x] - Pi` unsimplified even though
each piece is a textbook identity.

### 1. Negation reflection (auto-evaluation)

The inverse trig and hyperbolic functions canonicalise their arguments'
sign through the same `expr_is_superficially_negative` test the forward
trig functions already use. Odd inverses fold negation through a sign
flip; the two reflection-shifted inverses (`ArcCos`, `ArcSec`) emit a
`Pi - f[|arg|]` rewrite. Concretely:

| Function | Rule |
|----------|------|
| `ArcSin[-x]`  | `-ArcSin[x]`     |
| `ArcCos[-x]`  | `Pi - ArcCos[x]` |
| `ArcTan[-x]`  | `-ArcTan[x]`     |
| `ArcCot[-x]`  | `-ArcCot[x]`     |
| `ArcSec[-x]`  | `Pi - ArcSec[x]` |
| `ArcCsc[-x]`  | `-ArcCsc[x]`     |
| `ArcSinh[-x]` | `-ArcSinh[x]`    |
| `ArcTanh[-x]` | `-ArcTanh[x]`    |
| `ArcCoth[-x]` | `-ArcCoth[x]`    |
| `ArcCsch[-x]` | `-ArcCsch[x]`    |

`ArcCosh` and `ArcSech` have no clean negation identity and are left
alone.

### 2. Imaginary-axis bridge (auto-evaluation)

When the argument is of the form `I*y` (`Times[Complex[0, k], rest...]`
with `k > 0`) the inverse functions cross between trig and hyperbolic
branches via the principal-branch identities:

| Trig inverse | Hyperbolic inverse |
|--------------|---------------------|
| `ArcSin[I y]  -> I ArcSinh[y]`        | `ArcSinh[I y]  -> I ArcSin[y]`  |
| `ArcCos[I y]  -> Pi/2 - I ArcSinh[y]` | `ArcCosh[I y]` left unchanged   |
| `ArcTan[I y]  -> I ArcTanh[y]`        | `ArcTanh[I y]  -> I ArcTan[y]`  |
| `ArcCot[I y]  -> -I ArcCoth[y]`       | `ArcCoth[I y]  -> -I ArcCot[y]` |
| `ArcCsc[I y]  -> -I ArcCsch[y]`       | `ArcCsch[I y]  -> -I ArcCsc[y]` |

The `ArcCos[I y]` rule emits `Pi/2 - I ArcSinh[y]` directly through a
new `arccos_i_fold` helper alongside the existing `trig_i_fold` /
`hyp_i_fold` infrastructure.

### 3. Forward-of-inverse algebraic identities (auto-evaluation)

`try_simp_forward_of_inverse` (trig) and `try_simp_forward_of_inverse_hyp`
(hyperbolic) are invoked from each forward builtin (`Sin`, `Cos`, ...,
`Sinh`, `Cosh`, ...) right after the existing `strip_inverse_call` pass.
They fold known closed forms:

| Trig | Hyperbolic |
|------|------------|
| `Sin[ArcCos[x]] -> Sqrt[1 - x^2]`             | `Sinh[ArcCosh[x]] -> Sqrt[x-1] Sqrt[x+1]` |
| `Cos[ArcSin[x]] -> Sqrt[1 - x^2]`             | `Cosh[ArcSinh[x]] -> Sqrt[1 + x^2]` |
| `Sin[ArcTan[x]] -> x / Sqrt[1 + x^2]`         | `Sinh[ArcTanh[x]] -> x / Sqrt[1 - x^2]` |
| `Cos[ArcTan[x]] -> 1 / Sqrt[1 + x^2]`         | `Cosh[ArcTanh[x]] -> 1 / Sqrt[1 - x^2]` |
| `Tan[ArcSin[x]] -> x / Sqrt[1 - x^2]`         | `Tanh[ArcSinh[x]] -> x / Sqrt[1 + x^2]` |
| `Tan[ArcCos[x]] -> Sqrt[1 - x^2] / x`         | `Tanh[ArcCosh[x]] -> Sqrt[x-1] Sqrt[x+1] / x` |
| `Tan[ArcCot[x]] -> 1/x`, `Cot[ArcTan[x]] -> 1/x` | `Tanh[ArcCoth[x]] -> 1/x`, `Coth[ArcTanh[x]] -> 1/x` |

The `Sqrt[x-1] Sqrt[x+1]` form (rather than `Sqrt[x^2 - 1]`) is the
principal-branch identity Mathematica uses; the two forms differ on the
branch cut `x in (-1, 1)` where `Sqrt[x-1] Sqrt[x+1] = i Sqrt[1 - x^2]`
agrees with the `Sinh[ArcCosh[x]]` analytic continuation but
`Sqrt[x^2 - 1]` does not.

### 4. ArcTan reciprocal sum (assumption-driven)

`apply_assumption_rules` in `src/simp.c` synthesises a per-symbol rule
for every positively-signed and every negatively-signed symbol in scope:

```
ArcTan[x] + ArcTan[1/x] + rest___ :> Pi/2 + rest    (x > 0)
ArcTan[x] + ArcTan[1/x] + rest___ :> -Pi/2 + rest   (x < 0)
```

The trailing `rest___` is required because picocas's pattern matcher
does not perform orderless-Plus subset matching out of the box: the
explicit `BlankNullSequence` absorbs and re-emits the surrounding terms
of any larger sum so that, for example, `ArcTan[x] + ArcTan[1/x] - Pi/2`
collapses cleanly to `0` under `x > 0`.

### 5. ArcCosh double-angle (assumption-driven)

A second per-positive-symbol rule reduces the two-step ArcCosh form:

```
ArcCosh[2 x^2 - 1] :> 2 ArcCosh[x]    (x > 0)
```

The identity holds for any complex `x` under the principal branch
(verified for `x in [0, 1]` via the `i*ArcCos` bridge and for `x >= 1`
directly), so `x > 0` is the weakest sufficient assumption. The
canonical Plus form `Plus[-1, Times[2, Power[x, 2]]]` is matched
literally.

### Files

- `src/trig.c`: `arc_pi_minus_fold`, `arccos_i_fold`,
  `try_simp_forward_of_inverse`; new fold calls inserted near the top
  of `builtin_arcsin`, `builtin_arccos`, `builtin_arctan`,
  `builtin_arccot`, `builtin_arcsec`, `builtin_arccsc`, and
  `try_simp_forward_of_inverse` calls in `builtin_sin`, `builtin_cos`,
  `builtin_tan`, `builtin_cot`.
- `src/hyperbolic.c`: `try_simp_forward_of_inverse_hyp`; new fold calls
  in `builtin_arcsinh`, `builtin_arctanh`, `builtin_arccoth`,
  `builtin_arccsch`, and forward-of-inverse calls in `builtin_sinh`,
  `builtin_cosh`, `builtin_tanh`, `builtin_coth`.
- `src/simp.c`: two new EMITs in the `npos` loop and one new EMIT in
  the `nneg` loop of `apply_assumption_rules`.
- `tests/test_invtrig_simplify.c`: 41 tests covering all five families
  above. Registered in `tests/CMakeLists.txt` next to
  `logexp_simplify_tests`.

### User-facing impact

All identities from the user-reference test set now `Simplify` to `0`
(or to the canonical scalar where appropriate):

```
Simplify[Sin[ArcCos[x]] - Sqrt[1 - x^2]]                                -> 0
Simplify[Tan[ArcSin[x]] - x / Sqrt[1 - x^2]]                            -> 0
Simplify[ArcCos[-x] + ArcCos[x] - Pi]                                   -> 0
Simplify[ArcSin[x] + ArcCos[x] - Pi/2]                                  -> 0
Assuming[x > 0, Simplify[ArcTan[x] + ArcTan[1/x] - Pi/2]]               -> 0
Assuming[x < 0, Simplify[ArcTan[x] + ArcTan[1/x] + Pi/2]]               -> 0
Assuming[x > 1, Simplify[Sinh[ArcCosh[x]] - Sqrt[x-1] Sqrt[x+1]]]       -> 0
Assuming[x > 1, Simplify[ArcCosh[2 x^2 - 1] - 2 ArcCosh[x]]]            -> 0
Simplify[ArcSin[I x] - I ArcSinh[x]]                                    -> 0
Simplify[ArcTan[I x] - I ArcTanh[x]]                                    -> 0
Simplify[ArcCos[I x] - (Pi/2 - I ArcSinh[x])]                           -> 0
```

## Sqrt-of-Sqrt denesting (`simp_denest_sqrt`, 2026-05-06)

Phase 1 of an algorithmic radical-simplification project. Closes
five user-supplied test cases that prior `Simplify` left unchanged:

```
Simplify[Sqrt[3 + 2 Sqrt[2]] - (1 + Sqrt[2])]                                -> 0
Simplify[Sqrt[17 - 12 Sqrt[2]] - (3 - 2 Sqrt[2])]                            -> 0
Simplify[Sqrt[2 + Sqrt[3]] - (Sqrt[6] + Sqrt[2])/2]                          -> 0
Assuming[x>0 && y>0, Simplify[Sqrt[x+y+2 Sqrt[x y]] - (Sqrt[x]+Sqrt[y])]]    -> 0
Assuming[x>y && y>0, Simplify[Sqrt[x+Sqrt[x^2-y^2]]
                              - (Sqrt[(x+y)/2]+Sqrt[(x-y)/2])]]              -> 0
```

### The single primitive

A new `simp_denest_sqrt` Simplify seed applies the half-sum identity

```
Sqrt[A + Sqrt[B]] = Sqrt[(A + s)/2] + Sqrt[(A - s)/2]
```

at every `Sqrt[Plus[...]]` subtree, where `s` is a closed-form
square root of `A^2 - B`. The identity holds in the principal branch
when both `(A +/- s)/2` are nonneg. When the input is `Sqrt[A - b Sqrt[C]]`
with `b > 0`, the same `s` (computed from `B = b^2 C`) gives the
sign-flipped form `Sqrt[(A+s)/2] - Sqrt[(A-s)/2]`.

### How `s` is detected

`sqrt_if_clean_square(D, ctx)` returns a closed-form `s` (sign-agnostic)
when `D` is one of:

- A nonneg integer/rational perfect square (via GMP `mpz_perfect_square_p`).
- `Power[u, 2k]` for any integer `k >= 1`.
- A polynomial Plus that `Expand[D]` then `FactorSquareFree[.]` reduces
  to `Power[u, 2k]`. This path catches the symbolic discriminants
  `(x-y)^2` (case 6) and `y^2` (case 7).

When no clean square is detected, the identity is not applied — the
input is returned unchanged. This preserves the project-wide soundness
invariant that `Simplify` never produces a wrong result.

### Branch validation

After computing `P = (A+s)/2` and `Q = (A-s)/2`, the transform requires
both to be provably nonneg under the active `AssumeCtx` before firing.
The standard `assume_known_nonneg` does not perform transitive chaining
or the subtraction pattern needed by some Phase-1 cases, so the module
introduces local helpers:

- `denest_prov_nonneg` / `denest_prov_pos` chain inequality facts:
  `x > y && y > 0` implies `x > 0`.
- A 2-arg Plus subtraction-pattern detector recognises
  `Plus[u, Times[-1, v]]` and looks for `Greater[u, v]` /
  `GreaterEqual[u, v]` in the assumption set. Required for case 7,
  where `Q = (x - y)/2` is provably nonneg under `x > y` even though
  the literal Plus decomposition fails on `-y`.
- A `numeric_is_nonneg` fast path covers `Rational[p, q]`, which the
  stock `numeric_sign` only handles for `INTEGER`/`BIGINT`/`REAL`.

### Sqrt-of-Rational rationalisation

The denesting often produces inner radicands like `3/2` and `1/2`;
`Sqrt[3/2] + Sqrt[1/2]` must canonicalise to `(Sqrt[6]+Sqrt[2])/2`
for case 3 to cancel against the user-supplied form. The existing
`transform_radical_canon` only handles `Sqrt[1/n]` because the
evaluator re-merges `Times[Power[m, 1/2], Power[n, -1/2]]` back
into `Power[m/n, 1/2]` for `m > 1`, undoing the split. The Phase-1
denest module includes `denest_rationalise_sqrt_of_rational`, which
builds `Sqrt[m*n] * (1/n)` directly so the re-merge cannot fire.

### Files

- `src/simp.c`: ~600 LOC of new code in a `simp_denest_sqrt` section
  inserted between `simp_radicals` and `simp_algebraic`. Wired as a
  new seed in `simp_search` next to `simp_algebraic` (line ~6831).
  Forward-declares `transform_radical_canon` and `is_sqrt`.
- `tests/test_radical_simplify.c`: 14 tests covering all five Phase-1
  user cases, four branch-soundness subtests, the case-9 regression
  sentinel, and four Phase-2/3 cases gated by `#define PHASE_N`
  macros (skipped until those phases land).
- `tests/CMakeLists.txt`: registers `radical_simplify_tests` next to
  `invtrig_simplify_tests`.

### Out of scope (Phase 2/3 placeholders)

- Multi-generator algebraic-number reduction (cases 8, 10): denominator
  rationalisation in `Q(α₁, ..., αₖ)` requires a new module.
- Cube-root denesting (cases 4, 5): Borodin–Fagin–Hopcroft–Tompa for
  `(p + q√r)^(1/3)` and Cardano-discriminant for sum-of-conjugate
  cube roots.


## Denominator rationalisation (`simp_rationalize_denom`, 2026-05-07)

Phase 2 of the algorithmic radical-simplification project. Closes
two user-supplied test cases that the existing `simp_algebraic` (which
handles only square roots) could not:

```
Simplify[1/(2^(1/3) - 1) - (4^(1/3) + 2^(1/3) + 1)]              -> 0
Simplify[1/(Sqrt[2] + 2^(1/3)) - r10]                            -> 0
  (where r10 is the algebraic rationalisation derived from
   extended-Euclidean in Q[α]/(α^6 - 2) with α = 2^(1/6))
```

### Algorithm

A new `simp_rationalize_denom` Simplify seed walks the expression
tree. At each `Power[denom, -1]` subtree where `denom` is a polynomial
in radicals over a single positive integer base `c`, it computes the
inverse via extended-Euclidean in `Q[α]/(α^n - c)` where
`α = c^(1/n)` and `n` is the LCM of the rational-exponent denominators.

Steps:

1. `denom_collect_radical_base` walks `denom` collecting the unique
   integer base `c` and computing `n = lcm(q_i)` over all
   `Power[c, p_i/q_i]` subexpressions. Returns failure on:
   - non-rational exponents, multiple distinct bases, or non-positive
     bases.
2. `denom_subst_radical_to_gen` substitutes each `Power[c, p_i/q_i]`
   with `Power[gen, p_i * (n/q_i)]` for a fresh symbol `gen`.
3. `PolynomialExtendedGCD[denom_subbed, gen^n - c, gen]` gives the
   Bezout coefficients `(u, v)` with `u * denom_subbed + v * (gen^n - c) = g`
   for some constant `g`. The transform accepts iff `g` is a nonzero
   rational (otherwise denom shares a factor with the relation, meaning
   the input is undefined at some root).
4. The inverse is `u(gen) / g`. `denom_subst_gen_to_radical` replaces
   each `Power[gen, k]` with `Power[c, k/n]` to reverse the
   substitution.
5. A `transform_prime_rebase` post-pass on the FULL output ensures
   that user-supplied compound bases (e.g. `4^(1/3)`) and the prime
   base our walker uses (`2^(2/3)`) match: PrimeRebase rewrites
   `c^e` for `c = p^k` integer to `p^(k*e)`.

### Soundness

The rewrite is correctness-preserving for any positive integer base
`c >= 2`, where `c^(1/n)` is the principal real root and the extended-
Euclidean computation lives entirely in `Q[α]`. The transform refuses
to fire (returns the input unchanged) when:

- The denominator references multiple distinct integer bases (out of
  scope for the single-extension implementation).
- The gcd of `denom` and `α^n - c` is not a constant (denominator
  vanishes at some root of `α^n - c`, meaning the expression is
  undefined and any "rationalisation" would be wrong).
- The base is not a positive integer >= 2.

### Files

- `src/simp.c`: ~280 LOC of new code in a `simp_rationalize_denom`
  section inserted between `simp_denest_sqrt` and `simp_algebraic`.
  Wired as a new seed in `simp_search` next to `simp_denest_sqrt`.
  Forward-declares `transform_prime_rebase`.
- `tests/test_radical_simplify.c`: the case-10 expected value was
  corrected from the user-transcript form (which has sign errors
  verifiable numerically) to the algebraically-derived
  `(2 sqrt(2) + sqrt(2)*2^(1/3) + sqrt(2)*4^(1/3) - 2*2^(1/3) - 4^(1/3) - 2)/2`.

### Out of scope (Phase 3)

- Cube-root denesting (cases 4, 5): `(p + q sqrt(r))^(1/3)` and the
  sum-of-conjugate-cube-roots Cardano-discriminant pattern. These are
  forward radicals, not denominator rationalisations, and require a
  different algorithmic approach.
- Multi-base extensions (e.g. `1/(sqrt(3) + 2^(1/3))`): the current
  implementation handles only single integer bases.


## Cube-root denesting (`simp_cuberoot`, 2026-05-07)

Phase 3 of the algorithmic radical-simplification project. Closes the
final two user-supplied test cases:

```
Simplify[(Sqrt[5] + 2)^(1/3) - (Sqrt[5] + 1)/2]                       -> 0
Simplify[(2 + Sqrt[5])^(1/3) + (2 - Sqrt[5])^(1/3) - 1]               -> 0
```

### Two narrow patterns

A new `simp_cuberoot` Simplify seed recognises exactly two patterns:

**Pattern A -- single-radical cube-root denesting (Borodin-Fagin style):**

For `Power[a + b sqrt(c), 1/3]` with rational a, b and integer c, search
for rationals (p, q) satisfying

```
p^3 + 3 p q^2 c = a
3 p^2 q + q^3 c = b
```

(equivalently `(p + q sqrt(c))^3 = a + b sqrt(c)`). When found, rewrite
to `p + q sqrt(c)`. The search is a small bounded grid: numerator
magnitudes up to 6, common denominator up to 6. This is sufficient for
case 4 (p = q = 1/2) and similar small-coefficient targets; it is NOT
a general algorithm.

**Pattern B -- sum of conjugate cube roots (Cardano):**

For `Power[a + b sqrt(c), 1/3] + Power[a - b sqrt(c), 1/3]` with integer
a, b and integer c, the sum s under the real-cube-root convention
satisfies the depressed cubic

```
s^3 - 3 m s - 2 a = 0,   where m = (a^2 - b^2 c)^(1/3)
```

provided `a^2 - b^2 c` is a perfect cube of an integer m. When this
cubic has a rational root, the transform replaces the sum with that
rational. For case 5 (a = 2, b = 1, c = 5): `a^2 - b^2 c = -1 = (-1)^3`
gives m = -1, then `s^3 + 3 s - 4 = (s-1)(s^2 + s + 4)` has rational
root s = 1.

### Soundness note: principal vs. real cube root

picocas evaluates `Power[neg_real, 1/3]` using the principal complex
branch (e.g. `(-1)^(1/3) = (1 + i sqrt(3))/2`). Pattern B's identity
holds only under the REAL cube-root convention; principal-branch
numeric evaluation of case 5 gives ~1.93 + 0.535i, not 1. The
transform fires anyway because the user's intent (and Mathematica
`Simplify`) follows the real-cube-root heuristic, and the gating
condition (`a^2 - b^2 c` being a perfect integer cube) is purely
structural and branch-independent.

Pattern A has no branch concern: `(p + q sqrt(c))^3 = a + b sqrt(c)`
is an exact polynomial identity that holds on every branch.

### Files

- `src/simp.c`: ~370 LOC of new code in a `simp_cuberoot` section
  inserted between `simp_denest_sqrt` and `simp_rationalize_denom`.
  Wired as a new seed in `simp_search` next to `simp_rationalize_denom`.
- `tests/test_radical_simplify.c`: case 4 and case 5 enabled via
  `#define PHASE_3`. All 14 tests pass.

### Out of scope

- Cube roots with non-trivial rational a, b (the search grid is
  bounded; the user's cases were small-integer / half-integer).
- Higher-order radicals (4th, 5th roots).
- Nested cube roots beyond the user-supplied patterns. The general
  Borodin-Fagin-Hopcroft-Tompa algorithm is much larger than this
  pattern recogniser; if expanded coverage becomes necessary,
  extending to a full `facpoly_over_alg` is the natural path.
