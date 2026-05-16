# Functional Programming in Mathilda

In an [earlier note](./Mathilda_patterns_blog.md) I walked through pattern matching in Mathilda. This time I would like to tour the functional-programming side of the language: the handful of combinators that, together with pattern rules, do most of the heavy lifting in Mathematica-flavoured code.

Mathilda aims to be a faithful recreation of Mathematica's core programming language in a very small C99 codebase [1]. Functional programming is part of that core: `Map`, `Apply`, `Fold`, `Nest`, `FixedPoint`, `Select`, and friends let you express computations declaratively, without explicit loops or mutable state. Every example below was run in Mathilda; the output is the actual `Out[]` line the interpreter produced.

## Map and Apply — the two workhorses

`Map[f, list]` applies `f` to each element of a list. The `/@` operator is shorthand for the same thing:

```mathematica
In[1]:= Map[f, {a, b, c}]
Out[1]= {f[a], f[b], f[c]}

In[2]:= f /@ {a, b, c}
Out[2]= {f[a], f[b], f[c]}
```

In combination with a pure function (`#` for the argument, `&` to close the function) we can square a list of numbers without ever writing a loop:

```mathematica
In[3]:= Map[#^2 &, Range[5]]
Out[3]= {1, 4, 9, 16, 25}
```

`Apply[f, expr]` — or `f @@ expr` — replaces the head of `expr` with `f`. It is the natural complement to `Map`. Where `Map` distributes a function over the arguments, `Apply` folds the arguments into a single call:

```mathematica
In[4]:= Plus @@ Range[10]
Out[4]= 55

In[5]:= Times @@ Range[1, 10]
Out[5]= 3628800
```

Because every Mathilda expression is a tagged tree, `Apply` can convert between heads:

```mathematica
In[6]:= Apply[List, a + b + c]
Out[6]= {a, b, c}

In[7]:= Apply[Plus, {a, b, c, d}]
Out[7]= a + b + c + d
```

`Map` and `Apply` compose cleanly. The sum of the first 100 squares is a `Map` followed by an `Apply`:

```mathematica
In[8]:= Apply[Plus, Map[#^2 &, Range[100]]]
Out[8]= 338350
```

Both take an optional level specification. `Apply[Plus, list, {1}]` sums at level 1 only, i.e. it turns an outer list of inner lists into an outer list of row-sums:

```mathematica
In[9]:= Apply[Plus, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, {1}]
Out[9]= {6, 15, 24}

In[10]:= Map[f, {{a, b}, {c, d}}, {2}]
Out[10]= {{f[a], f[b]}, {f[c], f[d]}}
```

## Select — filtering

`Select[list, pred]` keeps only the elements for which `pred[element]` is `True`. With `PrimeQ` we get the primes below 20 in a single line:

```mathematica
In[11]:= Select[Range[20], PrimeQ]
Out[11]= {2, 3, 5, 7, 11, 13, 17, 19}
```

Any predicate — built-in or pure-function — works. Here are the multiples of 7 under 100:

```mathematica
In[12]:= Select[Range[100], Mod[#, 7] == 0 &]
Out[12]= {7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 98}
```

And twin primes under 50, using a predicate that combines two `PrimeQ` checks:

```mathematica
In[13]:= Select[Range[2, 50], PrimeQ[#] && PrimeQ[# + 2] &]
Out[13]= {3, 5, 11, 17, 29, 41}
```

`Select[list, pred, n]` stops after the first `n` matches:

```mathematica
In[14]:= Select[Range[1000], PrimeQ, 3]
Out[14]= {2, 3, 5}
```

## Nest and NestList — iterated application

`Nest[f, x, n]` applies `f` to `x` a total of `n` times; `NestList` returns the full history including the seed:

```mathematica
In[15]:= Nest[Sqrt, x, 3]
Out[15]= Sqrt[Sqrt[Sqrt[x]]]

In[16]:= NestList[f, x, 4]
Out[16]= {x, f[x], f[f[x]], f[f[f[x]]], f[f[f[f[x]]]]}
```

Newton iteration for $\sqrt{2}$ is a pure function iterated a few times:

```mathematica
In[17]:= NestList[(# + 2/#)/2 &, 1.0, 6]
Out[17]= {1.0, 1.5, 1.41667, 1.41422, 1.41421, 1.41421, 1.41421}
```

Pascal's triangle falls out of `NestList` with a very small kernel — repeatedly multiply by `(1 + x)` and expand. The rows of Pascal's triangle are the coefficients:

```mathematica
In[18]:= NestList[Expand[(1 + x) #] &, 1, 5]
Out[18]= {1, 1 + x, 1 + 2 x + x^2, 1 + 3 x + 3 x^2 + x^3, 1 + 4 x + 6 x^2 + 4 x^3 + x^4, 1 + 5 x + 10 x^2 + 10 x^3 + 5 x^4 + x^5}
```

The Fibonacci numbers have a classic functional encoding: iterate the map `{a, b} -> {b, a + b}` starting from `{0, 1}`:

```mathematica
In[19]:= NestList[{#[[2]], #[[1]] + #[[2]]} &, {0, 1}, 10]
Out[19]= {{0, 1}, {1, 1}, {1, 2}, {2, 3}, {3, 5}, {5, 8}, {8, 13}, {13, 21}, {21, 34}, {34, 55}, {55, 89}}

In[20]:= Map[First, %]
Out[20]= {0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55}
```

## Fold and FoldList — reducing with an accumulator

`Fold[f, x, list]` threads an accumulator through a list, applying `f[acc, elt]` for each element. `FoldList` returns the running trail of accumulators. Sums and products are one-liners:

```mathematica
In[21]:= Fold[Plus, 0, Range[1000]]
Out[21]= 500500

In[22]:= Fold[Times, 1, Range[10]]
Out[22]= 3628800
```

Cumulative sums are a `FoldList`:

```mathematica
In[23]:= FoldList[Plus, 0, Range[5]]
Out[23]= {0, 1, 3, 6, 10, 15}
```

`Fold` also handles reductions that aren't `Plus` or `Times`. For example, building a number out of its digits is `10 acc + digit`:

```mathematica
In[24]:= Fold[10 #1 + #2 &, 0, {4, 5, 1, 6, 7, 8}]
Out[24]= 451678
```

Horner's form for evaluating polynomials is the same pattern with `x` instead of `10`. Taking coefficients `{1, -3, 2, -1, 5}` (highest power first) we get $x^4 - 3x^3 + 2x^2 - x + 5$:

```mathematica
In[25]:= Fold[x #1 + #2 &, 0, {1, -3, 2, -1, 5}]
Out[25]= 5 + x (-1 + x (2 + x (-3 + x)))
```

Matrix powers are elegant too: fold `Dot` over a list of $n$ copies of the matrix starting from the identity. Here we compute $\left(\begin{smallmatrix}1&1\\0&1\end{smallmatrix}\right)^5 = \left(\begin{smallmatrix}1&5\\0&1\end{smallmatrix}\right)$:

```mathematica
In[26]:= Fold[#1 . #2 &, IdentityMatrix[2], Table[{{1, 1}, {0, 1}}, {5}]]
Out[26]= {{1, 5}, {0, 1}}
```

Finally, truncated harmonic partial sums $1 + 1/2 + 1/3 + \dots$ are a `FoldList` over reciprocals:

```mathematica
In[27]:= FoldList[#1 + 1/#2 &, 0, Range[6]]
Out[27]= {0, 1, 3/2, 11/6, 25/12, 137/60, 49/20}
```

## FixedPoint and FixedPointList — iterating to convergence

`FixedPoint[f, x]` keeps applying `f` until the result stops changing (or, optionally, until a user-supplied `SameTest -> s` predicate says so). `FixedPointList` returns the entire history.

Newton's method for $\sqrt{2}$ now needs no hand-chosen iteration count:

```mathematica
In[28]:= FixedPoint[(# + 2/#)/2 &, 1.0]
Out[28]= 1.41421
```

The golden ratio emerges from the continued fraction $x \mapsto 1 + 1/x$:

```mathematica
In[29]:= FixedPointList[1 + 1/# &, 1.0]
Out[29]= {1.0, 2.0, 1.5, 1.66667, 1.6, 1.625, 1.61538, 1.61905, 1.61765, 1.61818, 1.61798, 1.61806, 1.61803, 1.61804, 1.61803, 1.61803, ...}
```

The same trick works with $x \mapsto \sqrt{x + 1}$:

```mathematica
In[30]:= FixedPointList[Sqrt[# + 1] &, 1.0]
Out[30]= {1.0, 1.41421, 1.55377, 1.59805, 1.61185, 1.61612, 1.61744, 1.61785, 1.61798, 1.61802, 1.61803, 1.61803, ...}
```

`FixedPoint` also works symbolically. Repeatedly expanding `(1 + x)^3` stabilises to the fully expanded form:

```mathematica
In[31]:= FixedPoint[Expand[#] &, (1 + x)^3]
Out[31]= 1 + 3 x + 3 x^2 + x^3
```

And a rule-based Euclidean GCD: keep applying `{a, b} -> {b, Mod[a, b]}` until the second element is zero:

```mathematica
In[32]:= FixedPoint[# /. {a_, b_} /; b != 0 -> {b, Mod[a, b]} &, {28, 21}]
Out[32]= {7, 0}
```

## NestWhile and NestWhileList — conditional iteration

`NestWhile[f, x, test]` iterates `f` starting at `x` as long as `test` yields `True` on the most recent result. `NestWhileList` records the whole trajectory. Halving until odd is a one-liner:

```mathematica
In[33]:= NestWhile[#/2 &, 1024, EvenQ]
Out[33]= 1

In[34]:= NestWhileList[#/2 &, 2^20, EvenQ]
Out[34]= {1048576, 524288, 262144, 131072, 65536, 32768, 16384, 8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1}
```

A more famous example is the Collatz conjecture: starting at `n`, halve if even, otherwise `3 n + 1`, and keep going while the value is above one. The sequence from 27 is notoriously long:

```mathematica
In[35]:= collatz = NestWhileList[If[EvenQ[#], #/2, 3 # + 1] &, 27, # > 1 &];

In[36]:= Length[collatz]
Out[36]= 112

In[37]:= Max[collatz]
Out[37]= 9232
```

112 steps, reaching a maximum of 9232 before falling to 1 — the same trajectory we saw via pattern rules in the [earlier post](./Mathilda_patterns_blog.md).

## Combinatorial building blocks: Tuples, Permutations, Outer, Inner

`Tuples[list, n]` enumerates all length-`n` tuples drawn from `list` with replacement — the $n$-th Cartesian power:

```mathematica
In[38]:= Tuples[{H, T}, 3]
Out[38]= {{H, H, H}, {H, H, T}, {H, T, H}, {H, T, T}, {T, H, H}, {T, H, T}, {T, T, H}, {T, T, T}}
```

Counting the five-flip sequences with exactly three heads is `Tuples` composed with `Select` and `Count`:

```mathematica
In[39]:= Length[Select[Tuples[{H, T}, 5], Count[#, H] == 3 &]]
Out[39]= 10
```

`Tuples` with a single list-of-lists argument takes one element from each sub-list — the cross product of several alphabets:

```mathematica
In[40]:= Tuples[{{1, 2}, {x, y}, {A, B}}]
Out[40]= {{1, x, A}, {1, x, B}, {1, y, A}, {1, y, B}, {2, x, A}, {2, x, B}, {2, y, A}, {2, y, B}}
```

Combining `Tuples` and `Select` gives a brute-force search for Pythagorean triples with legs up to 20:

```mathematica
In[41]:= Select[Tuples[Range[20], 3], #[[1]]^2 + #[[2]]^2 == #[[3]]^2 && #[[1]] < #[[2]] &]
Out[41]= {{3, 4, 5}, {5, 12, 13}, {6, 8, 10}, {8, 15, 17}, {9, 12, 15}, {12, 16, 20}}
```

`Permutations` is the without-replacement cousin. There are $6! = 720$ permutations of `{1, 2, 3, 4, 5, 6}`:

```mathematica
In[42]:= Length[Permutations[Range[6]]]
Out[42]= 720
```

`Outer[f, a, b]` builds a generalized outer product: the matrix whose $(i, j)$ entry is `f[a[[i]], b[[j]]]`. With `Times` it's the ordinary outer product; with `List` it's the index grid:

```mathematica
In[43]:= Outer[Times, {1, 2, 3}, {x, y, z}]
Out[43]= {{x, y, z}, {2 x, 2 y, 2 z}, {3 x, 3 y, 3 z}}

In[44]:= Outer[List, {1, 2, 3}, {1, 2, 3}]
Out[44]= {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}, {{3, 1}, {3, 2}, {3, 3}}}
```

`Inner[f, u, v, g]` is the corresponding generalized inner product — the contraction of `Outer[f, u, v]` along the diagonal via `g`. The dot product is the special case `Inner[Times, u, v, Plus]`:

```mathematica
In[45]:= Inner[Times, {a, b, c}, {x, y, z}, Plus]
Out[45]= a x + b y + c z
```

## Through and Distribute — algebraic helpers

`Through` applies an expression that is a sum (or product) of functions to a common argument, distributing the argument inwards:

```mathematica
In[46]:= Through[(f + g + h)[x]]
Out[46]= f[x] + g[x] + h[x]

In[47]:= Through[(p + q)[x, y]]
Out[47]= p[x, y] + q[x, y]
```

`Distribute` is the symbolic counterpart — it expands `f[a + b, c + d]` into `f[a, c] + f[a, d] + f[b, c] + f[b, d]`:

```mathematica
In[48]:= Distribute[f[a + b, c + d]]
Out[48]= f[a, c] + f[a, d] + f[b, c] + f[b, d]
```

## Putting it all together

These combinators compose. Here is the sum of primes under 100, written three different ways — first with `Apply` + `Select`, then with `Fold` over a filtered range, then as a single-line `Apply` over `Select`:

```mathematica
In[49]:= Apply[Plus, Select[Range[100], PrimeQ]]
Out[49]= 1060

In[50]:= Apply[Times, Select[Range[20], PrimeQ]]
Out[50]= 9699690
```

And the number partition problem — how many ways can the values in `Range[3]^3` sum to 5? — is just `Tuples` + `Select`:

```mathematica
In[51]:= Select[Tuples[Range[3], 3], Total[#] == 5 &]
Out[51]= {{1, 1, 3}, {1, 2, 2}, {1, 3, 1}, {2, 1, 2}, {2, 2, 1}, {3, 1, 1}}
```

With just a small vocabulary — `Map`, `Apply`, `Select`, `Fold`, `Nest`, `FixedPoint`, `NestWhile`, `Tuples`, and `Outer` — we can express a surprisingly wide swath of computations in a declarative style, with no explicit loops and no mutable variables. Combined with pattern matching, this is most of what you need to program comfortably in the Mathematica style.

[1] https://en.wikipedia.org/wiki/Derive_(computer_algebra_system)
