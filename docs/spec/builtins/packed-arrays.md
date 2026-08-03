# Packed arrays

A **packed list** is an ordinary `List` that the system stores as a dense
machine-precision buffer instead of one expression node per element. It is
invisible: same `Head`, same printed form, same elements, same ordering, same
pattern matches. Only `NDArrayQ` tells you it is packed.

This is distinct from [`NDArray[...]`](linear-algebra.md), which is a *visibly*
different value with `Head` `NDArray`. Both use the same storage; they differ
only in what they claim to be.

```mathematica
In[1]:= a = ToNDArray[{1., 2., 3.}]
Out[1]= {1., 2., 3.}

In[2]:= {Head[a], ListQ[a], AtomQ[a], NDArrayQ[a]}
Out[2]= {List, True, False, True}

In[3]:= a === {1., 2., 3.}
Out[3]= True

In[4]:= b = NDArray[{1., 2., 3.}]
Out[4]= NDArray[{1., 2., 3.}]

In[5]:= {Head[b], ListQ[b], AtomQ[b], b === {1., 2., 3.}}
Out[5]= {NDArray, False, True, False}
```

## Why

A list of *n* elements is *n* separate expression nodes, and the evaluator
sweeps every argument of every call on every pass. Holding a large list is
therefore not free even when nothing is being computed with it:

| on a 10^6-element list of Reals | ordinary `List` | packed |
|---|---|---|
| `Length[x]` | 30.6 ms | 2 µs |
| `x[[7]]` | 32.3 ms | 82 µs |
| `Total[x]` | 95.8 ms | 4.0 ms |
| `Sin[x]` | 320 ms | 16.8 ms |

`Length` costs 30 ms because the argument sweep visits 10^6 nodes and rebuilds
the call, not because `Length` is doing anything — about 30 ns per element per
evaluation pass, charged every time the list is an argument to anything.

## The contract

Packing changes **representation and nothing else** — not a value, not an
element's head, not a printed form, not an ordering. Where the packed form
cannot express the exact answer, the operation abandons and the ordinary
list implementation runs instead: slower on that path, never different.

That rules out several conversions a purely numeric view would allow:

```mathematica
In[6]:= DataType[ToNDArray[{1, 2, 3}]]        (* exact integers stay exact *)
Out[6]= "int64"

In[7]:= Head[ToNDArray[{1, 2, 3}][[2]]]
Out[7]= Integer

In[8]:= NDArrayQ[ToNDArray[{1, 2.5}]]         (* mixed exact/inexact declines: *)
Out[8]= False                                 (* one buffer cannot hold both heads *)

In[9]:= ToNDArray[{9007199254740993}][[1]]    (* exact beyond 2^53 *)
Out[9]= 9007199254740993
```

## When packing happens

A list packs automatically when a **producer** builds one that is rectangular,
holds nothing but uniformly exact or uniformly inexact machine numbers, and has
at least **4 elements in total** (the product of its dimensions, so a
3 x 100000 matrix qualifies, and so does a 2 x 2). `ToNDArray` ignores the threshold.

| producer | how |
|---|---|
| `Range` (real and exact-integer) | writes the buffer directly |
| `ConstantArray` | writes the buffer directly |
| `RandomReal` | writes the buffer directly |
| `Table` (machine-real compiled body) | writes the buffer directly |
| `Subdivide` (a machine-real interval) | writes the buffer directly |
| `IdentityMatrix` | writes the buffer directly (`int64`) |
| `DiagonalMatrix` | writes the buffer directly (`int64` or `float64`) |
| `UnitVector` | writes the buffer directly (`int64`, or `float64` under `WorkingPrecision -> MachinePrecision`) |
| `Table` (all other branches) | offered after building |
| `Array`, `RandomInteger` | offered after building |
| `Sort`, `Select` | offered after building |
| `NestList`, `FoldList`, `NestWhileList`, `FixedPointList` | offered after building |

A producer writes a buffer only where every element it would otherwise have
boxed has the **same head**, which is what keeps the representation invisible.
For a producer that *invents* elements — the zeros in `DiagonalMatrix`, the
interior points of `Subdivide` — that is a question about the input's exactness,
and the answer is the ordinary numeric-tower one: **a machine `Real` anywhere
makes the whole result machine-real**, so `DiagonalMatrix[{1, 2, 3.}]` is
`{{1., 0., 0.}, {0., 2., 0.}, {0., 0., 3.}}` and `Subdivide[0, 1., 4]` is
`{0., 0.25, 0.5, 0.75, 1.}`. One dtype, so both pack.

Only *machine* `Real` is contagious, and the exceptions are what the rule is
for: an MPFR entry keeps the exact zeros (`DiagonalMatrix[{1.\`30, 2}]` has an
exact `0`), an all-exact input stays exact and packs as `int64`, and a symbolic
entry keeps the result unpacked while the zeros around it still turn `Real`
(`DiagonalMatrix[{a, 1.}]` is `{{a, 0.}, {0., 1.}}`).

The direct producers never build the expression nodes at all, which is where most
of the win is: packing `Range[1., 10^6]` after the fact costs 340 ms + 52 ms,
while writing 10^6 doubles costs under 1 ms.

The threshold is 4 — the element count of a 2 x 2 matrix, so the rule is *any
matrix packs*. It was 250 until 2026-08-02, chosen for blast radius rather than
for cost (the break-even is around two elements), and the margin turned out to
be expensive: a 6 x 6 is 36 elements, so it never packed, and the LAPACK path is
unreachable from an unpacked list. `Det` of a 6 x 6 cost 102.8 ms against
0.189 ms packed. Lowering it was swept against pattern matching, rule
application, `Table`, `Expand`, `Solve`, `D`, `Integrate`, `Sort` and
`Simplify`: no regression at 4, while at 2 the linear-algebra numbers stop
improving and `Integrate` gives back 7%.

A **literal** list is a separate matter and never packs, at any size:
`NDArrayQ[Table[1., {400}]]` is `True` and `NDArrayQ[{1., 2., 3., 4.}]` is
`False`. Packing is opt-in per producer, and a list written out in the source
has no producer to opt in — so `Inverse[{{4., 1.}, {1., 3.}}]` takes the
ordinary path whatever the threshold is. `ToNDArray` packs it on request.

Set `MATHILDA_PACK_MIN` in the environment to override the threshold for a
session.

Nested producers give a genuine rank-N array rather than a list of buffers:
`Table[i j, {i, 300}, {j, 300}]` is one rank-2 packed array.

## Exact integers

An all-`Integer` list packs to an `int64` buffer, and every operation on one is
either exactly what the ordinary list gives or falls back to it.

```mathematica
In[6]:= b = Range[1000]; {DataType[b], Head[b[[1]]]}
Out[6]= {"int64", Integer}

In[7]:= {Total[b], Mean[b], Median[b]}
Out[7]= {500500, 1001/2, 1001/2}

In[8]:= Total[Range[10^6]^3]          (* promotes past int64, not wrapped *)
Out[8]= 250000500000250000000000

In[9]:= {Take[Range[300] 5/2, 2], Take[Range[300] 2.5, 2], Take[Range[300]^-1, 2]}
Out[9]= {{5/2, 5}, {2.5, 5.}, {1, 1/2}}
```

`Sin[Range[3]]` is `{Sin[1], Sin[2], Sin[3]}` — symbolic, as it is for the
ordinary list. Where a buffer cannot hold the exact answer — an overflow, a
`Rational`, a radical, a symbolic result — the operation abandons and the list
implementation runs.

The same rule holds for a **visible** `NDArray`, and it has to be enforced in a
different place. A packed list is protected by the transparency gate, which
materialises an `int64` buffer for any head that has not declared itself exact
on one; a visible array is deliberately *not* gated, so the guard lives in the
kernel engine instead. An `int64` array handed to a kernel with no exact integer
arm declines and the list path answers:

```mathematica
In[10]:= Sin[NDArray[{1, 2, 3}, DataType -> "int64"]]
Out[10]= {Sin[1], Sin[2], Sin[3]}

In[11]:= Floor[NDArray[{1, 2, 3}, DataType -> "int64"]]   (* has an exact arm *)
Out[11]= NDArray[{1, 2, 3}]
```

Before 2026-08-01 the second rule was missing and `Out[10]` was
`NDArray[{0, 0, 0}]` — the double result truncated into the integer slot. Every
`real_closed` kernel did it. `Mod`, `Quotient`, `Floor`, `Ceiling`, `Round`,
`IntegerPart`, `Sign`, `UnitStep` and `Abs` have exact integer arms and keep
their buffer fast path.

## `ToNDArray`

`ToNDArray[list]` returns `list` stored as a dense buffer, ignoring the size
threshold that automatic packing applies. `ToNDArray[list, DataType -> "..."]`
forces the element type.

Returns `list` unchanged — not an error — when it cannot be packed:

| declines when | because |
|---|---|
| not rectangular | a buffer has one shape |
| empty | rank-0 has no buffer form |
| mixed exact and inexact (`{1, 2.5}`) | one dtype cannot give element 1 an `Integer` head and element 2 a `Real` head |
| any `Rational`, `BigInt`, arbitrary-precision, `Complex`, or symbolic element | no machine representation, or (for `Complex`) no faithful round trip yet |

Inferred dtypes: all-`Integer` → `"int64"`, all-`Real` → `"float64"`. An
explicit `DataType` may widen an exact list to a float buffer, but never rounds
an inexact one into an integer buffer:

```mathematica
In[10]:= DataType[ToNDArray[{1, 2, 3}, DataType -> "float64"]]
Out[10]= "float64"

In[11]:= NDArrayQ[ToNDArray[{1., 2.5}, DataType -> "int64"]]
Out[11]= False
```

`ToNDArray` on an already-packed list is a no-op; on an `NDArray[...]` it
restates the same values as a `List`.

## `ToPackedArray`

`ToPackedArray` is `ToNDArray` under Mathematica's name for the same operation —
the same builtin registered twice, not a rule that rewrites to the other, so it
costs no extra evaluation pass and cannot be shadowed. Every form and every
option is identical.

```mathematica
In[10]:= ToPackedArray[{1., 2., 3.}] === ToNDArray[{1., 2., 3.}]
Out[10]= True
```

## `FromNDArray`

`FromNDArray[expr]` undoes buffer storage: a packed list becomes an ordinary
list of separate elements, and an `NDArray[...]` becomes the nested list of its
entries. Anything else is returned unchanged.

```mathematica
In[12]:= NDArrayQ[FromNDArray[ToNDArray[{1., 2., 3.}]]]
Out[12]= False

In[13]:= FromNDArray[NDArray[{1., 2.}]]
Out[13]= {1., 2.}
```

`Normal` does the same for both forms.

## `FromPackedArray`

`FromPackedArray` is `FromNDArray` under Mathematica's name, on the same terms as
`ToPackedArray` above. It undoes both forms of buffer storage — a packed `List`
and an explicit `NDArray[...]`.

## What stays packed

Operations with a buffer-level implementation keep their result packed;
everything else produces an ordinary list with the same value.

Packed in, packed out: `Plus`, `Times`, `Power`, `Subtract`, `Divide`, `Dot`, the elementary and
special functions (`Sin`, `Exp`, `Gamma`, …), `Total`, `Mean`, `Min`, `Max`,
`MinMax`, `Median`, `Variance`, `Accumulate`, `Sort`, `Reverse`, `Transpose`, `Flatten`,
`Take`, `Drop`, `Partition`, `RotateLeft`/`RotateRight`, `Riffle`, `Join`,
`Differences`, `Ratios`, `Clip`, `Ramp`, `First`, `Last`, `Most`, `Rest`, `Part`,
`Extract`, `Append`, `Prepend`, `Catenate`, `TakeLargest`, `TakeSmallest`,
`Map`, `Select`, `TakeWhile`, `FoldList`, `Outer`, `MapThread`, `Inner`,
`Union`, `Intersection`, `Complement`, `Commonest`, `Rescale`, `MatrixPower`,
`PseudoInverse`, `LeastSquares`, `Mod`, `Quotient`, `RandomSample`,
`RandomChoice`, and the integer heads `GCD`, `LCM`, `DivisorSigma`, `EulerPhi`,
`MoebiusMu`, `IntegerLength`, `PowerMod`, `Prime`.

Four more read the buffer without returning one, because their answer is not a
uniform machine array: `Positive`/`Negative`/`NonNegative`/`NonPositive` answer
with a list of `True`/`False` (there is no boolean buffer), `IntegerDigits` with
a ragged list of digit lists, and `Counts` with an `Association`. They still
read the elements straight out of the buffer rather than materialising it.

`Subtract` and `Divide` are on that list without a buffer path of their own:
each rewrites to `Plus`/`Times`/`Power` without ever reading an element, so all
they needed was for the gate to stop firing one step before the head that can
use the buffer. `Divide` is deliberately *not* exact on an integer buffer —
`Range[10]/2` is a list of exact `Rational`s, so an integer argument
materialises and the ordinary list gives them.

A **symbolic** operand no longer costs a packed list its answer. `{1., 2., 3.}
+ x` threads to `{1. + x, 2. + x, 3. + x}` at every length; the same expression
used to come back unevaluated once the list was long enough to pack. A
*visible* `NDArray[...]` still declines with `NDArray::sym`, which is its
contract — naming that head asks for a purely numeric object.

`Ratios` keeps the buffer only for inexact data. `Ratios[{1, 2, 3}]` is
`{2, 3/2}` — exact `Rational`s, which no buffer holds — so an integer argument
takes the ordinary path. `TakeLargest` and `TakeSmallest` likewise decline an
integer buffer, because ordering it through `double` would compare two integers
past 2^53 equal; they use a bounded heap on real data, so the cost is
`O(n log k)` rather than a full sort.

`Append` and `Prepend` refuse an element that does not belong to the buffer's
dtype rather than coercing it. `Append[{1., 2.}, 0]` is the mixed
`{1., 2., 0}` — an exact `0` after two `Real`s, as `Mathematica` gives — and no
uniform buffer holds that, so the result is an ordinary list. `Append[{1., 2.},
0.]` stays packed.

`RandomSample` and `RandomChoice` draw from the same generator sequence the
ordinary path uses, so `SeedRandom[n]` gives the same answer whether the
argument is packed or not.

`Mod` and `Quotient` are worth naming explicitly because until 2026-08-01 they
were the break in the chain for integer work: both had a kernel that computed in
`double`, so neither could keep an `int64` buffer, and
`Mod[Range[10^6] 7919, 1000]` handed 10⁶ separate `Integer`s to whatever came
next. The heads downstream — `Union`, `Tally`, `DeleteDuplicates`, `Sort`,
`BinCounts` — were already buffer-capable and simply never received one. They
have exact integer arms now, and `Quotient` narrows: `Quotient[{1., 2., 3.}, 3]`
is the exact `{0, 0, 1}`, as `Mathematica` gives, while `Mod` keeps the
argument's exactness (`Mod[{1., 2., 3.}, 3]` is `{1., 2., 0.}`).

`ListConvolve` and `ListCorrelate` read a packed kernel, a packed list, or both,
and produce a packed real result; exact data, a symbolic kernel and a custom
`g`/`h` pair take the ordinary path and give the same answer.

`Fold` and `FoldList` use the buffer when the operator is one of `Plus`,
`Times`, `Max`, `Min` — written either as the bare symbol or as the equivalent
pure function (`Max[#1, #2] &`) — and the seed matches the buffer's exactness.
Any other operator still runs on the buffer through the numeric compiler, so a
general linear recurrence such as an exponential moving average stays packed
too:

```
In[1]:= v = RandomReal[{0, 1}, 10^6];
        NDArrayQ[FoldList[Max, First[v], Rest[v]]]
Out[1]= True
```

`Clip` keeps the buffer when both bounds are `Real`, and when an exact bound
clips nothing. It does **not** when an exact bound is actually reached, because
`Clip` returns the bound itself and the result is then a mixture of exact and
inexact numbers, which no uniform buffer holds:

```
In[2]:= Clip[{-2., 0., 2.}, {-1., 1.}]
Out[2]= {-1., 0., 1.}                     (* packed *)

In[3]:= Clip[{-2., 0., 2.}, {-1, 1}]
Out[3]= {-1, 0., 1}                       (* exact bounds reached; ordinary list *)
```

`First` and `Last` of a rank-1 packed array return a scalar with the head the
element has — an exact `Integer` from an integer buffer — and of a higher-rank
array return the sub-array with the leading axis dropped. `Rest` and `Most` of a
single-row array are `{}` and take the ordinary path.

`Outer` and `MapThread` cover the machine cases and hand the rest back:
`Outer[f, a, b]` uses the buffer for `Plus`/`Subtract`/`Times`/`Min`/`Max` on
real data, and `MapThread[f, arr]` for `Plus`/`Times`/`Min`/`Max` on real **or
integer** data. Any other `f`, or three or more arrays for `Outer`, takes the
ordinary path and gives the same answer. An integer `MapThread` of `Plus` or
`Times` that overflows `int64` abandons the whole array so the ordinary path
answers exactly, in GMP — never a wrapped sum.

`Part` gathers through a packed index list as well as a plain one:
`x[[idx]]` where `idx` is an integer packed list — which is what `Flatten`,
`Range`, `RandomInteger` and any arithmetic on them produce — reads the buffer
directly. This is the operation every sparse-matrix and graph kernel is built
out of, and until it read the index buffer it degraded the whole `Part` and
materialised *both* arrays. Positions must be integers and rank 1; a `Real`
position, an out-of-range one, a rank-2 index array or an empty list all take
the ordinary path and give the ordinary answer.

`Union`, `Intersection` and `Complement` merge **integer** buffers directly —
a sort and a linear merge, where the ordinary path allocates an expression per
element and sorts through `expr_compare`. Real data takes the ordinary path on
purpose: `0.` and `-0.` compare equal and print differently, so which of two
equal elements survives is a question the two paths must not be allowed to
answer differently. A `SameTest` option, symbolic elements, or a mixture of
packed and plain operands also take the ordinary path.

`Ramp` keeps a `Real` buffer packed. Its zero carries the *argument's*
exactness, so unlike `Clip` there is no bound whose head could leak into the
answer and nothing to gate. An integer buffer materialises.

`Plus` and `Times` also thread a lower-rank operand across the leading axes on
the buffer, which is what `matrix - rowVector` means:

```
In[1]:= m = RandomReal[{0, 1}, {3, 100000}];
        NDArrayQ[m - {0.1, 0.2, 0.3}]
Out[1]= True
```

`Tally` reads the buffer without producing one: its result is a list of
`{value, count}` pairs, which does not nest into a buffer, but it hashes the
machine words rather than materialising an expression per element. Only the
one-argument form — `Tally[list, test]` has to show the elements to `test`.
(An `int64` tally is the exception: both halves of each pair are integers, so it
comes back as a rank-2 packed matrix.)

`Commonest` counts with the same routine and *does* produce a buffer, because
every element it returns was copied out of one — so `Commonest[packed]` and
`Commonest[packed, n]` keep the input's dtype and stay packed. `Commonest[a, n]`
and `Commonest[a, UpTo[n]]` are both on the fast path; a count that is neither an
`Integer` nor `UpTo[Integer]` takes the ordinary path, as does a rank ≥ 2
argument, whose *elements* are rows rather than machine words.

Two limits worth knowing:

- **A list of packed lists becomes one array, or none.** Building `{p, q}` from
  two packed vectors of the same length and element class gives a single rank-2
  packed list, not two buffers inside an ordinary list — that second shape is not
  representable and is what the transparency gate exists to prevent. When the
  rows do not agree (ragged, or one exact and one inexact) the result is an
  ordinary list of ordinary lists, as before. Either way the value is the same:

  ```
  In[2]:= p = Range[1., 300.]; {NDArrayQ[{p, p}], Dimensions[{p, p}]}
  Out[2]= {True, {2, 300}}

  In[3]:= NDArrayQ[{p, Range[1., 299.]}]          (* ragged: declines *)
  Out[3]= False
  ```

  This matters most for a function that returns several arrays: without it every
  one of them is materialised at the `return`, and the caller pays for it on the
  *next* operation rather than on any it can see.
- Every other head materialises its packed arguments and runs normally. That is
  what makes `Count`, `Cases`, `Position`, `Level`, `ReplaceAll`, `Insert`,
  `Append`, pattern matching, and user-defined rules correct on a packed list
  without any of them knowing packing exists.

## Across a compiled call

A `CompiledFunction` borrows a packed argument and returns a packed result, so
packing survives a compiled call rather than being undone at each end.

```mathematica
In[13]:= f = Compile[{{u, _Real, 1}}, u^2 + 1.];
         r = f[Range[1., 200000.]]; {NDArrayQ[r], Head[r]}
Out[13]= {True, List}
```

The result's presentation follows its input, because a compiled call has to answer
with the head the interpreter would give for the same value:

| argument | result |
|---|---|
| plain `List` | plain `List` |
| packed `List` | packed `List`, at any size |
| `NDArray[...]` | `NDArray[...]` |

"At any size" is the *derived* rule, not the producer rule: `Sin[packedList]` is
packed however short it is, and a compiled call is no different. A body that
**builds** its array (`ConstantArray`, `Table`, `NestList`) has nothing to inherit
and follows the producer rule instead — the threshold and `$AutoArrayPacking`.

A **complex** result is never packed, for the same reason automatic packing
refuses complex lists: `Complex[re, 0.]` is not a form the evaluator produces.

A dtype mismatch costs one O(n) cast rather than a fall back to the interpreter,
so `Range[n]` (which infers `int64`) works at a `_Real` parameter. The cast
declines rather than rounds — an `int64` magnitude past 2^53 has no exact
`double`, and narrowing a float into an `_Integer` slot is a value change — and
the interpreter then gives the answer. Overflow inside a compiled body still
abandons and promotes to a bigint.

`RuntimeAttributes -> {Listable}` threads by rank and repacks by re-sniffing the
element type, so an integer-valued body over a packed integer list comes back as
`Integer`s.

### Which array heads compile

A head with a buffer path in the interpreter should also have one inside
`Compile[]`, and the reason is not speed alone: the compilable subset is a
**cliff**, so a single unlowered head sends the *whole* body to the interpreter
and takes every other head in it along. `Compile[{{v, _Real, 1}},
Total[v]/Mean[v]]` used to lose the compiled `Total` because of the `Mean`.

These are **delegated** — the VM calls the interpreter's own NDArray entry
point, so the compiled answer is bit-identical, rounding included:

| shape | heads |
|---|---|
| array → array | `Reverse` `Sort` `Ordering` `Accumulate` `Flatten` `Transpose` `Diagonal` `Take` `Drop` `Differences` `Ratios` `Most` `Rest` `Clip` `RotateLeft` `RotateRight` `MovingAverage` `MovingMedian` `TakeLargest` `TakeSmallest` `Inverse` `Normalize` `MatrixPower` `ReverseSort` `ConjugateTranspose` `PseudoInverse` |
| array → scalar | `Total` `Length` `Mean` `Median` `Variance` `StandardDeviation` `RootMeanSquare` `Max` `Min` `Tr` `Det` `MatrixRank` `Norm` |
| two arrays → array | `Dot` (matrix) `LinearSolve` `Cross` `LeastSquares` `ListConvolve` `ListCorrelate` `Join` |
| two arrays → scalar | `Dot` (vector·vector inner product) |
| elementwise | every registered kernel, including the narrowing ones (`Floor`, `Ceiling`, `Round`, `Sign`, `IntegerPart`, `UnitStep`) and the exact-integer ones (`Mod`, `Quotient`, `GCD`, `LCM`, `DivisorSigma`, `MoebiusMu`, `EulerPhi`, `IntegerLength`) |

The statistics reductions are **real element type only**, because their exact
answer is not a machine number: `Mean[{1, 2}]` is `3/2` and `Variance[Range[10]]`
is `55/6`. An integer vector declines and the interpreter answers exactly. `Max`
and `Min` are the exceptions — they *select* an element, so an integer vector
gives an `Integer`. The linear-algebra heads carry the same kind of gate per
head: `Inverse`, `Normalize`, `MatrixPower`, `Dot`, `LinearSolve`, `Cross` accept
real or complex (an int matrix inverse is exact Rationals → declines);
`LeastSquares`/`PseudoInverse` are real only; `ReverseSort`/`ConjugateTranspose`/
`Join` preserve any dtype. `Dot` delegates through a BLAS-first path
(`dgemm`/`dgemv`/`ddot`), so a compiled matrix product is as fast as the REPL's.

Anything a delegated head cannot handle comes back as a materialised `List`,
which the VM reads as "not a buffer" and hands to the interpreter — so
`Mod[v, 0]` and `Ratios` of an integer vector are as correct compiled as
interpreted, just not fast. `make check-compile-coverage` reports what is still
unlowered.

## Turning it off — `$AutoArrayPacking`

`$AutoArrayPacking` is `True` by default. Setting it to `False` stops every
producer from packing; setting it back to `True` resumes.

```mathematica
In[11]:= $AutoArrayPacking = False; NDArrayQ[Range[1., 300.]]
Out[11]= False

In[12]:= $AutoArrayPacking = True; NDArrayQ[Range[1., 300.]]
Out[12]= True
```

It does **not** disable `ToNDArray` / `ToPackedArray`, nor the explicit
`NDArray[...]` head — those are the user asking, not the system guessing. Only
`True` or `False` is accepted; anything else is refused with an
`$AutoArrayPacking::flagset` message and the symbol is rolled back to the live
state, so reading it back never lies about which path is running.

The environment variable `MATHILDA_NO_PACK` starts a session with it off, and
`$AutoArrayPacking` then reads back `False`.

Because packing never changes an answer, a run with and without it should agree
exactly. That equivalence is what the differential test suite checks: it flips
the same switch from C and diffs every output.

The companion switch for the other invisible optimisation is
`$AutoCompilation` — see [`control-flow.md`](control-flow.md#compile--compiledfunction).

## See also

- [`control-flow.md`](control-flow.md#compile--compiledfunction) —
  `$AutoCompilation`, the same kind of switch for the bytecode compiler that runs
  behind `Plot`, `Table` and friends.
- [`linear-algebra.md`](linear-algebra.md) — the visible `NDArray[...]` type,
  `DataType`, `NDArrayQ`, and the buffer-level linear algebra.
- `docs/design/packed_arrays.md` — the representation, the transparency gate,
  and the exactness audit.
