# Linear Algebra

## NDArray
A first-class, visibly-distinct dense machine-precision N-dimensional array,
modelled on numpy's `ndarray`.
- `NDArray[list]`: packs a rectangular, machine-precision nested list of any
  rank into a dense `EXPR_NDARRAY` value, storing a flat row-major (C-order)
  numeric buffer plus rank/dims/dtype directly on the node instead of a nested
  `List[List[...]]` tree. N-dimensional indexing is computed against the flat
  buffer.
- `NDArray[list, DataType -> "float32"]`: packs at a chosen element data type
  (dtype), analogous to numpy's `dtype`. Four types, mapping one-to-one onto
  BLAS's `s`/`d`/`c`/`z` precisions (the suffix is the per-component bit width):
  `"float64"` (`double`, the default), `"float32"` (`float`), `"complex64"`
  (2×`double`), `"complex32"` (2×`float`). Real dtypes accept `Integer`/`Real`
  leaves; complex dtypes additionally accept `Complex[...]` (and bare reals).
  With no `DataType` option the dtype is inferred from the data: complex leaves
  give `"complex64"`, otherwise `"float64"` (`"float32"`/`"complex32"` require an
  explicit option). `Options[NDArray]` is `{DataType -> "float64"}`; `DataType[a]` gives an
  array's dtype as a string. dtype is part of an array's identity, so
  `NDArray[{1,2},DataType->"float32"]` is not `SameQ` to the `"float64"` one.
- A rectangular list with a zero-length axis builds an **empty array** of that
  shape: `NDArray[{}]` has `Dimensions` `{0}` and `Length` `0`;
  `NDArray[{{}, {}}]` has `Dimensions` `{2, 0}`. With no elements to infer from,
  the dtype is the requested one (`"float64"` by default). Elementwise
  operations keep it empty, and reductions answer as on the empty List:
  - `Total[NDArray[{}]]` is `0.`;
  - `Max` and `Min` are `-Infinity` and `Infinity`;
  - `Median[NDArray[{}]]` stays unevaluated, as `Median[{}]` does.

  Only this constructor builds one: automatic packing never packs an empty
  List, so `Range[0]` and an empty `Select` result stay ordinary Lists.
- A visible `NDArray` meeting a plain List in a Listable call threads as its
  List does whenever the List cannot join the buffer: `NDArray[{1., 2., 3.}] +
  {1, 2}` emits `Thread::tdlen`, and `NDArray[{1., 2.}] + {x, y}` is
  `{1. + x, 2. + y}`. A same-shape numeric List still takes the buffer path
  (`NDArray[{1., 2.}] + {10, 20}` is `NDArray[{11., 22.}]`).
- `NDArray[...]` is always exactly what it says: `Head[NDArray[{{1, 2}, {3,
  4}}]]` is `NDArray`, never `List`; `MatrixQ`, `VectorQ` and `ListQ` all report
  `False` for it; and it always prints as `NDArray[{{1.0, 2.0}, {3.0, 4.0}}]`,
  never as bare `{{1.0, 2.0}, {3.0, 4.0}}`.
- The *invisible* form of the same storage is a **packed list** -- an ordinary
  `List` that happens to be held as a dense buffer, which `Head`, `ListQ`,
  printing and pattern matching cannot distinguish from any other `List`. That
  is what `ToNDArray` produces and what automatic packing uses; see
  [`packed-arrays.md`](packed-arrays.md). `NDArrayQ` is the one predicate that
  answers `True` for both forms, since it asks about storage rather than about
  what the value claims to be.
- **Every linear-algebra routine accepts `NDArray` input** (see
  `src/linalg/ndlinalg.c`). The heavy/common consumers — `Det`, `Inverse`,
  `LinearSolve`, `MatrixRank`, `Tr`, `Norm`, `Normalize`, `Cross` — take a fast
  C-level path over the flat buffer and return a closed-system `NDArray` (or a
  scalar `Real`/`Complex`); `Diagonal` (`src/ndstruct.c`) likewise slices the
  k-th diagonal straight off the buffer; the real `Det`/`Inverse`/`LinearSolve`/`MatrixRank`
  path uses an in-house partial-pivot `double` LU, so `Det[NDArray[...]]` is
  O(n³) even without LAPACK. The factorisations, predicates, constructors, and
  lattice routines (`Eigenvalues`, `QRDecomposition`, `SingularValueDecomposition`,
  `PositiveDefiniteMatrixQ`, `DiagonalMatrix`, `LatticeReduce`, ...) convert the
  `NDArray` argument to a `List` and re-evaluate, routing through their existing
  numeric kernels.
- `Dot`, elementwise `Plus`/`Times`/`Power` recognize `NDArray` operands and use
  a fast C-level path that loops directly over the raw buffers, per dtype,
  skipping symbolic construction per element. The result dtype follows a
  promotion lattice (complex dominates real; 64-bit component dominates 32-bit),
  so e.g. `float32 + float64 -> float64` and any complex operand yields a complex
  result. The float64 real case keeps a dedicated tight `double` loop (no
  regression). `Power` promotes to complex when a real base with a non-integer
  exponent leaves the real axis (`NDArray[{-1.0}]^0.5` is complex); integer
  exponents on complex bases are exact (`I^2 == -1`).
- **numpy-style scalar broadcasting.** A numeric scalar combines with an array
  elementwise: `1 + NDArray[{1,2,3}]`, `3 * NDArray[...]`, `2 ^ NDArray[...]`,
  and `-A` / `A - B`. Scalars are "weak": an `Integer`/`Real` scalar keeps the
  array's float width (`float32 + 1 -> float32`), while a `Complex` scalar moves
  it onto the complex axis at the same width (`float32 + I -> complex32`).
- **Elementary and special functions map over `NDArray` at C speed.** `Sin`,
  `Cos`, `Exp`, `Log`, `Sqrt`, the inverse trig/hyperbolic functions, `Abs`/`Re`/
  `Im`/`Arg`/`Conjugate`/`Sign`, and special functions (`Gamma`, `Erf`, `Erfc`,
  `LogGamma`, `BesselJ`/`BesselY`, `Beta`, …) loop the flat buffer directly (no
  per-element `Expr`). A real input that leaves the real axis promotes to a
  complex array (`Log[NDArray[{-1.,4.}]]`, `ArcSin` outside `[-1,1]`); projections
  (`Abs`/`Re`/`Im`/`Arg`) yield a real array even from complex input. Two-arg
  forms `Log[b, NDArray[...]]`, `ArcTan[NDArray[...], y]`, `BesselJ[n,
  NDArray[...]]`, etc. broadcast the scalar. Any element the C kernel cannot
  produce at machine precision (a pole, e.g. `Gamma[NDArray[{0.,1.}]]`, or a
  function with no machine kernel yet, e.g. `FresnelC`, `AiryAi`, `ProductLog`,
  `Zeta`) degrades faithfully to the `List` path, so the result always matches
  `f[{...}]`. Registration is one line per head via
  `symtab_set_ndarray_unary_kernel` / `_binary_kernel` (`src/ndkernels.c`).
- Machine precision only. Any input or fast-path result that would need a
  non-machine-precision entry (a symbol, an exact/rational number, a BigInt,
  an MPFR value) is left unpacked/unevaluated instead of forcing a lossy
  conversion -- `NDArray[...]` degrades to an ordinary nested `List` rather
  than hiding a representation change from the user.
- The linalg BLAS/LAPACK bridges accept all four dtypes as input and return
  results as dtyped `NDArray`s (real results as a real array, complex results
  as a `complex64` array), keeping `NDArray` a closed system under those ops.
- `Normal[NDArray[...]]` converts back to the equivalent nested `List`.

**Introspection** (numpy correspondence):
- `Dimensions[NDArray[...]]` is the `shape` -- reads `rank`/`dims` directly
  (O(1)), rather than recursively probing shape as it does for a nested `List`.
- `ArrayDepth`-style rank is given by `Depth[NDArray[...]]`, which is `ndim + 1`
  (matching the equivalent nested `List`; the `+1` counts the atom level).
- `Length[NDArray[...]]` is the leading-axis length (numpy `len`, = `shape[0]`).
- `NDArrayQ[expr]` is `True` iff `expr` is an `NDArray` value.
- `DataType[NDArray[...]]` is the element dtype string (numpy `.dtype`);
  unevaluated on a non-array.
- `MatrixQ`/`VectorQ` answer per rank (2 / 1).
- `Part[NDArray[...], i1, ...]` (`arr[[...]]`) indexes the flat buffer directly
  and stays packed. Subscripts are 1-based with negatives counting from the end;
  omitted trailing axes are implicitly `All`. An integer subscript drops that
  axis, `All` / `Span[a,b,c]` (step and reversal supported) / a list of
  positions keep it. Indexing every axis with integers gives a scalar leaf
  (`Real`, or `Complex` for a complex dtype); otherwise the result is a
  sub-`NDArray` whose shape is the kept axes — so `mat[[i]]` is a packed row,
  `mat[[All, j]]` a packed column, `mat[[All, a;;b]]` a packed sub-matrix, all
  preserving dtype. `arr[[0]]` is the head `NDArray`. Out-of-range or too-many
  subscripts leave the `Part[...]` unevaluated, as in the reference. (Exotic specs
  such as an empty span or `Key[...]` fall back to the general `List` `Part`.)

**Features**:
- `Protected`.
- A ragged (non-rectangular) `list` — unequal sublist shapes, or a mix of
  list and non-list siblings — can never form an array, so `NDArray[list]`
  prints a one-line `NDArray::ragged` warning and stays unevaluated. An empty
  list, a non-machine-precision entry (e.g. a symbol), or a non-list argument
  stays unevaluated silently (the symbolic case may become packable after
  further evaluation).
- `Dot[NDArray[a], NDArray[b]]` contracts the trailing axis of `a` with the
  leading axis of `b` over raw doubles for rank <= 2 operands, giving a new
  `NDArray` (or a bare machine `Real` for a vector.vector contraction). Falls
  back to converting through `Normal` and using the generic tensor path for
  higher-rank operands or a rank mismatch; a genuine shape mismatch (inner
  dimensions disagree) prints `Dot::dotsh` and leaves the call unevaluated.
- `NDArray[a] + NDArray[b]` / `NDArray[a] * NDArray[b]` compute elementwise
  `+`/`*` over raw doubles when both operands are `NDArray` values of
  identical shape. When the operands are all `NDArray` values but of
  disagreeing shape, a one-line `NDArray::shape` warning is printed (naming the
  two shapes) and the sum/product is left unevaluated, mirroring `Dot::dotsh`.
  A mixed `NDArray` + scalar/other operand set instead falls through to the
  generic symbolic `Plus`/`Times` path, treating the `NDArray` as an opaque
  term. numpy-style broadcasting (scalar/array, shape-compatible) is not yet
  implemented.
- Because an `NDArray` is purely numeric, combining one with a **symbolic**
  operand (a bare symbol or any non-numeric expression) can never be carried
  out elementwise. `Plus`/`Times`/`Power` print a one-line `NDArray::sym`
  warning and leave the expression unevaluated: `NDArray[{1., 3.}] + a`,
  `c NDArray[{1., 3.}]`, `NDArray[{1., 3.}]^n`. A numeric scalar operand
  (Integer/Real/Rational/Complex) still broadcasts silently and is unaffected.

```mathematica
In[1]:= NDArray[{{1, 2}, {3, 4}}]
Out[1]= NDArray[{{1.0, 2.0}, {3.0, 4.0}}]

In[2]:= Dimensions[NDArray[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]]
Out[2]= {2, 2, 2}

In[3]:= Depth[NDArray[{{1, 2}, {3, 4}}]]
Out[3]= 3

In[4]:= Dot[NDArray[{{1, 2}, {3, 4}}], NDArray[{{5, 6}, {7, 8}}]]
Out[4]= NDArray[{{19.0, 22.0}, {43.0, 50.0}}]

In[5]:= NDArray[{{1, 2}, {3, 4}}] + NDArray[{{5, 6}, {7, 8}}]
Out[5]= NDArray[{{6.0, 8.0}, {10.0, 12.0}}]

In[6]:= NDArrayQ[NDArray[{1, 2, 3}]]
Out[6]= True

In[7]:= NDArray[{{1, x}, {3, 4}}]
Out[7]= NDArray[{{1, x}, {3, 4}}]
```


## HermitianMatrixQ
Tests whether a matrix is explicitly Hermitian (self-adjoint).
- `HermitianMatrixQ[m]`: `True` if `m == ConjugateTranspose[m]` under
  structural equality, `False` otherwise.
- `HermitianMatrixQ[m, SameTest -> f]`: entries `m[i,j]` and
  `Conjugate[m[j,i]]` are treated as equal iff
  `f[m[i,j], Conjugate[m[j,i]]]` evaluates to `True`.
- `HermitianMatrixQ[m, Tolerance -> t]`: entries are accepted when
  `Abs[m[i,j] - Conjugate[m[j,i]]] <= t`.
- Diagonal entries must satisfy the same test, i.e. be self-conjugate
  (purely real for numeric matrices).

**Features**:
- `Protected`.
- Default test is structural: it accepts (Conjugate[a], a) / (a, Conjugate[a])
  symbolic pairs without requiring our `Conjugate` builtin to fold
  `Conjugate[Conjugate[x]] -> x`.
- Returns `False` (rather than leaving unevaluated) on non-matrix, non-square,
  ragged, empty, or higher-rank tensor inputs.
- Unknown options and non-Rule trailing arguments leave the call unevaluated.

```mathematica
In[1]:= HermitianMatrixQ[{{1, 3 + 4 I}, {3 - 4 I, 2}}]
Out[1]= True

In[2]:= HermitianMatrixQ[{{0, a, b}, {Conjugate[a], 1, c},
            {Conjugate[b], Conjugate[c], -1}}]
Out[2]= True

In[3]:= HermitianMatrixQ[{{1, 2 I}, {2 I, 3}}]
Out[3]= False

In[4]:= HermitianMatrixQ[{{1.0, 2.0 + 0.01 I}, {2.0 - 0.02 I, 1.5}},
            Tolerance -> 0.1]
Out[4]= True
```


## DiagonalMatrixQ
Tests whether a matrix has nonzero entries only on a particular diagonal.
- `DiagonalMatrixQ[m]`: `True` if every off-(main-diagonal) entry of `m`
  is zero, `False` otherwise.
- `DiagonalMatrixQ[m, k]`: `True` if every entry off the `k`-th diagonal
  (`m[[i, j]]` with `j - i != k`) is zero.  Positive `k` selects a
  superdiagonal above the main diagonal; negative `k` selects a
  subdiagonal below it.
- `DiagonalMatrixQ[m, Tolerance -> t]` (and `DiagonalMatrixQ[m, k,
  Tolerance -> t]`): off-diagonal entries `e` are taken to be zero when
  `Abs[e] <= t` evaluates to `True`.

**Features**:
- `Protected`.
- Works for rectangular matrices, not only square -- only the entry-zero
  predicate and the shape constraints matter.
- Default test is structural: only literal numeric zeros (`Integer 0`,
  `Real 0.0`, `BigInt 0`) count as zero.  Symbolic off-diagonal entries
  fail the test, so the predicate is conservative.
- Returns `False` (rather than leaving unevaluated) on non-list, scalar,
  vector, ragged, or higher-rank tensor inputs.  `{}` is rejected; an
  `n`-by-`0` matrix (e.g. `{{}, {}}`) is vacuously diagonal and returns
  `True`.
- Zero positional arguments emits a reference-compatible

  ```
  DiagonalMatrixQ::argt: DiagonalMatrixQ called with 0 arguments; 1 or 2 arguments are expected.
  ```

  to `stderr` and leaves the call unevaluated.

- More than two positional arguments (or any non-`Rule` junk in the
  option region) emits

  ```
  DiagonalMatrixQ::nonopt: Options expected (instead of <expr>) beyond position 2 in DiagonalMatrixQ[...]. An option must be a rule or a list of rules.
  ```

  to `stderr` and leaves the call unevaluated.

```mathematica
In[1]:= DiagonalMatrixQ[{{a, 0, 0}, {0, b, 0}, {0, 0, c}}]
Out[1]= True

In[2]:= DiagonalMatrixQ[{{1, 0, 0}, {0, 0, 2}, {3, 0, 0}}]
Out[2]= False

In[3]:= DiagonalMatrixQ[{{0, a, 0}, {0, 0, b}, {0, 0, 0}}, 1]
Out[3]= True

In[4]:= DiagonalMatrixQ[{{0, 0, 0}, {a, 0, 0}, {0, b, 0}}, -1]
Out[4]= True

In[5]:= DiagonalMatrixQ[{{1, 0, 0}, {0, 2, 0}}]
Out[5]= True

In[6]:= DiagonalMatrixQ[{{1, 0}, {0, 2}, {0, 0}}]
Out[6]= True

In[7]:= DiagonalMatrixQ[{{1., 10^-12, 0}, {0, 2., 10^-13}, {0, 0, 3.}},
            Tolerance -> 10^-12]
Out[7]= True

In[8]:= DiagonalMatrixQ[IdentityMatrix[5]]
Out[8]= True
```


## UpperTriangularMatrixQ
Tests whether a matrix is upper triangular relative to a chosen diagonal.
- `UpperTriangularMatrixQ[m]`: `True` if every entry strictly below the
  main diagonal of `m` is zero, `False` otherwise.
- `UpperTriangularMatrixQ[m, k]`: `True` if every entry `m[[i, j]]` with
  `j - i < k` is zero.  Positive `k` selects a superdiagonal above the
  main diagonal (the test becomes stricter -- the main diagonal must be
  zero too); negative `k` selects a subdiagonal below it (the test
  becomes more permissive -- the first `|k|` subdiagonals are allowed to
  be nonzero).
- `UpperTriangularMatrixQ[m, Tolerance -> t]` (and
  `UpperTriangularMatrixQ[m, k, Tolerance -> t]`): sub-diagonal entries
  `e` are taken to be zero when `Abs[e] <= t` evaluates to `True`.

**Features**:
- `Protected`.
- Works for rectangular matrices, not only square -- only the entry-zero
  predicate and the shape constraints matter.
- Default test is structural: only literal numeric zeros (`Integer 0`,
  `Real 0.0`, `BigInt 0`) count as zero.  Symbolic sub-diagonal entries
  fail the test, so the predicate is conservative.
- Returns `False` (rather than leaving unevaluated) on non-list, scalar,
  vector, ragged, or higher-rank tensor inputs.  `{}` is rejected; an
  `n`-by-`0` matrix (e.g. `{{}, {}}`) is vacuously upper triangular and
  returns `True`.
- Zero positional arguments emits a reference-compatible

  ```
  UpperTriangularMatrixQ::argt: UpperTriangularMatrixQ called with 0 arguments; 1 or 2 arguments are expected.
  ```

  to `stderr` and leaves the call unevaluated.

- More than two positional arguments (or any non-`Rule` junk in the
  option region) emits

  ```
  UpperTriangularMatrixQ::nonopt: Options expected (instead of <expr>) beyond position 2 in UpperTriangularMatrixQ[...]. An option must be a rule or a list of rules.
  ```

  to `stderr` and leaves the call unevaluated.

```mathematica
In[1]:= UpperTriangularMatrixQ[{{a, b, c}, {0, e, f}, {0, 0, g}}]
Out[1]= True

In[2]:= UpperTriangularMatrixQ[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[2]= False

In[3]:= UpperTriangularMatrixQ[{{0, 1, 2}, {0, 0, 3}, {0, 0, 0}}, 1]
Out[3]= True

In[4]:= UpperTriangularMatrixQ[{{1, 2, 3}, {4, 5, 6}, {0, 7, 9}}, -1]
Out[4]= True

In[5]:= UpperTriangularMatrixQ[{{1, 2, 3}, {0, 4, 5}}]
Out[5]= True

In[6]:= UpperTriangularMatrixQ[{{1, 2}, {0, 4}, {0, 0}}]
Out[6]= True

In[7]:= UpperTriangularMatrixQ[{{1., 2., 3.}, {10^-12, 4., 5.}, {0, 10^-13, 6.}},
            Tolerance -> 10^-12]
Out[7]= True

In[8]:= UpperTriangularMatrixQ[IdentityMatrix[5]]
Out[8]= True
```


## SquareMatrixQ
Tests whether a matrix has the same number of rows and columns.
- `SquareMatrixQ[m]`: `True` if `Dimensions[m] == {n, n}` for some `n >= 1`,
  `False` otherwise.

**Features**:
- `Protected`.
- Pure shape test: no element predicate or option is consulted.
- Works for symbolic as well as numerical matrices (`{{a,b},{c,d}}` is
  square; entries are not evaluated).
- Returns `False` (rather than leaving unevaluated) on non-list, scalar,
  vector, empty (`{}`, `{{}}`), ragged, rectangular, or higher-rank
  tensor inputs.
- Exactly one argument is accepted; any other count emits a
  reference-compatible

  ```
  SquareMatrixQ::argx: SquareMatrixQ called with N arguments; 1 argument is expected.
  ```

  to `stderr` and leaves the call unevaluated.

```mathematica
In[1]:= SquareMatrixQ[{{1, 2}, {3, 4}}]
Out[1]= True

In[2]:= SquareMatrixQ[{{1, 2, 3}, {4, 5, 6}}]
Out[2]= False

In[3]:= SquareMatrixQ[{1, 2, 3}]
Out[3]= False

In[4]:= SquareMatrixQ[{{1}, {2, 3}}]
Out[4]= False

In[5]:= SquareMatrixQ[{{a, b, c}, {d, e, f}, {g, h, i}}]
Out[5]= True
```


## SymmetricMatrixQ
Tests whether a matrix is explicitly symmetric (`m == Transpose[m]`).
- `SymmetricMatrixQ[m]`: `True` if every off-diagonal pair satisfies
  `m[[i,j]] == m[[j,i]]` under structural equality, `False` otherwise.
- `SymmetricMatrixQ[m, SameTest -> f]`: entries `m[i,j]` and `m[j,i]`
  are treated as equal iff `f[m[i,j], m[j,i]]` evaluates to `True`.
- `SymmetricMatrixQ[m, Tolerance -> t]`: entries are accepted when
  `Abs[m[i,j] - m[j,i]] <= t`.

**Features**:
- `Protected`.
- Default test is structural via `expr_eq`; the diagonal is exempt
  (trivially symmetric).
- Uses `m^T == m` for both real- and complex-valued matrices, so a
  complex symmetric matrix need not be Hermitian (and vice versa).
- Returns `False` (rather than leaving unevaluated) on non-matrix,
  non-square, ragged, empty, or higher-rank tensor inputs.
- Unknown options and non-`Rule` trailing arguments leave the call
  unevaluated.

```mathematica
In[1]:= SymmetricMatrixQ[{{1, 2}, {2, 3}}]
Out[1]= True

In[2]:= SymmetricMatrixQ[{{a, b, c}, {b, d, e}, {c, e, f}}]
Out[2]= True

In[3]:= SymmetricMatrixQ[{{1 + I, 2 - 3 I}, {2 - 3 I, 2 - 3 I}}]
Out[3]= True

In[4]:= SymmetricMatrixQ[{{1, 3 + 4 I}, {3 - 4 I, 2}}]   (* Hermitian, not symmetric *)
Out[4]= False

In[5]:= SymmetricMatrixQ[{{1, Log[x^2]}, {2 Log[x], 2}},
            SameTest -> (Simplify[#1 - #2, x > 0] == 0 &)]
Out[5]= True
```


## PositiveDefiniteMatrixQ
Tests whether a matrix is explicitly positive definite.
- `PositiveDefiniteMatrixQ[m]`: `True` if `Re[Conjugate[x] . m . x] > 0`
  for every nonzero vector `x`, `False` otherwise.

**Features**:
- `Protected`.
- Equivalent to: the Hermitian part `(m + ConjugateTranspose[m]) / 2`
  has only positive eigenvalues, equivalently admits a Cholesky
  factorisation with a real positive diagonal.
- On numeric matrices the test is performed by attempting Cholesky on
  the Hermitian part.  When `USE_LAPACK` is available the routine
  dispatches to BLAS/LAPACK's `dpotrf` (real) or `zpotrf` (complex);
  otherwise an in-house Cholesky is used.  Either returns `info == 0`
  iff the matrix is positive definite.
- Builds the Hermitian part regardless of whether the input is
  Hermitian, so e.g. a real matrix with non-symmetric entries is
  tested via `(m + m^T) / 2`.
- For symbolic or otherwise non-coercible entries the predicate
  conservatively returns `False`; "explicitly positive definite" is
  not proved symbolically.
- Returns `False` (rather than leaving unevaluated) on non-matrix,
  non-square, ragged, empty, or higher-rank tensor inputs.
- Exactly one argument is accepted; any other count emits a
  reference-compatible

  ```
  PositiveDefiniteMatrixQ::argx: PositiveDefiniteMatrixQ called with N arguments; 1 argument is expected.
  ```

  to `stderr` and leaves the call unevaluated.

```mathematica
In[1]:= PositiveDefiniteMatrixQ[{{5, -1}, {-1, 4}}]
Out[1]= True

In[2]:= PositiveDefiniteMatrixQ[{{2.3, -1.2}, {0.6, 3.7}}]
Out[2]= True

In[3]:= PositiveDefiniteMatrixQ[{{1, 2 I}, {-I, 4}}]
Out[3]= True

In[4]:= PositiveDefiniteMatrixQ[{{Pi, -5, 2}, {E, -3, -3}, {5, Sqrt[2], 5}}]
Out[4]= False

In[5]:= PositiveDefiniteMatrixQ[{{1, a}, {b, 2}}]
Out[5]= False

In[6]:= PositiveDefiniteMatrixQ[Table[1/(i + j - 1), {i, 8}, {j, 8}]]
Out[6]= True
```


## NegativeDefiniteMatrixQ
Tests whether a matrix is explicitly negative definite.
- `NegativeDefiniteMatrixQ[m]`: `True` if `Re[Conjugate[x] . m . x] < 0`
  for every nonzero vector `x`, `False` otherwise.

**Features**:
- `Protected`.
- Equivalent to: `-m` is positive definite, i.e. the negated Hermitian
  part `-(m + ConjugateTranspose[m]) / 2` has only positive eigenvalues
  and admits a Cholesky factorisation with a real positive diagonal.
- On numeric matrices the test is performed by attempting Cholesky on
  the negated Hermitian part.  When `USE_LAPACK` is available the
  routine dispatches to BLAS/LAPACK's `dpotrf` (real) or `zpotrf`
  (complex); otherwise an in-house Cholesky is used.  Either returns
  `info == 0` iff the matrix is negative definite.
- Builds the Hermitian part regardless of whether the input is
  Hermitian, so e.g. a real matrix with non-symmetric entries is
  tested via `-(m + m^T) / 2`.
- For symbolic or otherwise non-coercible entries the predicate
  conservatively returns `False`; "explicitly negative definite" is
  not proved symbolically.
- Returns `False` (rather than leaving unevaluated) on non-matrix,
  non-square, ragged, empty, or higher-rank tensor inputs.
- Exactly one argument is accepted; any other count emits a
  reference-compatible

  ```
  NegativeDefiniteMatrixQ::argx: NegativeDefiniteMatrixQ called with N arguments; 1 argument is expected.
  ```

  to `stderr` and leaves the call unevaluated.

```mathematica
In[1]:= NegativeDefiniteMatrixQ[{{-5, 1}, {1, -4}}]
Out[1]= True

In[2]:= NegativeDefiniteMatrixQ[{{-2.3, -1.2}, {0.6, -3.7}}]
Out[2]= True

In[3]:= NegativeDefiniteMatrixQ[{{-1, 2 I}, {-I, -4}}]
Out[3]= True

In[4]:= NegativeDefiniteMatrixQ[{{Pi, -5, 2}, {E, -3, -3}, {5, Sqrt[2], 5}}]
Out[4]= False

In[5]:= NegativeDefiniteMatrixQ[{{-1, a}, {b, -2}}]
Out[5]= False

In[6]:= NegativeDefiniteMatrixQ[Table[-1/(i + j - 1), {i, 8}, {j, 8}]]
Out[6]= True
```


## Dot (.)
Gives products of vectors, matrices, and tensors.
- `a . b` or `Dot[a, b]`
- `a . b . c` or `Dot[a, b, c]`

**Features**:
- `Flat`, `OneIdentity`, `Protected`.
- Contracts the last index in `a` with the first index in `b`.
- Applying `Dot` to a rank `n` tensor and a rank `m` tensor gives a rank `m+n-2` tensor.
- Scalar product of two vectors yields a scalar; the empty inner product
  `{} . {}` is the empty sum `0`.
- Product of a matrix and a vector yields a vector.
- Product of two matrices yields a matrix.
- When arguments are not lists, `Dot` remains unevaluated, and prints infix as
  `a . b` (OutputForm and InputForm; `FullForm` keeps `Dot[a, b]`), at the
  parser's precedence: `(a b) . c`, `(-a) . b`, `(1/2) . a`, `a^2 . b`,
  `(a . b)^2`, `a . b + c`.
- Gives an error message `Dot::dotsh` if the shapes of the inputs are not compatible.
- **Arbitrary precision.** Vectors and matrices of MPFR numbers (one operand
  all MPFR, the other MPFR or exact numbers, zeros included) contract as fused dot
  products: each output element is one dot product of the operands' mantissas
  (FLINT's `arb_approx_dot`, else `mpfr_dot`), rounded once, with large
  products split over threads. The answer has the smaller operand precision,
  as `Times` / `Plus` give. A 256 x 256 matrix times a vector at 800 digits
  takes ~0.02 s (0.35 s through per-element `Times` / `Plus` in v0.6225).
  An exact zero in the other operand is an MPFR zero at its partner's
  precision (as `0 x` is), so `m . DiagonalMatrix[w]` of an MPFR `m` is the
  fused path too: 64 x 64 at 20 digits 7 ms, was 1.4 s through `Times` /
  `Plus`.
- **Compilable.** Over machine arrays `Dot` lowers inside `Compile[]` and
  auto-compilation: matrix shapes through a BLAS-first path (`A_NDFN2`), the
  `vector·vector` inner product to a scalar (`V_NDFN2`), real and complex.
  Together with `Inverse`, `LinearSolve`, `Cross`, `LeastSquares`, `Normalize`,
  `MatrixPower`, `PseudoInverse`, and `ConjugateTranspose` — see
  [packed-arrays.md](packed-arrays.md#which-array-heads-compile).

```mathematica
In[1]:= {a, b, c} . {x, y, z}
Out[1]= a x + b y + c z

In[2]:= {{a, b}, {c, d}} . {x, y}
Out[2]= {a x + b y, c x + d y}

In[3]:= {{a, b}, {c, d}} . {{1, 2}, {3, 4}}
Out[3]= {{a + 3 b, 2 a + 4 b}, {c + 3 d, 2 c + 4 d}}
```

### Dot: exact rational operands
Two operands of rank <= 2 whose entries are all integers or rationals are
multiplied as one FLINT `fmpq_mat_mul` (a common denominator, then a
multimodular integer product) instead of a `Times` and a `Plus` evaluated per
element; the first non-rational entry declines to the generic path.  The
21 x 101 by 101 x 21 normal matrix of a degree-20 polynomial fit took 154 ms
per element-wise and takes 2 ms.

### Dot: numeric contagion
An operand mixing machine Reals and Complexes (the usual shape of a complex
factor), or exact numbers beside a machine Real, is packed for the BLAS path:
the product of two 166 x 166 such matrices takes 0.07 s (was 18 s).

## Inner
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

> **Packed arrays.** `Inner[Times, a, b, Plus]` **is** a `Dot`, so with two
> buffers it runs as one — the same BLAS path `a . b` takes. Any other
> operator pair has no kernel and takes the ordinary path with the same
> answer.

## Outer
Gives the generalized outer product of the `listi`, forming all possible combinations of the lowest-level elements in each of them, and feeding them as arguments to `f`.
- `Outer[f, list1, list2, ...]`
- `Outer[f, list1, list2, ..., n]`: treats as separate elements only sublists at level `n`.
- `Outer[f, list1, list2, ..., n1, n2, ...]`: treats as separate elements only sublists at level `ni` in the corresponding `listi`.

**Features**:
- `Protected`.
- Applying `Outer` to two tensors of ranks $r$ and $s$ gives a tensor of rank $r+s$.
- The heads of all `listi` must be the same, but need not necessarily be `List`.
- **Compiled cells.** When every `listi` is a flat List of machine reals, or of
  MPFR reals at one precision, and applying `f` is a pure substitution -- a
  `Function[{x, y, ...}, body]`, or a symbol whose only definition is
  `f[x_, y_, ...] := body` with plain distinct blanks (no condition, not
  `Orderless` / `Flat`) -- the body is compiled once (auto-compilation,
  machine or at that working precision) and run per cell. The values are the
  interpreter's bit for bit; a cell the compiled code declines (complex,
  non-finite) is evaluated by the interpreter. `Outer[K, pts, pts]` with
  `K[x_, y_] := E^(x + y + x^2 + y^2 + x y + x^2 y^2)` on 64 points: 1.3 ms
  machine (was 40 ms), 47 ms at 20 digits (was 260 ms).

```mathematica
In[1]:= Outer[f, {a, b}, {x, y, z}]
Out[1]= {{f[a, x], f[a, y], f[a, z]}, {f[b, x], f[b, y], f[b, z]}}

In[2]:= Outer[Times, {1, 2, 3, 4}, {a, b, c}]
Out[2]= {{a, b, c}, {2 a, 2 b, 2 c}, {3 a, 3 b, 3 c}, {4 a, 4 b, 4 c}}

In[3]:= Outer[g, f[a, b], f[x, y, z]]
Out[3]= f[f[g[a, x], g[a, y], g[a, z]], f[g[b, x], g[b, y], g[b, z]]]
```
## Det
Gives the determinant of a square matrix.
- `Det[m]`

**Features**:
- `Protected`.
- Evaluates the determinant of a square matrix symbolically or numerically using Laplace expansion.
- **Symbolic matrices** (v0.6124) use cofactor expansion with each minor memoised by its column set, and the minor of an exact-zero entry is skipped. A dense n×n costs n·2ⁿ products instead of n!, and a sparse one far fewer. A 10×10 symbolic diagonal matrix went from 1.6 s to 41 µs. Answers are unchanged.
- **Expanded answer** (v0.6139). A symbolic determinant is returned expanded, as the reference returns it: `Det[{{a, b, c}, {d, e, f}, {g, h, k}}]` is `-(c e g) + b f g + ...` and `Det[{{a + b, c}, {d, e}}]` is `-(c d) + a e + b e`. A matrix with a sum in a denominator keeps the cofactor form. Expand would distribute every numerator over that denominator, and a 6×6 matrix of rational functions overflowed. Only `Det` itself expands; the adjugate in `Inverse` and `Cross` keep their forms.
- Returns a warning `Det::matsq` if `m` is not a non-empty square matrix.
- **FLINT acceleration** (when built with FLINT): a matrix whose entries are all integer or rational is computed exactly via `fmpq_mat_det` in polynomial time, avoiding the `O(n!)` Laplace expansion (e.g. a 12×12 Hilbert determinant is instant and exact). Symbolic matrices fall through to Laplace. The same kernel is exposed directly as `` FLINT`Det `` (see the FLINT` context section in *Structural Manipulation*).
- **Packed/NDArray fast path** (machine reals): the LU factorisation runs through LAPACK `dgetrf` / `zgetrf` (in-house partial-pivot `double` LU as fallback).
- **Overflow/underflow → arbitrary precision**: the determinant is the product of the LU pivots, which can exceed the IEEE-double range even when the matrix is ordinary — a 200×200 `RandomReal[{-10,10}]` matrix has `|Det| ≈ 10^340`. Rather than return `inf` (or `0` from a mid-product underflow), the pivot product is re-accumulated in a 53-bit-mantissa MPFR value whose exponent range is effectively unbounded, so the answer is finite and correct to machine precision (`Det[RandomReal[{-10,10},{200,200}]]` returns a `≈ -1.08×10^340` real). A genuinely singular matrix still returns `0`.
- **Arbitrary-precision matrices**: a genuine MPFR matrix (precision > machine) uses an `O(n^3)` MPFR LU determinant (`mpfr_det_dispatch`) rather than the `O(n!)` Laplace expansion, which previously hung for `n ≳ 12`.
- A non-square machine-real matrix emits `Det::matsq` once and stays
  unevaluated. It used to run to `$RecursionLimit`: the buffer path declined,
  re-evaluated `Det` on the List, and the List was packed back onto the buffer
  path.

```mathematica
In[1]:= Det[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= 0

In[2]:= Det[{{1.7, 7.1, -2.7}, {2.2, 8.7, 3.2}, {3.2, -9.2, 1.2}}]
Out[2]= 251.572

In[3]:= Det[{{a, b, c}, {d, e, f}, {g, h, i}}]
Out[3]= -c e g + b f g + c d h - a f h - b d i + a e i
```

## Cross
Gives the vector cross product of its arguments.
- `Cross[v1, v2, ...]`

**Features**:
- `Protected`.
- Returns the cross product or totally antisymmetric product of $n-1$ vectors of length $n$.
- Works for symbolic and numerical inputs.
- Outputs `Cross::nonn1` if list arguments are not equal-length vectors or if the number of arguments is not one less than their length. A non-list argument (`Cross[a[x], b[x]]`) leaves the call unevaluated silently, as in the reference (v0.6134).

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

## Norm
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
- **Packed/NDArray fast path**: an inexact vector or matrix routes through the
  LAPACK-backed `ndla_norm` (`dlange`/`zlange`, and an SVD for the induced
  matrix 2-norm); `Norm` is on the packed-array `AWARE` list. A machine-precision
  vector (even a small unpacked `List`) is routed here too, so it does not hit
  the overflow described next; exact vectors keep their exact symbolic answer
  (`Norm[{1,2}] == Sqrt[5]`) and MPFR vectors stay on the symbolic path.
- **Overflow-safe 2-/p-norm**: the vector 2-norm uses LAPACK `dnrm2`'s
  running-scale recurrence (and the `p`-norm scales by the max component), so a
  well-scaled vector whose true norm is representable does not overflow to `inf`
  or underflow to `0` — `Norm[{1e200, 1e200}]` gives `1.414×10^200` (was `inf`),
  and `Normalize` of the same vector gives `{0.707…, 0.707…}` (was the zero
  vector).
- **`Compile[]`**: both `Norm[v]`/`Norm[m]` and the two-argument forms lower to a
  machine `Real`. Over a declared array argument the compiled body delegates to
  `ndla_norm` (the `V_NORM` opcode carries the literal `p`), covering
  `Norm[v, 1]`, `Norm[v, Infinity]`, `Norm[v, k]`, `Norm[v, 2.5]`, and the matrix
  `Norm[m, 1]` / `Norm[m, Infinity]` / `Norm[m, "Frobenius"]`. An induced matrix
  `p`-norm with `p ∉ {1, 2, Infinity}` has no LAPACK path and stays interpreted.
- **Inline vector**: a `Norm[{e1, …, en}, p]` over a `List` literal (rather than a
  declared array) compiles by expanding to the scalar arithmetic above
  (`Sqrt[Σ Abs[ei]^2]`, `Σ Abs[ei]`, `Max[Abs[ei]]`, `(Σ Abs[ei]^p)^(1/p)`), so a
  `Norm[{f[x], …}]` body auto-compiles inside `Plot`/`Table`/`NIntegrate`/`FindRoot`.

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

### Norm: the spectral norm of an exact, symbolic or MPFR matrix
`Norm[m]` (`Norm[m, 2]`) of an exact or symbolic matrix is `Sqrt` of the
largest eigenvalue of the smaller Gram matrix `m . m^H` / `m^H . m`, in the
reference's form: `Norm[{{4, 1}, {0, 5}}]` is `Sqrt[21 + Sqrt[41]]`,
`Norm[{{1, 2, 3}, {4, 5, 6}}]` is `Sqrt[(91 + Sqrt[8065])/2]`, a 3 x 3 gives
`Sqrt[Root[...]]`, and a symbolic matrix a `Max` over the `Sqrt` of every
eigenvalue (it has no order). An MPFR matrix answers its largest singular value
(one-sided Jacobi) at the input precision; both were left unevaluated.

## Normalize
Gives the normalized form of a scalar, vector, or tensor.
- `Normalize[v]`
- `Normalize[z]`
- `Normalize[expr, f]`

**Features**:
- `Protected` (not `Listable` — it acts on the whole vector, not element-wise).
- `Normalize[v]` is `v / Norm[v]` when `v` is a vector or tensor; the empty list `{}` and any all-zero input is returned unchanged.
- `Normalize[z]` for a scalar (possibly complex) `z` is `z / Abs[z]`, with `Normalize[0]` returning `0`.
- `Normalize[expr, f]` is `expr / f[expr]`, again with the zero short-circuit.
- The short-circuit is triggered only by a literal numeric zero (Integer, Real, BigInt, MPFR) in the evaluated norm — a symbolic norm that happens to be zero stays in the symbolic division so the user can see what they wrote.
- Wrong arity (`0` or `≥ 3` arguments) prints `Normalize::argt` and leaves the call unevaluated.

```mathematica
In[1]:= Normalize[{1, 5, 1}]
Out[1]= {1/3/Sqrt[3], 5/3/Sqrt[3], 1/3/Sqrt[3]}

In[2]:= Normalize[{3, 4}]
Out[2]= {3/5, 4/5}

In[3]:= Normalize[3 + 4 I]
Out[3]= 3/5 + 4/5*I

In[4]:= Normalize[{x, y}]
Out[4]= {x/Sqrt[Abs[x]^2 + Abs[y]^2], y/Sqrt[Abs[x]^2 + Abs[y]^2]}

In[5]:= Normalize[{x, y}, f]
Out[5]= {x/f[{x, y}], y/f[{x, y}]}

In[6]:= Normalize[{0, 0, 0}]
Out[6]= {0, 0, 0}
```

## Tr
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

## Diagonal
Gives the elements on a diagonal of a matrix.
- `Diagonal[m]`
- `Diagonal[m, k]`

**Features**:
- `Protected`.
- `Diagonal[m]` gives the leading diagonal `{m[[1,1]], m[[2,2]], ...}`, of length
  `Min[rows, cols]` — so it works for a non-square `m`.
- `Diagonal[m, k]` gives the `k`-th diagonal: `k > 0` above the leading diagonal
  (superdiagonals), `k < 0` below (subdiagonals). An out-of-range `k` gives `{}`.
- Elements are copied verbatim, so symbolic, exact, machine and complex entries
  all flow through and an exact-integer diagonal stays exact.
- Generalises to higher rank: the diagonal of a rank-`n` array is rank-`(n-1)`,
  since each `m[[i, i+k]]` is itself an `(n-1)`-tensor.
- Machine-precision matrices (a packed `List` or a visible `NDArray`) slice the
  diagonal directly off the flat buffer (`ndstruct_diagonal`, `src/ndstruct.c`),
  and the rank-2 case lowers inside `Compile[]` (`ND_FNS` / `A_NDFN`).

```mathematica
In[1]:= Diagonal[{{a, b, c}, {d, e, f}, {g, h, i}}]
Out[1]= {a, e, i}

In[2]:= Diagonal[{{a, b, c}, {d, e, f}, {g, h, i}}, 1]
Out[2]= {b, f}

In[3]:= Diagonal[{{a, b, c}, {d, e, f}, {g, h, i}}, -1]
Out[3]= {d, h}

In[4]:= Diagonal[{{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}}]
Out[4]= {1, 6, 11}

In[5]:= Diagonal[{{a, b, c, d}, {e, f, g, h}, {i, j, k, l}, {m, n, o, p}}, -2]
Out[5]= {i, n}
```

## Inverse
Gives the inverse of a square matrix.
- `Inverse[m]`
- `Inverse[m, Method -> "<name>"]`

**Features**:
- `Protected`.
- **Symbolic matrices** (v0.6125). An exact matrix of order 6 or less that is not all-rational uses the adjugate over `Det`'s memoised cofactor expansion by default; an inexact matrix keeps `DivisionFreeRowReduction`. Each entry stays a quotient `num/den`, as in the reference, so `Inverse[{{a + b, c}, {c, d + e}}][[1, 1]]` is `(d + e)/(-c^2 + a d + b d + a e + b e)` instead of a sum of two fractions. A 3×3 with linear entries went from 0.13 s and LeafCount 33,264 to 8 ms and 3,753 (the reference gives 1,433).
- Works on both symbolic and numerical matrices.
- For matrices with approximate real or complex numbers, the inverse is generated to the maximum possible precision given the input.
- Issues `Inverse::sing` warning and returns unevaluated if the matrix is singular.
- Issues `Inverse::matsq` warning and returns unevaluated if the argument is not a non-empty square matrix.
- **FLINT acceleration** (when built with FLINT): an all-integer/rational matrix is inverted exactly via `fmpq_mat_inv` in polynomial time (default/division-free method); the inverse is unique so the result matches the classical path. Singular or symbolic matrices fall through unchanged. Exposed directly as `` FLINT`Inverse `` (see the FLINT` context section in *Structural Manipulation*).
- Satisfies the relation `a . Inverse[a] == Inverse[a] . a == IdentityMatrix[n]`.
- Satisfies the relation `Inverse[a . b] == Inverse[b] . Inverse[a]`.
- Accepts an optional `Method -> "<name>"` argument that selects the inversion algorithm. Shares the same method-name grammar as `RowReduce` and `LinearSolve`.
  - `Method -> Automatic` or `Method -> "Automatic"` (default) — alias for `"DivisionFreeRowReduction"`.
  - `Method -> "DivisionFreeRowReduction"` — Bareiss-like fraction-free Gauss-Jordan elimination on the augmented matrix `[A | I]`. Best choice for exact integer / rational / symbolic input — never produces a denominator larger than necessary.
  - `Method -> "OneStepRowReduction"` — classical Gauss-Jordan on `[A | I]` with one division per pivot per row entry. Each entry is canonicalised via `Together` so symbolic cancellations are still detected. Fast on numeric matrices.
  - `Method -> "CofactorExpansion"` — adjugate / determinant formula `A^-1[i,j] = (-1)^(i+j) Det[M[j,i]] / Det[A]`, with `Det` computed via Laplace cofactor expansion. Time complexity is `O(n! n^2)`; intended for small `n` or closed-form symbolic inverses of small matrices. Emits `Inverse::sing` if `Det[A]` is structurally zero.
- Unknown method names emit `Inverse::method` and the call remains unevaluated.
- A singular or non-square machine-real matrix emits `Inverse::sing` /
  `Inverse::matsq` once and stays unevaluated; it no longer recurses to
  `$RecursionLimit` (see `Det`).

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

In[11]:= Inverse[{{1, 2}, {3, 4}}, Method -> "CofactorExpansion"]
Out[11]= {{-2, 1}, {3/2, -1/2}}

In[12]:= Inverse[{{a, b}, {c, d}}, Method -> "OneStepRowReduction"] // Together
Out[12]= {{d/(a d - b c), -b/(a d - b c)}, {-c/(a d - b c), a/(a d - b c)}}
```

> Implementation lives in `src/linalg/inv.c` (registered by
> `matinv_init`).  The algorithm — fraction-free Gauss-Jordan
> elimination on the augmented matrix `[A | I]` — is unchanged.

### Inverse: conditioning and arbitrary precision
- A machine matrix whose LAPACK 1-norm rcond is below `10^6 eps` (probed:
  `2.215*^-10` warns, `2.23*^-10` does not) gives `Inverse::luc` and the
  inverse; an exactly singular one `Inverse::sing`. Every call warns.
- An MPFR (real or complex) matrix is inverted by Gauss-Jordan at the input
  precision (`src/linalg/rref_numeric.c`); it used to drop to MachinePrecision
  because the augmented identity was written as machine `1.`/`0.`. Entries the
  elimination never touched stay exact `0`, as in the reference.
- One machine number makes the whole matrix machine (contagion), so a mixed
  exact/machine matrix takes the LAPACK path.

## PseudoInverse
Gives the Moore-Penrose pseudoinverse of a rectangular (or
rank-deficient) matrix.
- `PseudoInverse[m]`
- `PseudoInverse[m, Tolerance -> t]`

**Features**:
- `Protected`.
- Works on rectangular and rank-deficient matrices over the integers,
  rationals, machine-precision reals, MPFR reals, exact complex
  (`Complex[a, b]` entries), and inexact complex.
- For non-singular square matrices, `PseudoInverse[m] == Inverse[m]`.
- Computes a full-rank decomposition `m = B . C` (with `B` `m x r` and
  `C` `r x n`) by row-reducing `m` to identify the rank `r` and the
  pivot columns, then returns
  `PseudoInverse[m] = ConjugateTranspose[C] . Inverse[C . ConjugateTranspose[C]] . Inverse[ConjugateTranspose[B] . B] . ConjugateTranspose[B]`.
- For the `m x n` zero matrix, returns the `n x m` zero matrix.
- **A machine-real matrix takes the SVD instead.** `A⁺ = V Σ⁺ Uᵀ` from one thin
  LAPACK `gesdd`, with the singular values at or below the cutoff zeroed. The
  rationalise-and-row-reduce pipeline described above is the right answer for an
  exact or symbolic matrix and an unusable one for a machine matrix — it did not
  finish on a 300 × 300 in 180 s, where the SVD path takes 23 ms. An exact
  matrix is unaffected and still answers exactly.
- **An MPFR matrix takes the arbitrary-precision SVD** (v0.6143): `A⁺ = V Σ⁺ Uᴴ`,
  the rank decided by the precision. The automatic cutoff is
  `max(m, n) × 2^(1−bits) × σ₁`, the precision analogue of the machine rule. Rationalising
  an MPFR matrix instead turns entries known to 24 digits, such as
  `Sin[2 N[Pi, 24]] ≈ 1e-24`, into exact nonzero numbers, making a rank-deficient
  matrix full rank and its pseudoinverse meaningless (a 24-digit `Fit` with a
  vanishing basis column returned `-637534208 Sin[x/2]`). The exact pipeline remains
  the fallback when the SVD does not evaluate.
- A `Real` matrix small enough not to be held as a
  [packed list](packed-arrays.md) is rationalised at machine precision,
  computed exactly to preserve rank, and numericalised back.
- The `Tolerance` option accepts `Automatic` (default), a non-negative
  number, or a non-negative `Rational`. It takes effect on the SVD path, where
  there are singular values to truncate: `Automatic` is
  `max(m, n) × $MachineEpsilon × σ₁` — LAPACK's own rank criterion, and NumPy's
  default `rcond` — and an explicit `t` is a fraction of the largest singular
  value `σ₁`. The exact pipeline has no singular values and ignores it; the MPFR SVD path honours it as the machine one does.
- Issues `PseudoInverse::matrix` warning and returns unevaluated if the
  argument is not a non-empty rank-2 tensor.
- Returns unevaluated when an unknown option name is supplied or
  `Tolerance` receives a negative value.
- Satisfies the Moore-Penrose identities
  `m . p . m == m` and `p . m . p == p` for `p = PseudoInverse[m]`.

```mathematica
In[1]:= PseudoInverse[{{1,2},{3,4}}]
Out[1]= {{-2, 1}, {3/2, -1/2}}

In[2]:= PseudoInverse[{{1,2},{3,4}}] == Inverse[{{1,2},{3,4}}]
Out[2]= True

In[3]:= PseudoInverse[{{1,2,3},{4,5,6},{7,8,9}}]
Out[3]= {{-23/36, -1/6, 11/36}, {-1/18, 0, 1/18}, {19/36, 1/6, -7/36}}

In[4]:= PseudoInverse[{{0,0,0},{0,0,0}}]
Out[4]= {{0, 0}, {0, 0}, {0, 0}}

In[5]:= PseudoInverse[{{2,3},{2,2},{3,1},{4,3}}]
Out[5]= {{-29/134, -2/67, 22/67, 17/134}, {49/134, 8/67, -21/67, -1/134}}

In[6]:= PseudoInverse[{{1.25, 3.2, 3.2}, {7.9, -1.4, 5.1}, {0, 0, 0}}]
Out[6]= {{-0.0385185, 0.0966633, 0.}, {0.210183, -0.0659894, 0.}, {0.117363, 0.0282303, 0.}}
```

### PseudoInverse: machine complex
A machine complex matrix (and any small machine matrix, packed without the
packing threshold) uses the complex SVD (`zgesdd`, `P = V S^+ U^H` by
`zgemm`); it ran the rationalised exact pipeline and a 20 x 20 did not finish.

## MatrixPower
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
- **A machine-real matrix stays a [packed list](packed-arrays.md) throughout**,
  so each squaring is one `Dot` and reaches BLAS `dgemm`: `MatrixPower[A, 4]` on
  a 300 x 300 costs 0.69 ms, the two matrix products it is. An *integer* matrix
  is materialised first and stays exact, because the exact fourth power of an
  integer matrix is an integer matrix and a float64 buffer cannot hold one past
  2^53.
- **Non-integer exponents** (rational, real, complex or symbolic `n`) give the
  principal power, `exp(n log m)`, through the [MatrixFunction](#matrixfunction)
  machinery (`src/linalg/matfunc.c`): an exact matrix with an exact or symbolic
  `n` uses its Jordan form (`MatrixPower[{{2, 1}, {0, 2}}, s]` is
  `{{2^s, 2^(-1 + s) s}, {0, 2^s}}`); a machine matrix uses the complex Schur
  form -- the Bjorck--Hammarling square root for `n == 1/2`, Schur--Parlett with
  the closed-form derivatives of `z^n` otherwise; an MPFR matrix keeps its
  precision (Denman--Beavers square roots, `exp(n log m)`). An exact matrix with
  a machine `n` is computed in machine arithmetic. A real matrix whose power is
  real (no eigenvalue on the negative real axis) gets a real answer; a negative
  eigenvalue gives the principal complex power. A real exponent always takes
  this route, even an integer-valued one (`MatrixPower[m, 2.]`), as in the
  reference; the integer path is unchanged.
- `MatrixPower::zegev` for a negative power of a singular machine matrix;
  `MatrixPower::nosol` when a zero eigenvalue sits in a Jordan block
  (`MatrixPower[{{0., 1.}, {0., 0.}}, 0.5]`); the exact nilpotent case stays
  unevaluated without a message, as in the reference.
- Compile: `MatrixPower[m, s]` lowers for a real or complex rank-2 `m` and an
  integer OR real `s` (the A_NDFN trailing argument may be a Real register);
  when a real matrix's power is complex the compiled call declines to the
  interpreter, which gives the complex answer.

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

In[11]:= MatrixPower[{{2, 1}, {1, 2}}, 1/2]
Out[11]= {{1/2 + Sqrt[3]/2, -1/2 + Sqrt[3]/2}, {-1/2 + Sqrt[3]/2, 1/2 + Sqrt[3]/2}}

In[12]:= MatrixPower[{{1., 2.}, {3., 4.}}, 0.5]
Out[12]= {{0.553689 + 0.464394 I, 0.806961 - 0.212426 I}, {1.21044 - 0.31864 I, 1.76413 + 0.145754 I}}
```

### MatrixPower: negative powers of a singular matrix
`MatrixPower[m, -k]` says `MatrixPower::sing` (not `Inverse::sing`) when `m` is
singular, and for a machine matrix already when its 1-norm rcond is below
`100 eps` (probed: `2.22*^-14` is singular, `2.3*^-14` is not); between that and
`10^6 eps` it answers with `Inverse::luc`, as the reference does. An MPFR matrix
keeps its precision through the inverse.

## MatrixExp
The matrix exponential.
- `MatrixExp[m]`: the matrix exponential of the square matrix `m`.
- `MatrixExp[m, v]`: the matrix exponential of `m` applied to the vector `v`.

**Features**:
- `Protected`. Options `Method -> Automatic` (accepted; one algorithm per
  precision class) and `Parameter -> t` (computes `MatrixExp[t m]`).
- **SparseArray** `m` (v0.6188, `src/sparse/sparse_la_func.c`): a machine
  SparseArray answers a **SparseArray**, the reference's head, computed on the
  dense form within the 2^27-element densify cap. An exact, MPFR or symbolic
  one answers a List through the dense path. `MatrixExp[s, v]` for a machine
  SparseArray runs the Al-Mohy–Higham truncated Taylor action on SpMV and
  never forms the exponential, so it scales with the SpMV: the 90000 x 90000
  2-D Laplacian takes 76 ms. The answer is a dense vector, as in the reference.
- **Exact and symbolic** `m`: the closed form over the Jordan form
  `m = S J S^-1`, written as `sum f^(k)(lambda)/k! C_(lambda,k)` with each
  coefficient matrix simplified on its own, so the answer is a short sum of
  `E^lambda` terms. A nilpotent block gives the finite series
  (`MatrixExp[{{0, 1}, {0, 0}}]` is `{{1, 1}, {0, 1}}`). A real exact matrix with
  complex eigenvalues is answered in real form via `ComplexExpand`
  (`MatrixExp[{{0, 1}, {-1, 0}}]` is `{{Cos[1], Sin[1]}, {-Sin[1], Cos[1]}}`).
- **Machine** `m` (real or complex): Pade scaling and squaring with the
  Al-Mohy--Higham (2009) degree selection (3/5/7/9/13) and `ell(A, m)`
  backward-error correction against overscaling. A real matrix gives a real
  packed result; a 1 x 1 matrix is the scalar `Exp` exactly. An overflowing
  machine result is retried at 53-bit MPFR.
- **MPFR** `m`: scaling and squaring with a truncated Taylor series whose degree
  and scaling are chosen for the working precision (guard bits for the
  squarings), so the answer keeps the input's precision.
- `MatrixExp[m, v]` for machine input never forms `exp(m)`: the Al-Mohy--Higham
  (2011) truncated-Taylor action with the `theta_m` table and the trace shift
  (falling back to `MatrixExp[m] . v` when that is cheaper). `v` must be a
  non-empty vector (`MatrixExp::vector`) of length `Length[m]`
  (`MatrixExp::lslc`).
- `MatrixExp::matsq` for a List that is not a non-empty square matrix; a
  non-List argument stays unevaluated silently.
- Packed / NDArray input reaches the kernels directly (on `pack.c`'s `AWARE`
  list; `ndla_matrixexp`); List, packed and NDArray answers are bit-identical.
- Compile: `MatrixExp[m]` (ND_FNS, rank 2) and `MatrixExp[m, v]` (ND_FN2S,
  rank rule `R2_MATVEC`) lower for real and complex operands.
- SparseArray: answered on the dense form (gate mode `SPG_DENSE`); the reference
  re-sparsifies a machine answer (see the sparse gate).

```mathematica
In[1]:= MatrixExp[{{1, 1}, {0, 1}}]
Out[1]= {{E, E}, {0, E}}

In[2]:= MatrixExp[{{1., 2.}, {3., 4.}}]
Out[2]= {{51.969, 74.7366}, {112.105, 164.074}}

In[3]:= MatrixExp[{{1., 2.}, {3., 4.}}, {1., 1.}]
Out[3]= {126.706, 276.179}

In[4]:= MatrixExp[{{1, -2}, {2, 1}}]
Out[4]= {{E Cos[2], -E Sin[2]}, {E Sin[2], E Cos[2]}}
```

## MatrixLog
The principal matrix logarithm.
- `MatrixLog[m]`: the principal logarithm of the square matrix `m`.

**Features**:
- `Protected`. Option `Method -> Automatic`.
- **Exact and symbolic**: the Jordan-form closed form with `Log` and its
  derivatives (`MatrixLog[{{1, 1}, {0, 1}}]` is `{{0, 1}, {0, 0}}`;
  `MatrixLog[{{-1, 0}, {0, 2}}]` is `{{I Pi, 0}, {0, Log[2]}}`). A zero
  eigenvalue: `MatrixLog::fnand` (diagonalisable) or `MatrixLog::fnanc` (in a
  Jordan block), unevaluated.
- **Machine**: inverse scaling and squaring on the complex Schur form
  (Al-Mohy--Higham 2012): triangular square roots until `||T - I||` is within the
  Pade bound, the Gauss--Legendre partial-fraction Pade approximant of
  `log(I + X)`, then the diagonal and superdiagonal recomputed from closed forms
  (the unwinding-number formula). A real matrix goes through the real Schur form
  (triangularised by Givens rotations), so a negative real eigenvalue lands on
  `+I Pi` deterministically. A real matrix with no eigenvalue on the negative
  real axis gets a real answer (`MatrixLog[{{0., 1.}, {-1., 0.}}]` is real);
  a negative eigenvalue gives a complex answer. Singular: `MatrixLog::zegev`.
- **MPFR**: inverse scaling and squaring with Denman--Beavers square roots and
  the `atanh` series, at the input's precision; an eigenvalue on the branch cut
  is handled by rotating the matrix off the cut (`log A = i theta + log(e^(-i
  theta) A)`), located by a machine eigenvalue pass.
- Compile: lowers for real and complex rank-2 operands; a real matrix whose log
  is complex declines to the interpreter.

```mathematica
In[1]:= MatrixLog[{{2, 1}, {1, 2}}]
Out[1]= {{Log[3]/2, Log[3]/2}, {Log[3]/2, Log[3]/2}}

In[2]:= MatrixExp[MatrixLog[{{1., 2.}, {3., 4.}}]]
Out[2]= {{1., 2.}, {3., 4.}}
```

## MatrixFunction
A primary matrix function.
- `MatrixFunction[f, m]`: the matrix generated by the scalar function `f` (a
  symbol or pure function) at the square matrix `m`.

**Features**:
- `Protected`. Option `Method -> Automatic`.
- `Exp`, `Log` and (machine / MPFR) `Sqrt` use the `MatrixExp`, `MatrixLog` and
  square-root kernels, so `MatrixFunction[Exp, m] === MatrixExp[m]`.
- **Exact and symbolic**: the Jordan form with the derivatives of `f` (by `D`),
  so an undefined `f` keeps its derivatives:
  `MatrixFunction[f, {{1, 1}, {0, 1}}]` is `{{f[1], f'[1]}, {0, f[1]}}`, and a
  non-analytic one too (`Abs'[1]`). A singular value: `MatrixFunction::nosol`.
- **Machine**: Schur--Parlett (Davies--Higham 2003): eigenvalues clustered at
  distance 0.1, the Schur form reordered (adjacent Givens swaps) so clusters are
  contiguous, each cluster by a Taylor series about its mean with `f`'s
  derivatives from `D`, the off-diagonal blocks by the block Parlett recurrence.
  A derivative that is not numeric: `MatrixFunction::drvnnum`. An `f` whose
  values are not numbers at all (an undefined `g`) gives the
  eigen-decomposition form `V . DiagonalMatrix[g /@ vals] . Inverse[V]`.
- **MPFR**: through the exact engine on the rationalised input, numericalised
  back to the input precision (Exp / Log / Sqrt use the MPFR kernels).
- `MatrixFunction::argr` for a wrong argument count, `::matsq` for a bad matrix.
- Not compiled: the first argument is a function value, which no compiled
  register carries (`tools/compile_coverage.py` EXEMPT).

```mathematica
In[1]:= MatrixFunction[Sin, {{1, 1}, {0, 1}}]
Out[1]= {{Sin[1], Cos[1]}, {0, Sin[1]}}

In[2]:= MatrixFunction[Sin, {{1., 2.}, {3., 4.}}]
Out[2]= {{-0.465581, -0.148424}, {-0.222637, -0.688218}}
```

### MatrixFunction: normal matrices
When the Schur form of a machine matrix is diagonal to rounding (every
Hermitian, real symmetric or unitary matrix), `f` is applied to the eigenvalues
directly: Schur-Parlett's 0.1 clustering otherwise made a single Taylor-series
atom of a tightly clustered spectrum. `MatrixFunction[Sin, s]` of a 200 x 200
SPD matrix takes 0.05 s (was 0.87 s; scipy `funm` 0.066 s).
At arbitrary precision `MatrixFunction[f, m]` uses the Parlett recurrence on
the complex Schur form at twice the input precision (it rationalised and
worked symbolically: a 6 x 6 30-digit matrix ran past any time limit); a
repeated eigenvalue in a coupled block (a defective matrix) keeps the exact
path, so `MatrixFunction[Sin, N[{{1, 1}, {0, 1}}, 30]]` has `Cos[1]` above the
diagonal.

### MatrixFunction: clustered spectra (Higham–Davies truncation)
A non-normal machine matrix whose eigenvalues chain into one large
Schur–Parlett block (any random 100 x 100 matrix scaled small) is summed as a
Taylor series about the block's mean eigenvalue, and the series stops only when
the Higham–Davies tail bound `||M^(k+1)/(k+1)!|| mu max_j |f^(k+1+j)|/j!` falls
below `u ||F||` (`mu = ||(I - |N|)^-1 e||`, derivatives at the mean and the
farthest eigenvalue). Derivatives of `f` are built incrementally and cached with
no fixed order limit. Before v0.6221 the series required at least as many terms
as the block size and the derivative cache stopped at order 64, so such a matrix
came back unevaluated: `MatrixFunction[Cos, RandomReal[1, {100, 100}]/10]`
now takes 0.015 s (scipy `cosm` 0.006 s), with errors of order `1e-15` against
`(MatrixExp[I m] + MatrixExp[-I m])/2`.

## RowReduce
Gives the row-reduced form of the matrix `m`.
- `RowReduce[m]`
- `RowReduce[m, Method -> "<name>"]`
- `RowReduce[m, ZeroTest -> f]`

**Features**:
- `Protected`.
- Uses fraction-free division logic to perform exact algorithmic reduction across numerical, rational, and symbolics expressions natively avoiding division errors.
- **FLINT acceleration** (when built with FLINT): an all-integer/rational matrix is reduced exactly via `fmpq_mat_rref`; the reduced row echelon form is unique so it matches the classical result, computed in polynomial time. Symbolic matrices fall through. Exposed directly as `` FLINT`RowReduce `` (see the FLINT` context section in *Structural Manipulation*).
- Lives in `src/linalg/linsolve.c`; the helper primitives (Laplace cofactor determinant, exact polynomial division, tensor flatten / dimensions) are exposed from `src/linalg/util.c` and `src/linalg/det.c`.
- Accepts an optional `Method -> "<name>"` argument:
  - `Method -> Automatic` or `Method -> "Automatic"` (default) — alias for `"DivisionFreeRowReduction"`.
  - `Method -> "DivisionFreeRowReduction"` — Bareiss-like fraction-free Gauss-Jordan. Best for exact integer / rational / symbolic input — never produces a denominator larger than necessary.
    When every entry is a rational function over Q of symbols, with radicals allowed, the elimination runs natively on FLINT `fmpz_mpoly` (`flint_mpoly_mat_rref`, src/poly/flint_bridge.c, since v0.6108). The entries' field is first mapped onto a purely transcendental one:
    - a radical of a symbol, `s^(p/q)`, makes `s` stand for its own m-th root;
    - a radical of a compound base `B`, such as `Sqrt[(a - 1) b^2 + c]`, becomes a fresh generator `g`, with a variable `v` in which `B` is linear eliminated: `B = kappa v + B0` gives `v = (g^m - B0)/kappa`.

    Rows are cleared to `Z[t]`, the fraction-free Gauss-Jordan runs with exact `fmpz_mpoly` divisions, and the pivot rows are normalised in lowest terms and mapped back. Radicals of numbers (`Sqrt[2]`), numeric-constant symbols (`Pi`), and a compound base linear in no variable (`Sqrt[a^2 + b^2]`) are not free generators; those matrices take the expression-level path below. On a 17 x 6 parametric ansatz system that path did not finish; the native one takes 1 s.
    Its exact divisions (`exact_div_wrapper`, src/linalg/util.c, shared with `LinearSolve`, `NullSpace`, `MatrixRank` and `Inverse`) work over radical entries: when plain exact division fails, each radical base `B` with exponent denominator `m` is rewritten as `g^m` (the substitution `Together` and `Cancel` use), so `c`, `Sqrt[c]` and `c^(3/2)` become powers of one generator. The division then runs in a free polynomial ring and maps back; the quotient is checked in the original ring. Entries with `AlgebraicNumber` coefficients, or with a radical of a constant (`Sqrt[2]`), keep the old path: neither is a free generator (`g = Sqrt[2]` still satisfies `g^2 = 2`, which the free ring forgets). Before v0.6107 they were three independent variables, so a quotient such as `c^2/Sqrt[c]` was missed, and the unreduced fraction left in its place defeated the pivot zero tests. A consistent system over `Q(a, b, c)(Sqrt[c])` then reduced as inconsistent.
  - `Method -> "OneStepRowReduction"` — classical Gauss-Jordan with one division per pivot per element. Each entry is canonicalised via `Together` so symbolic cancellations are still detected. Fast on numeric matrices.
  - `Method -> "CofactorExpansion"` — for a non-singular square matrix, returns the identity (verified via `Det[m] != 0` computed by Laplace cofactor expansion). On singular or rectangular input, falls back to `"DivisionFreeRowReduction"` and emits `RowReduce::cofnsq`.
- `ZeroTest -> f` supplies a predicate deciding when an entry is zero during the reduction (zero iff `f[entry]` is `True`). Both a body such as `(RootReduce[Together[#]] === 0 &)` and a predicate head such as `PossibleZeroQ` are accepted. With a `ZeroTest`, an exact Gauss-Jordan RREF consults the predicate at every pivot (a structural RREF would pivot on an algebraic zero the predicate rejects) and keeps entries canonical via `Together`. Shares the machinery with `NullSpace` (`matsol_rref_with_zerotest`); `Method` and `ZeroTest` may be combined.
- Unknown method names emit `RowReduce::method` and the call remains unevaluated.
- The pivot 1s of a machine-real matrix are Reals, even where a pivot row
  divides to exactly 1: `RowReduce[{{1., 2., 1.}, {2., 4., 2.}}]` is
  `{{1., 2., 1.}, {0., 0., 0.}}`.

```mathematica
In[1]:= RowReduce[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= {{1, 0, -1}, {0, 1, 2}, {0, 0, 0}}

In[2]:= RowReduce[{{a, b, c}, {d, e, f}, {a+d, b+e, c+f}}]
Out[2]= {{1, 0, (-b f + c e)/(a e - b d)}, {0, 1, (-a f + c d)/(-a e + b d)}, {0, 0, 0}}

In[3]:= RowReduce[{{2, 1, 0}, {0, 3, 1}, {1, 0, 2}}, Method -> "CofactorExpansion"]
Out[3]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}

In[4]:= RowReduce[{{a, b}, {c, d}}, Method -> "OneStepRowReduction"]
Out[4]= {{1, 0}, {0, 1}}
```

### RowReduce of an inexact matrix
A machine or MPFR matrix (real or complex) is reduced by a tolerance-aware
Gauss-Jordan with partial pivoting (`src/linalg/rref_numeric.c`): a candidate
pivot counts as zero when `|p| <= max(r, c) eps max|a_ij|` for a machine matrix
-- the reference's own rule, probed on `diag(1, e)` of orders 2, 3 and 6 -- and
64 times that at arbitrary precision. The answer has the reference's form:
exact `1` at the pivots, exact `0` below them and in the zero rows, an inexact
`0.` above a machine pivot (exact `0` above an MPFR one). A machine pivot below
`10^6 eps max|a_ij|` warns `RowReduce::luc`. A singular 30-digit matrix used to
row-reduce to the identity (the `10^-30` residue taken as a pivot) or to
`Overflow[]` entries.

```mathematica
In[1]:= RowReduce[{{1., 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= {{1, 0., -1.}, {0, 1, 2.}, {0, 0, 0}}
```

## NullSpace
Gives a basis for the null space of a matrix.
- `NullSpace[m]`
- `NullSpace[m, Method -> "<name>"]`
- `NullSpace[m, ZeroTest -> f]`

**Features**:
- `Protected`.
- Returns a list of linearly-independent vectors whose span equals
  `{v : m . v == 0}`. If `m` has full column rank the result is the
  empty list `{}`.
- Works on numerical (Integer / Rational / Real / MPFR / Complex),
  big-integer, and symbolic matrices.
- The matrix `m` may be square or rectangular.
- Basis vectors are returned with the **rightmost free column first**,
  matching the standard ordering.
- For exact integer / rational input each basis vector is scaled to
  clear integer denominators, so the result is integer-valued whenever
  the input is integer-valued. For symbolic input the basis vectors
  are left in their natural rational form.
- Internally calls `RowReduce[m, Method -> "<name>"]` and extracts the
  basis from the resulting RREF. The `Method` option therefore shares
  the same grammar as `RowReduce` / `LinearSolve` / `Inverse`:
  - `Method -> Automatic` (default) — alias for
    `"DivisionFreeRowReduction"`.
  - `Method -> "DivisionFreeRowReduction"` — Bareiss-like fraction-free
    Gauss-Jordan; best for exact integer / rational / symbolic input.
  - `Method -> "OneStepRowReduction"` — classical Gauss-Jordan with one
    division per pivot per element. Each entry is canonicalised via
    `Together`.
  - `Method -> "CofactorExpansion"` — identity-if-invertible path
    inside RowReduce; falls back to `"DivisionFreeRowReduction"` on
    singular / rectangular input.
- `ZeroTest -> f` supplies a predicate that decides when an entry is
  zero during the reduction (an entry is zero iff `f[entry]` is `True`).
  Both a body such as `(RootReduce[Together[#]] === 0 &)` and a predicate
  head such as `PossibleZeroQ` are accepted. With a `ZeroTest`, `NullSpace`
  runs its own exact Gauss-Jordan RREF that consults the predicate at every
  pivot (a structural RREF would pivot on an algebraic zero the predicate
  rejects), keeping entries canonical via `Together`; e.g.
  `NullSpace[{{1, Sqrt[2]}, {Sqrt[2], 2}}, ZeroTest -> PossibleZeroQ]` is
  `{{-Sqrt[2], 1}}`. `Method` and `ZeroTest` may be combined.
- Issues `NullSpace::matrix` and returns unevaluated if the argument
  is not a non-empty rank-2 tensor.
- Issues `NullSpace::method` and returns unevaluated for unknown
  method names.

```mathematica
In[1]:= NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= {{1, -2, 1}}

In[2]:= NullSpace[{{a, b}, {2 a, 2 b}}]
Out[2]= {{-(b/a), 1}}

In[3]:= NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}]
Out[3]= {}

In[4]:= NullSpace[{{a, b, c}, {c, b, a}, {0, 0, 0}}]
Out[4]= {{1, -((a + c)/b), 1}}

In[5]:= NullSpace[{{3, 2, 2, 4}, {2, 3, -2, 7}, {3, 2, 5, 7}}]
Out[5]= {{12, -23, -5, 5}}

In[6]:= NullSpace[IdentityMatrix[5]]
Out[6]= {}

In[7]:= m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}; m . First[NullSpace[m]]
Out[7]= {0, 0, 0}
```

> Implementation lives in `src/linalg/nullspace.c` (registered by
> `matnull_init`).  The basis extraction calls back into RowReduce via
> the evaluator, so any future improvement to RowReduce automatically
> propagates to NullSpace.

### NullSpace at arbitrary precision
An MPFR matrix (and a machine complex one) uses the tolerance-aware elimination
and returns an orthonormal basis at the input precision
(`NullSpace[N[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, 30]]` is
`{{0.408..., -0.816..., 0.408...}}`, as in the reference); it answered `{}`.

## MatrixRank
Gives the rank of a matrix.
- `MatrixRank[m]`
- `MatrixRank[m, Method -> "<name>"]`
- `MatrixRank[m, Tolerance -> t]`
- `MatrixRank[m, Method -> "<name>", Tolerance -> t]`

**Features**:
- `Protected`.
- Returns a non-negative `Integer` equal to the number of linearly
  independent rows of `m` (equivalently, of linearly independent
  columns).
- Works on numerical (Integer / Rational / Real / MPFR / Complex),
  big-integer, and symbolic matrices, square or rectangular.
- **Two execution paths**:
  - *Exact path* (every leaf is exact, no `Tolerance`): an all-integer/rational
    matrix gets its rank directly from FLINT (`fmpq_mat_rref`) in polynomial
    time when built with FLINT (rank is basis-independent, so the value matches
    the classical count); otherwise routes
    through `RowReduce[m, Method -> "<name>"]` and counts the
    non-zero rows of the RREF, using `is_zero_poly` for structural
    zero. Honors the same `Method` grammar as NullSpace / RowReduce /
    LinearSolve / Inverse (`Automatic` / `"DivisionFreeRowReduction"`
    / `"OneStepRowReduction"` / `"CofactorExpansion"`). The FLINT kernel is
    also exposed directly as `` FLINT`MatrixRank `` (see *Structural
    Manipulation*).
  - *Numerical path* (any inexact leaf, or any explicit `Tolerance`):
    runs partial-pivot Gaussian forward-elimination over a portable
    `double`-complex kernel with tolerance-aware pivot selection: a
    column is skipped when its largest sub-pivot `|entry|` is `<= t`.
    `Method` does not affect this path.
    An inexact matrix that also holds free symbols gets its *generic*
    rank on this path: each symbol is replaced by a fixed complex value
    off every proper subvariety, rows and then columns are scaled to unit
    max-modulus, and the result is eliminated numerically
    (`MatrixRank[{{a, 2. a}, {1., 2.}}]` is `1`, `MatrixRank[{{a, 1.}, {1., a}}]`
    is `2`).
  - *Arbitrary-precision path* (every inexact leaf MPFR, default
    `Tolerance`): partial-pivot elimination at the matrix's own precision
    `p` bits, a pivot counting as zero at
    `64 max(rows, cols) 2^(1 - p) Max[|entries|]` (the numeric `RowReduce`
    rule). `MatrixRank[N[HilbertMatrix[40], 100]]` is `40` (condition number
    `~10^58`: far beyond machine precision, far within 100 digits); until
    v0.6226 it was rounded to doubles and answered `21`.
- **Tolerance** accepted forms:
  - `Tolerance -> Automatic` (default) — `max(rows, cols) *
    MachineEpsilon * Max[|entries|]` for inexact matrices; `0`
    otherwise (so integer / rational input agrees with the exact
    path).
  - `Tolerance -> 0` — no tolerance; even arbitrarily small entries
    count.
  - `Tolerance -> <non-negative number>` — Integer, Real, Rational,
    or a `Power[10, k]`-style expression resolved via `N[...]`.
- Issues `MatrixRank::matrix` and returns unevaluated if the
  argument is not a non-empty rank-2 tensor.
- Issues `MatrixRank::opt` and returns unevaluated for an unknown
  `Method` value, an unknown option key, or a negative / symbolic
  `Tolerance`.

```mathematica
In[1]:= MatrixRank[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= 2

In[2]:= MatrixRank[{{a, b, c}, {d, e, f}, {g, h, i}}]
Out[2]= 3

In[3]:= MatrixRank[{{0, 5, 2, 4, 4}, {2, 5, 0, 4, 0}, {5, 1, 5, 4, 5}}]
Out[3]= 3

In[4]:= MatrixRank[{{1.25, 3.2, 3.2}, {7.9, -1.4, 5.1}, {1.1, 2.5, -1.5}}]
Out[4]= 3

In[5]:= MatrixRank[{{a, b}, {2 a, 2 b}}]
Out[5]= 1

In[6]:= m = {{1, 1, 1}, {0, 10^-10, 0}, {0, 0, 10^-20}};
        MatrixRank[m]
Out[6]= 3

In[7]:= MatrixRank[N[m]]
Out[7]= 2

In[8]:= MatrixRank[N[m], Tolerance -> 0]
Out[8]= 3

In[9]:= MatrixRank[N[m], Tolerance -> 10^-8]
Out[9]= 1
```

> Implementation lives in `src/linalg/matrank.c` (registered by
> `matrank_init`). The exact path delegates to RowReduce, so any
> improvement to RowReduce / its `Method` kernels propagates to
> MatrixRank.  The numerical path is self-contained: a portable
> `cplx_t = {re, im}` struct stands in for `double _Complex` to keep
> the build strictly C99.

## LatticeReduce
Gives an LLL-reduced basis for a lattice.
- `LatticeReduce[m]` — `m` an `n × d` matrix whose rows are the lattice
  generators.

**Features**:
- `Protected`.
- Returns an `n × d` matrix whose rows form a reduced basis of the same
  lattice (same `Z`- / `Z[i]`-module).
- Entries may be integers, rationals, Gaussian integers, or Gaussian
  rationals (`Complex[a, b]` with exact `a`, `b`).
- Exact Lenstra–Lenstra–Lovász reduction (Lovász parameter `δ = 3/4`)
  carried out entirely in GMP rational arithmetic — no floating point —
  so it is correct for both machine-size and arbitrary-precision (bignum)
  inputs. This matters for integer-relation finding, where a rounding
  error would yield a wrong relation.
- The Gram–Schmidt orthogonalisation uses the Hermitian inner product
  `⟨x, y⟩ = Σ x_k conj(y_k)`, so real and Gaussian lattices share one
  code path; size reduction rounds to the nearest Gaussian integer.
- The lattice — and hence `Abs[Det]` and every linear relation in the
  right null space — is preserved exactly.
- The rows must be linearly independent; a dependent generating set is
  reported with `LatticeReduce::dep`.
- Diagnostics: `LatticeReduce::argx` (wrong argument count),
  `LatticeReduce::matrix` (not a non-empty rectangular matrix),
  `LatticeReduce::latm` (an entry is not rational).

```mathematica
In[1]:= LatticeReduce[{{1, 0, 0, 1345}, {0, 1, 0, 35}, {0, 0, 1, 154}}]
Out[1]= {{0, 9, -2, 7}, {1, 1, -9, -6}, {1, -3, -8, 8}}

In[2]:= {w1, w2} = LatticeReduce[{{12, 2}, {13, 4}}]
Out[2]= {{1, 2}, {9, -4}}

(* Integer-relation finding: x0 + a1 x1 + a2 x2 == 0 for a1=2, a2=3 *)
In[3]:= a = {{1, 0, 0, -1}, {0, 1, 0, -2}, {0, 0, 1, -3}};
        b = LatticeReduce[a]
Out[3]= {{1, 0, 0, -1}, {-1, 1, 0, -1}, {-1, -1, 1, 0}}

In[4]:= b . {1, 2, 3, 1}    (* relations preserved *)
Out[4]= {0, 0, 0}

In[5]:= LatticeReduce[{{1, 2}, {3, 4.5}}]
LatticeReduce::latm: Matrix contains an entry that is not rational.
Out[5]= LatticeReduce[{{1, 2}, {3, 4.5}}]
```

> Implementation lives in `src/linalg/latticereduce.c` (registered by
> `linalg_init`). Each Gaussian-rational scalar is a pair of GMP `mpq_t`
> (`GRat`); the Gram–Schmidt data (`μ`, `|b*|²`) is computed once and
> maintained incrementally — updated in place on size reduction and on a
> Lovász swap via the conjugate-aware Cohen swap formulas — so no full
> recomputation is needed.

## HermiteDecomposition
Hermite normal form of an integer matrix, with the unimodular transform.
- `HermiteDecomposition[m]` — gives `{u, r}` where `u` is unimodular
  (`Abs[Det[u]] == 1`), `r` is the row Hermite normal form of `m`, and
  `u . m == r`. `r` is in echelon shape: pivots move strictly right as rows
  descend, every pivot is positive, and each entry *above* a pivot is reduced
  into `[0, pivot)`. Works for any `m × n` integer matrix (`m` need not equal
  `n`).
- Integer entries only; a non-integer entry leaves the call unevaluated with
  `HermiteDecomposition::intm`. Diagnostics: `HermiteDecomposition::argx`
  (wrong argument count), `HermiteDecomposition::matrix` (not a non-empty
  rectangular matrix).

```
In[1]:= HermiteDecomposition[{{2, 3, 1}, {4, 1, 5}, {6, 2, 0}}]
Out[1]= {{{0, 2, -1}, {1, 4, -3}, {1, 7, -5}}, {{2, 0, 10}, {0, 1, 21}, {0, 0, 36}}}

In[2]:= h = HermiteDecomposition[{{1, 2, 3}, {4, 5, 6}}]; h[[1]] . {{1, 2, 3}, {4, 5, 6}} == h[[2]]
Out[2]= True
```

> Implementation lives in `src/linalg/hnf.c` (registered by `linalg_init`).
> The row reduction eliminates each column with the extended-gcd `2×2`
> unimodular transform, recording every operation on `u`; pivots are then made
> positive and entries above them reduced. The same `linalg_hnf` primitive
> backs exact linear Diophantine **system** solving in `Solve[…, Integers]`
> (`src/solve/solveint_linear.c`).

## FindIntegerNullVector
Finds an integer relation among a list of numbers (PSLQ).
- `FindIntegerNullVector[{x1, …, xn}]` — integers `{a1, …, an}`, not all
  zero, with `a1 x1 + … + an xn == 0`.
- `FindIntegerNullVector[{x1, …, xn}, d]` — restricts the search to
  relations of norm `≤ d`.

**Features**:
- `Protected`.
- The `xi` may be **real or complex**, **exact or inexact**. For complex
  `xi` the `ai` are **Gaussian integers**.
- Built on the exact LLL machinery of `LatticeReduce`: the numbers are
  numericalised to a working precision `b` and embedded as the rows
  `(e_i | round(2^b x_i))` of an integer-relation lattice; the relation is
  read off the shortest reduced row.
- **Validation.** For exact `xi` the residual `a·x` is checked with
  `PossibleZeroQ`; a confidence guard (`n·log2(‖a‖²) < 1.35·b`) rejects
  large-coefficient artefacts that hold to only a few digits, forcing a
  precision increase. For inexact `xi` the relation holds to the precision
  of the input.
- **Precision** (`WorkingPrecision` option, default `Automatic`): inexact
  input uses the precision of the input; exact input starts at
  `$MachinePrecision` and escalates up to
  `$MachinePrecision + $MaxExtraPrecision` digits. An explicit digit count
  fixes the precision (no escalation).
- **Certified bound.** The rigorous LLL Gram–Schmidt bound
  `λ₁(L)² ≥ minᵢ |b*ᵢ|²` gives a lower bound
  `B = √(M² / (1 + n/4))` on the norm of *any* integer relation (`M² =
  minᵢ |b*ᵢ|²`), conditional on the numeric evaluation being correct to
  precision. `B` drives the no-relation diagnostics. Because this is the
  LLL bound (not the tighter PSLQ bound the reference reports), it is
  conservative: near the minimal-relation norm Mathilda may return the
  relation with a `lgrelb` message where the reference proves nonexistence.
- Relations are returned up to sign and (for complex input) up to a
  Gaussian-unit multiple — `{a}`, `{-a}`, and `{I a}` are all valid null
  vectors.
- Diagnostics: `FindIntegerNullVector::norel` (proven no relation with
  norm `≤ d`), `::lgrelb` (a relation was found but it exceeds `d`, and
  nonexistence is proven only below a smaller bound — the larger relation
  is returned), `::rnfb` (inexact/bounded: none found `≤ d`, proven none
  below a smaller bound), `::rnfu` (exact/unbounded: no relation found
  within the precision cap), `::ztest1` (the residual could not be
  proven zero and is assumed zero). When no vector is returned the call is
  left unevaluated.

```mathematica
In[1]:= FindIntegerNullVector[{Log[2], Log[4]}]
Out[1]= {-2, 1}

In[2]:= FindIntegerNullVector[{Pi, ArcTan[1/5], ArcTan[1/239]}]
Out[2]= {1, -16, 4}                 (* Machin's formula *)

In[3]:= a = Sqrt[2] + 3^(1/3); FindIntegerNullVector[a^Range[0, 6]]
Out[3]= {1, -36, 12, -6, -6, 0, 1}  (* minimal polynomial coefficients *)

In[4]:= FindIntegerNullVector[{1, 2 I + Sqrt[3], (2 I + Sqrt[3])^2}]
Out[4]= {-7, -4 I, 1}               (* Gaussian-integer relation *)

In[5]:= FindIntegerNullVector[{E, Pi}, 1000000]
FindIntegerNullVector::norel: There is no integer null vector for {E, Pi}.
Out[5]= FindIntegerNullVector[{E, Pi}, 1000000]
```

> Implementation lives in `src/linalg/latticereduce.c` alongside
> `LatticeReduce`, reusing its exact Gaussian-rational LLL kernel
> (`lll_reduce`, extended to report `minᵢ |b*ᵢ|²`).

## JordanDecomposition
Gives the Jordan decomposition of a square matrix.
- `JordanDecomposition[m]` — returns `{s, j}` where `j` is the Jordan
  canonical form of `m` and `s` is the similarity matrix, so that
  `m == s . j . Inverse[s]`. `j` is block-diagonal in Jordan blocks (an
  eigenvalue on the diagonal, `1`s on the superdiagonal); the columns of
  `s` are the generalized eigenvectors grouped into Jordan chains that
  match `j`'s blocks. A `1` on `j`'s superdiagonal marks a deficient
  eigenspace.

**Features**:
- `Protected`.
- Works on every input family:
  - exact integer / rational matrices (exact `s` and `j`);
  - free-symbolic matrices with closed-form eigenvalues (e.g. any
    2×2, or small matrices with a factorable characteristic
    polynomial);
  - complex matrices;
  - machine-precision Real matrices (the numeric eigensolver; a generic
    distinct spectrum is diagonalizable, so `j` is diagonal);
  - arbitrary-precision MPFR matrices (output at the input precision).
- Algorithm:
  - **Exact / symbolic.** The eigenvalues come from the characteristic
    polynomial (`eigen_char_poly_faddeev` → `Solve`); for each eigenvalue
    `λ` the Jordan block structure is read from the nullity sequence
    `dim ker (m - λ I)^k`, and a top-down chain-top selection (via a
    `MatrixRank`-based "extend a spanning set to a basis") builds the
    Jordan chains bottom-up. The generalized eigenvectors of `s` are any
    valid Jordan basis, not a canonical one.
  - **Numeric.** A matrix with a distinct spectrum is diagonalizable, so
    `s` is the numeric eigenvectors as columns and `j = DiagonalMatrix`
    of the eigenvalues. A **machine-real** matrix (a packed / `NDArray`
    argument is read straight off its float64 buffer, no delist) uses a
    LAPACK-`dgeev` kernel that solves once and builds `s`/`j` directly
    from the raw buffers — no intermediate boxed eigenvector list, no
    second solve (a 100×100 in ~9.6 ms). A complex-entry matrix (or
    `USE_LAPACK=0`) instead runs one `Eigenvectors` solve and recovers
    each eigenvalue from its eigenvector by an `O(n)` component ratio;
    arbitrary-precision input keeps the full-precision solves. A matrix
    with a repeated numeric eigenvalue (defective, or ambiguous at
    machine precision) is rationalised, decomposed by the exact core, and
    numericalised back, which recovers the genuine block structure.
- Result-fidelity notes: the exact printed `s` is not canonical (any
  Jordan basis satisfying `m == s . j . Inverse[s]` is valid); numeric
  block-splitting is not bit-reproducible from the reference (an ill-posed
  problem), so the numeric contract is the residual `m ≈ s . j . Inverse[s]`
  plus a block structure matching the numeric nullities.
- Limitation: an exact matrix with an *irrational, defective* eigenvalue
  whose generalized eigenspace cannot be spanned exactly (the
  `is_zero_poly` pivot test cannot decide the row reduction) is left
  unevaluated rather than returned wrong.
- Issues `JordanDecomposition::argx` for the wrong argument count, and
  `JordanDecomposition::matsq` for a non-square or empty matrix; both
  leave the call unevaluated.
- Not lowered by `Compile[]` (the result is a heterogeneous pair of
  matrices, like `QRDecomposition` / `SingularValueDecomposition`); a
  packed / `NDArray[…]` argument is accepted (delisted and re-evaluated).

```mathematica
In[1]:= JordanDecomposition[{{27, 48, 81}, {-6, 0, 0}, {1, 0, 3}}][[2]]
Out[1]= {{6, 0, 0}, {0, 12, 1}, {0, 0, 12}}

In[2]:= {s, j} = JordanDecomposition[{{-103, -191, -255}, {110, 190, 222}, {9, 9, 33}}];
        j
Out[2]= {{24, 0, 0}, {0, 48, 1}, {0, 0, 48}}

In[3]:= {s, j} = JordanDecomposition[{{-1.2, 2.7, 3.8}, {4.2, 4.4, 5.3}, {3.5, 7.6, 6.8}}];
        Diagonal[j]
Out[3]= {13.7715, -2.12527, -1.6462}
```

## SchurDecomposition
Gives the Schur decomposition of a numerical square matrix.
- `SchurDecomposition[m]` — returns `{q, t}` where `q` is orthonormal
  (unitary) and `t` is block upper-triangular (the Schur form), so that
  `m == q . t . ConjugateTranspose[q]`.
- `SchurDecomposition[{m, a}]` — the generalized (QZ) decomposition,
  `{q, s, p, t}` with `q, p` orthonormal and `s, t` upper-triangular, so that
  `m == q . s . ConjugateTranspose[p]` and `a == q . t . ConjugateTranspose[p]`.

**Features**:
- `Protected`.
- **Numerical only** — a generic matrix has no closed-form symbolic Schur
  form, so a non-numeric (symbolic) matrix is left unevaluated. Input families:
  - machine-precision Real / complex matrices, and `NDArray` / packed arrays
    (read straight off the buffer);
  - exact integer / rational matrices (computed at machine precision, since
    there is no exact Schur form);
  - arbitrary-precision MPFR matrices — supported for the **standard real**
    form (output at the input precision).
- Options:
  - `Pivoting -> True` returns `{q, t, d}` with an extra scaling / permutation
    matrix `d`, such that `m . d == d . q . t . ConjugateTranspose[q]` (the
    balancing transform, applies to the single-matrix form).
  - `RealBlockDiagonalForm -> True` (default) keeps `t` real with 2×2 blocks
    for complex-conjugate eigenvalue pairs (so `t` is block — not strictly —
    upper-triangular); `-> False` makes `t` complex upper-triangular.
  - `TargetStructure -> "Dense"` (default) `| "Structured"` — both return
    dense matrices (Mathilda has no distinct structured-matrix type; `t` is
    already triangular and `q` orthonormal, so nothing is lost).
- Algorithm:
  - **Machine.** LAPACK `dgees` (real) / `zgees` (complex) for the standard
    form and `dgges` / `zgges` (the QZ algorithm) for the generalized form.
    `RealBlockDiagonalForm -> False` on a real matrix loads it complex and
    uses the `z` driver. `Pivoting -> True` balances with `dgebal` and
    reconstructs `d` by back-transforming the identity with `dgebak`. A
    100×100 real matrix takes a few milliseconds.
  - **Arbitrary precision.** The standard real form reuses the in-house
    Hessenberg + Francis QR (`eigen_schur_real_mpfr`), which accumulates the
    orthogonal Schur vectors alongside the quasi-triangular Schur form.
- Result-fidelity notes: the Schur form is not unique (`q`, `t` are fixed only
  up to the eigenvalue ordering and unitary freedom within degenerate blocks),
  so the contract is the residual `m ≈ q . t . ConjugateTranspose[q]` with
  `q` orthonormal and `t` (quasi-)triangular, not bit-reproducible factors.
- Limitations: at arbitrary precision only the standard real form is native;
  a complex or generalized arbitrary-precision matrix falls back to machine
  precision (a complex MPFR QR / an MPFR QZ would be needed for full
  precision).
- Issues `SchurDecomposition::argx` for zero arguments and
  `SchurDecomposition::sqma` for an argument that is neither a non-empty
  square matrix nor a pair of square matrices; both leave the call
  unevaluated.
- Not lowered by `Compile[]` (the result is a heterogeneous tuple of
  matrices, like `QRDecomposition` / `JordanDecomposition`); a packed /
  `NDArray[…]` argument is accepted (on `pack.c`'s AWARE list).

```mathematica
In[1]:= {q, t} = SchurDecomposition[{{2.7, 4.8, 8.1}, {-.6, 0, 0}, {.1, 0, .3}}];
        Chop[{{2.7, 4.8, 8.1}, {-.6, 0, 0}, {.1, 0, .3}} - q . t . ConjugateTranspose[q]]
Out[1]= {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}

In[2]:= {q, t} = SchurDecomposition[{{1.81066, 0.31066, 1.5}, {-0.53033, 2.03033, 0.43934},
          {-0.96967, -0.53033, 2.56066}}];
        {UpperTriangularMatrixQ[t, -1], UpperTriangularMatrixQ[t]}
Out[2]= {True, False}

In[3]:= {q, s, p, t} = SchurDecomposition[{{{.5, 1}, {1.5, 2}}, {{2.5, 3}, {3.5, 4}}}];
        Chop[{{.5, 1}, {1.5, 2}} - q . s . ConjugateTranspose[p]]
Out[3]= {{0, 0}, {0, 0}}
```

### SchurDecomposition: exact input and complex arbitrary precision
- The reference treats Schur as a numerical decomposition: an exact numeric
  matrix (integers, rationals, Gaussian rationals, `Pi`) gives
  `SchurDecomposition::schurf` and a symbolic one `::schurn`, both unevaluated;
  a generalized pair with no inexact entry, or a symbolic one, gives `::exnum`;
  a non-square matrix `::matsq`. One inexact entry makes the matrix numeric.
- A complex MPFR matrix, or `RealBlockDiagonalForm -> False` at arbitrary
  precision, uses a complex Schur kernel (Givens Hessenberg reduction and the
  Wilkinson-shift QR iteration, 32 guard bits) and keeps the input precision.

## CholeskyDecomposition
Gives the Cholesky decomposition of a Hermitian (or real symmetric)
positive-definite matrix.
- `CholeskyDecomposition[m]` — the upper-triangular `c` with
  `ConjugateTranspose[c] . c == m`.

**Features**:
- `Protected` (the reference also sets `ReadProtected`; Mathilda never does).
- Every precision class:
  - **exact and symbolic** entries give an exact answer built with `Sqrt` and
    `Conjugate`, in the reference's form:
    `c[[i, i]] = Sqrt[m[[i, i]] - Sum c[[k, i]] Conjugate[c[[k, i]]]]` and
    `c[[i, j]] = (m[[i, j]] - Sum Conjugate[c[[k, i]]] c[[k, j]]) / c[[i, i]]`;
  - **machine** reals and complexes (a `List`, a packed array or an `NDArray`)
    use LAPACK `dpotrf` / `zpotrf` on the upper triangle (an in-house column
    Cholesky when LAPACK is absent), and answer a packed matrix with exact-zero
    bits below the diagonal. A complex factor is a packed complex array too, as
    the reference returns it (every machine linear-algebra result goes through
    `na_finish`, which boxed complex results into `Complex[...]` nodes until
    v0.6221);
  - **arbitrary precision** keeps the input precision, with exact `0` below the
    diagonal, as the reference does.
- Hermitian test: exact/symbolic `m` must equal `ConjugateTranspose[m]`, or equal
  `Transpose[m]` with no explicit complex entry (so `{{a, b}, {b, c}}` is
  accepted and `{{2, I}, {I, 2}}` is not); inexact `m` passes a relative
  `1e-12` test.
- Messages (reference text): `::matsq` (not a non-empty square matrix),
  `::herm` (not Hermitian or real symmetric), `::npdef` (exact, not positive
  definite), `::posdef` (inexact, not sufficiently positive definite),
  `::optx` (an option other than `TargetStructure`). A symbolic matrix with a
  non-positive numeric pivot, a vector and an atom are left unevaluated
  silently, as in the reference.
- Option `TargetStructure -> Automatic | "Dense" | "Structured"`: all give a
  dense `List`.
- `SparseArray` input is densified by the transparency gate (answer: a `List`,
  as in the reference).
- **Compile**: lowers for `_Real` and `_Complex` rank-2 operands (an `ND_FNS`
  row calling the same machine kernel, so `List`, packed, `NDArray` and
  `Compile[]` agree bit for bit); a non-positive-definite operand degrades to
  the interpreter, which reports it once.

```mathematica
In[1]:= CholeskyDecomposition[{{4, 2, 1}, {2, 3, 0}, {1, 0, 5}}]
Out[1]= {{2, 1, 1/2}, {0, Sqrt[2], -1/2*1/Sqrt[2]}, {0, 0, Sqrt[37/2]/2}}

In[2]:= CholeskyDecomposition[{{a, b}, {b, c}}]
Out[2]= {{Sqrt[a], b/Sqrt[a]}, {0, Sqrt[c - (b*Conjugate[b/Sqrt[a]])/Sqrt[a]]}}

In[3]:= CholeskyDecomposition[{{1, 2}, {2, 1}}]
CholeskyDecomposition::npdef: The matrix {{1, 2}, {2, 1}} is not positive definite.
Out[3]= CholeskyDecomposition[{{1, 2}, {2, 1}}]
```

## SingularValueList
Gives the singular values of a matrix.
- `SingularValueList[m]` — the nonzero singular values, in decreasing order.
- `SingularValueList[m, k]` — the first `k` (`k < 0`: the last `|k|`) of all
  `min(rows, cols)` values, zeros included; also `{k}`, `{i, j}`, `{i, j, s}`
  and `UpTo[k]` (up to `k` of all the values, zeros included, as for `k`).
- `SingularValueList[{m, a}]` — the generalized singular values (machine
  matrices only).
- A machine `SparseArray` with a count, and `Method -> "Arnoldi"`: see
  "Eigenvalues / Eigenvectors / Eigensystem / SingularValueList /
  SingularValueDecomposition: SparseArray and Method -> \"Arnoldi\"" after `Eigensystem`.

**Features**:
- `Protected`. Options `Tolerance -> Automatic`, `Method -> Automatic`
  (`"Arnoldi"` accepted; dense input runs the direct method), `Cubics`,
  `Quartics` (passed to `Eigenvalues` for exact input).
- Tolerance (reference semantics): a value `s` is dropped when
  `s < tol * s_max`. `Automatic` is `max(rows, cols) * $MachineEpsilon` for
  machine input, about `8 max(rows, cols) 2^-bits` at arbitrary precision, and
  exact zeros only for exact input. An all-zero matrix keeps its zeros
  (`{0, 0}`), as the reference does. With an explicit `k` every value is a
  candidate (`SingularValueList[{{1, 2}, {2, 4}}, 2]` is `{5, 0}`).
- Precision classes:
  - **machine** (List, packed, `NDArray`): LAPACK `gesdd` with `jobz = 'N'` —
    values only, no `U` or `V`; a packed real vector;
  - **arbitrary precision** (`svd_mpfr_values` in `svdecomp_mpfr.c`): a real
    matrix is reduced to upper bidiagonal form by Householder reflectors from
    both sides (fused `arb_approx_dot` inner products, rows / columns split
    over threads), and the bidiagonal's singular values come from the
    implicit-shift Golub-Kahan QR iteration.  The work runs at the input
    precision plus 32 guard bits; when the surviving values spread over more
    than 16 bits (sigma_max / sigma_min) the pass is repeated once with the
    spread added, so each value -- the smallest included -- carries the input
    precision in relative terms (the 1.2e-42 smallest value of the 50 x 50
    equispaced Vandermonde matrix at 100 digits: all 100).  A complex matrix
    (or a build without FLINT) uses one-sided Jacobi with `V` not accumulated;
  - **exact / symbolic**: `Sqrt` of the `Eigenvalues` of the smaller Gram
    matrix (`ConjugateTranspose[m] . m` or `m . ConjugateTranspose[m]`), each
    eigenvalue passed through `Together` when numeric, sorted decreasing.
- Generalized: LAPACK `ggsvd3` without vectors; the finite values
  `alpha / beta`, decreasing. Exact or MPFR pairs give `::gfargs`, as in the
  reference.
- Messages: `::matrix`, `::take` ("Cannot take singular values 1 through 3 out
  of the total of 2 singular values."), `::seqs`, `::tol`, `::emeth`,
  `::gfargs`, `::optx`.
- **Compile**: lowers for `_Real` and `_Complex` rank-2 operands (an `ND_FNS`
  row with the new `real_result` flag — a complex matrix's singular values are
  real), calling the interpreter's own `gesdd('N')` kernel, so all surfaces agree
  bit for bit. The option and `k` forms are not lowered.

```mathematica
In[1]:= SingularValueList[{{1, 2}, {3, 4}}]
Out[1]= {Sqrt[15 + Sqrt[221]], Sqrt[15 - Sqrt[221]]}

In[2]:= SingularValueList[{{1., 0.}, {0., 1.*^-3}}, Tolerance -> 0.01]
Out[2]= {1.}

In[3]:= SingularValueList[{{1, 2}, {2, 4}}, 2]
Out[3]= {5, 0}
```

## Minors
Gives the minors of a matrix.
- `Minors[m]` — the minors of order `n - 1`, `n = min(rows, cols)`.
- `Minors[m, k]` — the `k × k` minors: a `C(rows, k) × C(cols, k)` matrix
  whose `(r, c)` entry is the determinant of the rows in the `r`-th and the
  columns in the `c`-th `k`-subset, subsets in **lexicographic** order (the
  reference's convention: `Minors[{{a,b,c},{d,e,f},{g,h,i}}][[1]]` is
  `{a e - b d, a f - c d, b f - c e}`). `k = 0` gives `{{1}}`; `k < 0` or
  `k > n` gives `{}`; a non-integer `k` leaves the call unevaluated silently.
- `Minors[m, k, f]` — `f` applied to each `k × k` submatrix instead of `Det`.

**Features**:
- `Protected`. `::matrix` for an argument that is not a non-empty rectangular
  matrix.
- Machine matrices (List, packed, `NDArray`): one partial-pivot LU
  determinant per submatrix, a packed answer. Exact, symbolic and
  arbitrary-precision matrices: `Det` of each submatrix (FLINT for rational
  entries), at the input precision.
- **Compile**: `Minors[m]` lowers for `_Real` / `_Complex` rank-2 operands
  (rank 2 → rank 2, the same kernel); `Minors[m, k]` and `Minors[m, k, f]` are
  not lowered (they run in the interpreter).

```mathematica
In[1]:= Minors[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}]
Out[1]= {{-3, -6, -3}, {-6, -11, -4}, {-3, -2, 2}}

In[2]:= Minors[{{a, b, c}, {d, e, f}}, 2]
Out[2]= {{-(b d) + a e, -(c d) + a f, -(c e) + b f}}
```

### Minors: Det::luc
A machine minor whose infinity-norm rcond is below `10^6 eps` warns `Det::luc`,
as the reference's `Det` of each submatrix does.

## Adjugate
Gives the adjugate (classical adjoint) of a square matrix: the transposed
cofactor matrix, with `Adjugate[m] . m == Det[m] IdentityMatrix[n]`, defined
for a singular `m` too.

**Features**:
- `Protected` (the reference also sets `ReadProtected`; Mathilda never does).
- Machine matrices, well conditioned (LU `rcond >= sqrt(u)`): `det(m) m^-1`
  from one LU factorisation (`getrf` + `getri`), 0.0017 s at 200 x 200 (was
  0.0185 s through the SVD; numpy `det * inv` 0.0026 s).
- Machine matrices, near singular: the SVD form
  `adj(U S V^H) = det(U) det(V^H) V adj(S) U^H` with `adj(S)_k = Prod_{j≠k} s_j`
  (LAPACK `gesdd`), which stays accurate for a singular matrix where
  `Det[m] Inverse[m]` cannot be formed. One kernel for List, packed, `NDArray`
  and `Compile[]`.
- Exact rational, nonsingular, order ≥ 4: `Det[m] Inverse[m]` (FLINT).
  Otherwise (exact, symbolic, MPFR): cofactors through `Det`, at the input
  precision. `Adjugate[{{5}}]` is `{{1}}`.
- `::matsq` for a non-square argument; a second argument gives `::intnm`, as
  in the reference.
- **Compile**: lowers for `_Real` / `_Complex` rank-2 operands.

```mathematica
In[1]:= Adjugate[{{a, b}, {c, d}}]
Out[1]= {{d, -b}, {-c, a}}

In[2]:= Adjugate[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[2]= {{-3, 6, -3}, {6, -12, 6}, {-3, 6, -3}}
```

## Orthogonalize
Gram–Schmidt orthonormalisation.
- `Orthogonalize[{v1, v2, ...}]` — an orthonormal list under the inner product
  `u . Conjugate[v]`; a vector dependent on the earlier ones becomes a zero
  vector, so there are as many outputs as inputs.
- `Orthogonalize[{e1, e2, ...}, f]` — under the inner product `f[u, v]`; the
  `ei` need not be vectors (`Orthogonalize[{1, x, x^2}, Integrate[#1 #2, {x, -1, 1}] &]`
  gives normalised Legendre polynomials).

**Features**:
- `Protected`. Options `Method -> Automatic | "GramSchmidt" |
  "ModifiedGramSchmidt" | "Householder" | "Reorthogonalization"` and
  `Tolerance -> Automatic | t` (a residual shorter than `t` times its vector
  becomes a zero vector: `Orthogonalize[{{1., 2.}, {3., 4.}}, Tolerance -> 0.5]`
  is `{{0.447…, 0.894…}, {0., 0.}}`).
- Exact, symbolic and arbitrary-precision input take the reference's classical
  form `w_k = v_k - Sum_j (<v_k, w_j>/<w_j, w_j>) w_j`, `e_k = w_k/Norm[w_k]`
  whatever the Method (as the reference does: Householder on an exact matrix
  gives the Gram–Schmidt answer). An exact residual with (Gaussian) rational
  entries is scaled to a primitive integer vector before normalising, giving
  `(2 + 3 I)/Sqrt[26]` as in the reference. Mathilda's arithmetic does not
  merge `13/(14 Sqrt[13/14])`-type products, so some radicals print in an
  equivalent but different form (`13/Sqrt[182]` for the reference's
  `Sqrt[13/14]`).
- Machine numbers (default inner product): raw-double kernels — classical GS,
  modified GS, MGS with a full second pass (`Automatic`,
  `"Reorthogonalization"`), or Householder (LAPACK `geqrf` + `orgqr`, answer
  signs as LAPACK gives them). Default tolerance `64 max(k, n) $MachineEpsilon`;
  at arbitrary precision about `2^(8 - bits)`.
- Messages: `::orbop` (bad Method), `::ornhip` (Householder with a user inner
  product), `::ornhsl` (Householder on a dependent set), `::ornfa` (an inner
  product that does not return a number or symbol), `::ornt` (ragged first
  argument), `::orwdt` (an element unusable with Dot).
- **Compile**: `Orthogonalize[m]` lowers for `_Real` / `_Complex` rank-2
  operands (the Automatic machine kernel; List, packed, `NDArray` and
  `Compile[]` agree bit for bit). A user inner product or options run in the
  interpreter.

```mathematica
In[1]:= Orthogonalize[{{1, I}, {2, 3}}]
Out[1]= {{1/Sqrt[2], I/Sqrt[2]}, {(2 + 3 I)/Sqrt[26], (3 - 2 I)/Sqrt[26]}}

In[2]:= Orthogonalize[{{1, 2}, {3, 4}, {5, 6}}]
Out[2]= {{1/Sqrt[5], 2/Sqrt[5]}, {2/Sqrt[5], -(1/Sqrt[5])}, {0, 0}}
```

## QRDecomposition
Gives the QR decomposition of a matrix.
- `QRDecomposition[m]` — returns `{q, r}` such that
  `m == ConjugateTranspose[q] . r`, with `q` row-orthonormal
  (row-unitary for complex `m`) and `r` upper triangular.
- `QRDecomposition[m, Pivoting -> True]` — returns `{q, r, p}` where
  `p` is a permutation matrix such that `m . p == ConjugateTranspose[q] . r`.
  Pivoting orders the diagonal of `r` by decreasing magnitude.

**Features**:
- `Protected`.
- Computes the "thin" QR factorisation: when `m` has rank `r`, both
  `q` and `r` have `r` rows. For an `n x p` input, `q` has dimensions
  `r x n` and `r` has dimensions `r x p`.
- Works on every input family:
  - exact integer / rational matrices (output stays exact with
    `Sqrt[...]` in the column norms);
  - complex matrices (rows of `q` are unitary in the Hermitian
    inner product);
  - machine-precision Real matrices (output Real at machine
    precision);
  - arbitrary-precision MPFR matrices (output at the input
    precision via the shared rationalise → exact → numericalise
    pipeline);
  - free-symbolic matrices (closed-form symbolic output).
- Algorithm: dispatched on leaf precision.
  - **MachinePrecision inputs (`min_bits <= 53`)** use a LAPACK
    Householder kernel (`dgeqrf` / `dgeqp3` for real, `zgeqrf` /
    `zgeqp3` for complex, plus `dorgqr` / `zungqr` to form `q`).
    Wired through the four-tier autodetection ladder described in
    `src/linalg/lapack.h` (Apple Accelerate → pkg-config lapacke →
    system lapacke → graceful fall-back).
  - **MPFR inputs (`min_bits > 53`)** use a hand-rolled Householder
    kernel over column-major MPFR arrays (paired re/im planes for
    complex; no MPC dependency, same convention as the eigen
    kernels).  Column pivoting follows Businger-Golub; numerical
    rank uses the cutoff `|R[i,i]| < 2^(-bits/2) * |R[0,0]|`.
    Reconstruction error scales as `2^(-bits)`, matching the
    requested working precision.
  - **Exact / symbolic inputs** stay on the Modified Gram-Schmidt
    kernel, driven through the evaluator so symbolic-real, exact
    rational, and free-variable inputs share one code path.
  - Rank-deficient inputs produce a shorter `q` / `r` without error.
    With `Pivoting -> False` on a rank-deficient input the MPFR
    kernel bails to symbolic (which handles mid-stream rank
    deficiency cleanly); with `Pivoting -> True` the MPFR kernel
    truncates the output at the numerical rank.
  - When BLAS/LAPACK is unavailable at build time (`USE_LAPACK=0`),
    machine-precision inputs transparently route to the symbolic
    kernel; similarly `USE_MPFR=0` routes MPFR inputs to symbolic.
    Correctness is preserved across every combination; only
    performance changes.
- Issues `QRDecomposition::matrix` and returns unevaluated if the
  argument is not a non-empty rank-2 tensor.
- Issues `QRDecomposition::opts` and returns unevaluated for an
  unknown option key or value. `TargetStructure -> "Structured"` is
  reserved for a future release and currently leaves the call
  unevaluated.

```mathematica
In[1]:= QRDecomposition[{{1, 2}, {3, 4}}]
Out[1]= {{{1/Sqrt[10], 3/Sqrt[10]}, {3/Sqrt[10], -(1/Sqrt[10])}},
         {{Sqrt[10], 7 Sqrt[2/5]}, {0, Sqrt[2/5]}}}

In[2]:= {q, r} = QRDecomposition[{{1, 2}, {3, 4}, {5, 6}}];
        Transpose[q] . r
Out[2]= {{1, 2}, {3, 4}, {5, 6}}

In[3]:= {q, r} = QRDecomposition[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}];
        {Length[q], Length[r]}
Out[3]= {2, 2}      (* rank-deficient -- "thin" QR *)

In[4]:= {q, r} = QRDecomposition[{{1.2, 2.3, 3.4},
                                    {2.3, 4.5, 5.6},
                                    {3.2, 7.6, 6.5}}];
        Chop[Transpose[q] . r - {{1.2, 2.3, 3.4},
                                  {2.3, 4.5, 5.6},
                                  {3.2, 7.6, 6.5}}]
Out[4]= {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}

In[5]:= QRDecomposition[{{1, 2}, {3, 4}}, Pivoting -> True]
Out[5]= {q, r, {{0, 1}, {1, 0}}}    (* column 1 (larger norm) pivoted first *)
```

> Implementation is split across:
> - `src/linalg/qrdecomp.c` -- builtin entry, option parsing, and
>   `qr_dispatch` which routes to the precision-matched kernel.
>   Hosts the symbolic / fallback MGS pipeline.
> - `src/linalg/qrdecomp_machine.c` -- LAPACK fast path
>   (`qr_machine_dispatch`).  Loads the matrix into a column-major
>   double buffer (interleaved re/im pairs for complex), calls the
>   wrappers in `lapack.c` (`mat_lapack_dgeqp3`, etc.), then
>   reconstructs `q` / `r` / `p` as Mathilda lists.  Numerical rank
>   uses LAPACK's standard cutoff `max(m, n) * eps * |R[0,0]|`.
> - `src/linalg/qrdecomp_mpfr.c` -- MPFR Householder fast path
>   (`qr_mpfr_dispatch`).  Loads the matrix into column-major MPFR
>   arrays at `min_bits` precision (paired re/im planes for complex),
>   runs Householder reflections in place with Businger-Golub
>   pivoting, then reconstructs `q` / `r` / `p` as MPFR-precision
>   Mathilda lists.  Updates already-stored R rows in-step with
>   column swaps so R's column ordering stays consistent with the
>   pivoted A.  Reconstruction residuals scale as `2^(-bits)`.
> - `src/linalg/lapack.h` / `lapack.c` -- platform-papering Fortran
>   ABI wrappers shared across machine-precision linalg kernels.
>
> The MGS loop allocates a column-major Q buffer and a row-major R
> buffer of size `min(n,p) x max(n,p)`, frees the unused tail when
> the rank turns out to be smaller, and steals the in-use cells
> into the final `List[List[...]]` result.  For complex inputs the
> `q` entries are conjugated at construction time; real inputs (no
> Complex head, no `I` leaf) skip the conjugation to keep the
> printed form free of `Conjugate[Sqrt[...]]` residues that
> Mathilda's simplifier does not reduce.

### QRDecomposition: exact radicals, machine residual
- An integer / rational matrix is factored by Gram-Schmidt over Q with the
  square roots taken only at the end (`flint_mat_qr`), each entry in the
  reference's canonical radical form `k Sqrt[B/A]`:
  `QRDecomposition[{{1, 2}, {3, 4}}]` is
  `{{{1/Sqrt[10], 3/Sqrt[10]}, {3/Sqrt[10], -1/Sqrt[10]}}, {{Sqrt[10], 7 Sqrt[2/5]}, {0, Sqrt[2/5]}}}`.
  A 5 x 5 integer matrix took a minute; it takes about a millisecond.
- A machine QR drops only trailing rows of `r` that are negligible in their
  entirety; cutting at the first small diagonal discarded live rows (a wide
  65 x 75 matrix had `|q^H . r - m| = 0.84`).
- The LAPACK fast path is for machine precision only; an MPFR matrix keeps
  its precision through the Householder kernel.
- The exact zero matrix answers the reference's one row:
  `QRDecomposition[ConstantArray[0, {3, 4}]]` is `{{{1, 0, 0}}, {{0, 0, 0, 0}}}`.
- Known difference: a rank-deficient machine matrix is truncated to its
  numerical rank (q is rank x n); the reference returns the full
  min(m, n)-row factorisation.

## LUDecomposition
Gives the LU decomposition of a matrix.
- `LUDecomposition[m]` — returns `{lu, p, c}` where:
  - `lu` is the combined Doolittle factor matrix (same shape as `m`).
    The strictly-lower triangle is L with an implicit unit diagonal;
    the upper triangle is U.
  - `p` is the 1-indexed row-permutation vector of length
    `Length[m]` (rows) with `m[[p]] == l . u`.  For square `n x n`
    input, `l = LowerTriangularize[lu, -1] + IdentityMatrix[n]` and
    `u = UpperTriangularize[lu]`.  For rectangular `rows x cols`
    input with `steps = Min[rows, cols]`, `l` is `rows x steps` (unit
    lower diagonal on the leading `steps x steps` block) and `u` is
    `steps x cols` (upper).
  - `c` is an L∞-condition-number estimate for approximate numerical
    square inputs, or the exact Integer `0` for exact / symbolic /
    rectangular `m` (the condition number is mathematically undefined
    for non-square matrices).

**Features**:
- `Protected`.
- Any non-empty rectangular `rows x cols` matrix is accepted.  Empty
  matrices and non-matrix arguments emit `LUDecomposition::matsq` and
  the call is left unevaluated.
- Algorithm: Doolittle's elimination with partial row pivoting.
  - **MachinePrecision inputs (`min_bits <= 53`)** use the LAPACK
    fast path: `dgetrf` / `zgetrf` for the factorisation, plus
    `dgecon` / `zgecon` (with `dlange` / `zlange` for ‖A‖∞) for the
    condition estimate.
  - **MPFR inputs (`min_bits > 53`)** use a hand-rolled kernel over
    row-major MPFR arrays (paired re/im for complex). A real matrix (with
    FLINT) is eliminated in Crout's order: each entry of `L` and `U` is one
    fused dot product (`arb_approx_dot`, rounded once), and the independent
    dot products of each step run on worker threads -- the same pivots and
    permutation as the right-looking Doolittle loop, which complex matrices
    keep. A 256 x 256 matrix at 500 digits factors ~4.5x faster than the
    v0.6225 Doolittle loop (2.5 s against 11.3-13.5 s, back to back on one
    loaded 16-core machine; 3.5 s for the old loop unloaded).
    For real matrices the condition number is estimated by the
    Hager-Higham one-norm iteration on `A^{-T}` (LAPACK's `*lacn2`
    strategy; typically 2-5 triangular-solve pairs, each `O(n^2)`).
    For complex matrices the kernel falls back to the explicit
    inverse via `n` back-substitution pairs (`O(n^3)`).
  - **Exact / symbolic inputs** run Doolittle through the Mathilda
    evaluator, so integer / rational / Complex / Sqrt-bearing /
    free-symbolic entries all share one code path.  Pivot selection:
    when every entry of the active column is an exact numeric
    (`Integer` / `BigInt` / `Rational` / `Complex` of those) the pivot
    with the **smallest absolute value** is chosen — matching
    the reference (e.g. `LUDecomposition[{{1/2, 1/3}, {1/5, 1/7}}]`
    picks the `1/5` pivot, keeping intermediate `L` entries small).
    For any column containing a free symbol, `Sqrt`, or other
    non-exact-numeric leaf, the rule falls back to "first non-zero"
    — matching the reference example
    `LUDecomposition[{{a, b}, {c, d}}]` returning `p = {1, 2}`.
  - When BLAS/LAPACK is unavailable at build time (`USE_LAPACK=0`)
    or MPFR is unavailable (`USE_MPFR=0`) the corresponding fast
    path transparently routes to the symbolic kernel, which itself
    understands inexact input via the standard rationalise → exact
    → numericalise round-trip.
- Singular inputs emit `LUDecomposition::sing` and the factorisation
  completes with a zero on U's diagonal at the singular step
  (matching the reference's behaviour).
- Ill-conditioned numerical inputs emit `LUDecomposition::luc`: the
  factorisation completes but the reported `c` exceeds the
  precision-loss threshold (`1 / $MachineEpsilon` for machine input,
  `2^bits` for MPFR input), matching the reference's
  `LUDecomposition::luc` behaviour.  Both warnings are one-shot per
  process; subsequent calls are silent.

```mathematica
In[1]:= LUDecomposition[{{1, 1, 1}, {2, 4, 8}, {3, 9, 27}}]
Out[1]= {{{1, 1, 1}, {2, 2, 6}, {3, 3, 6}}, {1, 2, 3}, 0}

In[2]:= LUDecomposition[{{a, b}, {c, d}}]
Out[2]= {{{a, b}, {c/a, -((b c)/a) + d}}, {1, 2}, 0}

In[3]:= {lu, p, c} = LUDecomposition[{{1.6, 2.7, 3.6},
                                       {1.2, 3.2, 5.2},
                                       {3.3, 3.4, 6.5}}];
        c
Out[3]= <real condition estimate, ~20-30 for the example matrix>
```

> Implementation is split across:
> - `src/linalg/ludecomp.c` -- builtin entry, top-level dispatcher,
>   and the symbolic Doolittle core driven through the evaluator.
> - `src/linalg/ludecomp_machine.c` -- LAPACK fast path; loads to a
>   column-major double buffer, runs `dgetrf` / `zgetrf` then
>   `dgecon` / `zgecon`, builds `{lu, p, c}` as Mathilda lists.
> - `src/linalg/ludecomp_mpfr.c` -- MPFR Doolittle kernel; row-major
>   MPFR arrays (paired re/im for complex), largest-magnitude pivot
>   selection, ‖A‖∞ * ‖A⁻¹‖∞ condition estimate.

## LUBackSubstitution
Solves a linear system from an LU decomposition.
- `LUBackSubstitution[{lu, p, c}, b]` — the `x` with `m . x == b`, where
  `{lu, p, c} = LUDecomposition[m]`; `b` is a vector, or a matrix whose `k`-th
  solution column solves `m . x == b[[All, k]]`.

**Features**:
- `Protected`.
- Two triangular solves on `b[[p]]`: forward through the unit lower triangle
  of `lu`, back through its upper triangle. `O(n^2)` per right-hand side, so
  factoring once and solving many times costs one elimination.
- The answer has the precision of `lu` and `b` together, as arithmetic does:
  - **machine** (a machine number on either side): doubles read straight off
    a packed array / `NDArray` or a list, real or complex; the result is packed
    and keeps `b`'s presentation (an `NDArray` `b` gives an `NDArray`);
  - **arbitrary precision** (MPFR entries, no machine number): MPFR at the
    smaller of the two precisions, real or complex — the same solver an
    arbitrary-precision `LinearSolveFunction` applies;
  - **exact / symbolic**: through the evaluator, exact throughout (no
    simplification beyond evaluation).
- Messages (the call stays unevaluated): `LUBackSubstitution::matsq` when
  `lu` is not a non-empty square matrix (a rectangular decomposition),
  `LUBackSubstitution::lslc` when `b` does not have `Length[lu]` rows,
  `LUBackSubstitution::singlu` when `U` has a zero on its diagonal. A first
  argument that is not a three-element list, or a `p` that is not a
  permutation of `1..n`, leaves the call unevaluated silently.
- Not lowered by `Compile[]`: its first argument is the heterogeneous triple
  `LUDecomposition` returns (exempt for the same reason).
- Lives in `src/linalg/lubacksub.c`.

```
In[1]:= LUBackSubstitution[LUDecomposition[{{1, 2}, {3, 4}}], {5, 6}]
Out[1]= {-4, 9/2}

In[2]:= LUBackSubstitution[LUDecomposition[{{1, 2}, {3, 4}}], {{5, 1}, {6, 0}}]
Out[2]= {{-4, -2}, {9/2, 3/2}}

In[3]:= LUBackSubstitution[LUDecomposition[{{1, 2 I}, {3, 4}}], {1, 1}]
Out[3]= {7/13 + (4 I)/13, -2/13 - (3 I)/13}

In[4]:= h = N[HilbertMatrix[256], 500];
        x = LUBackSubstitution[LUDecomposition[h], N[ConstantArray[1, 256], 500]];
        Precision[x]
Out[4]= 500.
```

## SingularValueDecomposition
Gives the singular value decomposition of a matrix.
- `SingularValueDecomposition[m]` — returns `{u, sigma, v}` such that
  `m == u . sigma . ConjugateTranspose[v]`, with `u` (`n x n`) and `v`
  (`p x p`) orthonormal-columns matrices and `sigma` an `n x p`
  rectangular matrix carrying the singular values on its leading
  diagonal in non-increasing order.
- `SingularValueDecomposition[m, k]` — returns the truncated SVD using
  the `k` largest singular values (or `|k|` smallest if `k < 0`).
  Shape: `u` is `n x k`, `sigma` is `k x k` square, `v` is `p x k`.
- `SingularValueDecomposition[m, UpTo[k]]` — returns the truncated SVD
  using up to `k` of the largest singular values, clamping `k` to
  `MatrixRank[m]`.
- `SingularValueDecomposition[{m, a}]` — generalized SVD; returns
  `{{u, ua}, {sigma, sigma_a}, v}` such that
  `m == u . sigma . ConjugateTranspose[v]` and
  `a == ua . sigma_a . ConjugateTranspose[v]`.  Both matrices must
  share the same number of columns; truncation `k` is not supported in
  this form.
- A machine `SparseArray` with a count, and `Method -> "Arnoldi"`: see
  "Eigenvalues / Eigenvectors / Eigensystem / SingularValueList /
  SingularValueDecomposition: SparseArray and Method -> \"Arnoldi\"" after `Eigensystem`.

**Features**:
- `Protected`.
- Options: `Tolerance -> t` zeroes out singular values with
  `|sigma_i| < t`; `TargetStructure -> "Dense" | "Structured"`
  (`"Structured"` wraps `sigma` as `DiagonalMatrix[{...}]` -- which
  currently evaluates back to a dense matrix because Mathilda has no
  unevaluated structured-matrix head yet).
- Works on every input family:
  - exact integer / rational matrices (output stays exact, with
    `Sqrt[lambda]` forms for the singular values when the
    characteristic polynomial of `m^H . m` doesn't simplify further);
  - complex matrices (`u` and `v` are unitary in the Hermitian inner
    product);
  - machine-precision Real matrices (output Real at machine precision
    via LAPACK `dgesdd` / `zgesdd`);
  - arbitrary-precision MPFR matrices (output at the input precision
    via a one-sided Jacobi SVD).
- Algorithm: dispatched on leaf precision.
  - **MachinePrecision inputs (`min_bits <= 53`)** use LAPACK
    divide-and-conquer SVD (`dgesdd` for real, `zgesdd` for complex)
    with `jobz = 'A'` so the full square `u` and `V^H` are returned.
    `V^H` is conjugate-transposed in place to produce `v`.
  - **MPFR inputs (`min_bits > 53`)** use a one-sided Jacobi SVD: sweep
    over column pairs of `A`, apply the rotation that diagonalises
    `[<a_i,a_i>, <a_i,a_j>; <a_i,a_j>, <a_j,a_j>]` to columns `i` and
    `j` of both `A` and the accumulating `V`.  After convergence (off-
    diagonal norm below `n · 2^(-bits)` of the column norms -- v0.6143; it
    was `2^(-bits/2)`, which left `u` and `v` accurate to only half the
    displayed digits), column norms of `A` are the
    singular values and `A[:, i] / sigma_i` are the left singular
    vectors.  Orthonormal completion of the `n - rank` remaining
    columns of `u` uses Gram-Schmidt at the working precision.  Real
    and complex MPFR inputs are both supported: the complex path uses
    paired (re, im) MPFR buffers throughout and factors each 2x2
    rotation as a phase rotation (that makes the column inner product
    real-positive) followed by a real Jacobi rotation.
  - **Exact / symbolic inputs** compute the eigendecomposition of
    `m^H . m` (or `m . m^H`, whichever is larger so the null-space
    eigenvectors come back automatically), normalise the eigenvectors
    individually, derive the secondary singular vectors via
    `secondary = m . v_i / sigma_i`, and orthonormal-complete the
    smaller side via the existing `qr_symbolic_core` Gram-Schmidt
    pipeline.

    Because the kernel is an eigendecomposition, it inherits
    `Eigenvectors`' handling of irrational algebraic eigenvalues. A
    full-rank exact matrix whose gram has an **irreducible cubic**
    characteristic polynomial — `{{1,2,3},{4,5,6},{7,8,10}}`, for
    instance — used to come back unevaluated: the gram's eigenvalues are
    *casus irreducibilis*, `Eigenvectors` returned zero vectors, and
    normalising divided by a zero norm. Fixed 2026-08-03 by the adjugate
    route documented under `Eigenvectors`; it now reconstructs to
    `1.7e-34` in 0.81 s. **Radical depth, not rank deficiency, was the
    trigger** — that is a separate failure, fixed 2026-08-02.
  - When BLAS/LAPACK is unavailable at build time (`USE_LAPACK=0`),
    machine-precision inputs transparently route to the symbolic
    kernel; similarly `USE_MPFR=0` routes MPFR inputs to symbolic.
    Correctness is preserved across every combination; only
    performance changes.
- Issues `SingularValueDecomposition::matrix` and returns unevaluated
  if the argument is not a non-empty rank-2 tensor.
- Issues `SingularValueDecomposition::sval` and returns unevaluated
  for an out-of-range `k` or `UpTo[k]`.
- Issues `SingularValueDecomposition::opts` and returns unevaluated
  for an unknown option key or value.
- Issues `SingularValueDecomposition::matdims` and returns unevaluated
  when the generalized form `{m, a}` is called with mismatched column
  counts.
- Issues `SingularValueDecomposition::nogsymb` and returns unevaluated
  when the generalized form is called with free symbolic content (no
  closed-form symbolic kernel exists for the generalized GSVD).
- Issues `SingularValueDecomposition::gmpdwn` (one-shot) when the
  generalized form is called on high-precision MPFR input: the result
  is computed at machine precision via LAPACK because Mathilda does
  not yet have a native MPFR Paige/Van Loan kernel.

**Generalized SVD (`{m, a}`) algorithm**:
- Real or complex machine-precision input is forwarded to LAPACK
  `dggsvd3` (real) or `zggsvd3` (complex).  The output of LAPACK,
  `A = U . D1 . (0 R) . Q^H` and `B = V . D2 . (0 R) . Q^H`, maps to
  the reference's form with `u = U`, `ua = V`, `v = Q`, and
  `sigma = D1 . (0 R)` and `sigma_a = D2 . (0 R)` materialised as
  dense rectangular matrices.  Mathilda handles both LAPACK
  layouts (`M >= K+L` and `M < K+L`, the latter splits `R` across
  the destroyed `A` and `B` buffers).
- Exact-numeric symbolic input (Integer / Rational / Complex of
  numeric) is numericalised to 53-bit Reals and re-dispatched
  through the LAPACK path.  Inputs with free symbolic content emit
  `::nogsymb` and leave the call unevaluated.
- High-precision MPFR input is also routed through the LAPACK path
  (with the `::gmpdwn` one-shot precision-downgrade warning) until a
  native MPFR generalized kernel lands.

```mathematica
In[1]:= SingularValueDecomposition[{{1, 2}, {1, 2}}]
Out[1]= {{{1/Sqrt[2], 1/Sqrt[2]}, {1/Sqrt[2], -1/Sqrt[2]}},
         {{Sqrt[10], 0}, {0, 0}},
         {{1/Sqrt[5], -2/Sqrt[5]}, {2/Sqrt[5], 1/Sqrt[5]}}}

In[2]:= sv = SingularValueDecomposition[{{1.2, 3.4}, {5.6, 7.8}, {9.0, 1.2}}];
        sv[[2]]
Out[2]= {{12.4778, 0.0}, {0.0, 5.65202}, {0.0, 0.0}}

In[3]:= SingularValueDecomposition[{{1.1, 2, 5}, {3, -11, 4.2}}, 1]
        (* "thin" SVD with only the top singular value *)
Out[3]= shape: {2 x 1, 1 x 1, 3 x 1}

In[4]:= SingularValueDecomposition[N[{{1, 2}, {3, 4}}, 30]]
        (* MPFR Jacobi -- output at the input precision *)

In[5]:= SingularValueDecomposition[{{1.0, 0}, {1.0, 10^-14}}, Tolerance -> 10^-15]
        (* explicit Tolerance preserves the tiny singular value *)
Out[5]= {..., {{1.41421, 0.}, {0., 7.07107e-15}}, ...}
```

> Implementation is split across:
> - `src/linalg/svdecomp.c` -- builtin entry, positional + option
>   parser (`svd_parse_args`), top-level dispatcher (`svd_dispatch`),
>   symbolic dispatcher + core (`svd_symbolic_core`), and the shared
>   post-processing helper `svd_apply_postprocess` (truncation /
>   tolerance / TargetStructure).
> - `src/linalg/svdecomp_machine.c` -- LAPACK fast path
>   (`svd_machine_dispatch`).  Loads to a column-major double buffer,
>   calls `mat_lapack_dgesdd` / `mat_lapack_zgesdd`, builds `u`,
>   rectangular `sigma`, and `v = ConjugateTranspose[V^H]` as Mathilda
>   lists.
> - `src/linalg/svdecomp_mpfr.c` -- MPFR one-sided Jacobi kernel
>   (`svd_mpfr_dispatch`).  Column-major MPFR arrays at the input
>   precision; orthonormal completion of `u` via Gram-Schmidt; emits
>   MPFR zeros (not Real 0.0) for off-diagonal sigma entries so
>   downstream Dot products preserve the working precision.

### SingularValueDecomposition: wide MPFR matrices
The one-sided Jacobi kernel needs rows >= columns; a wide MPFR matrix is
decomposed through its conjugate transpose. It overran `U` and crashed
(`SingularValueDecomposition[N[{{1, 2, 3}, {4, 5, 6}}, 30]]`, and
`LeastSquares` on it); the factors now agree with the reference digit for
digit, signs included.

## IdentityMatrix
Generates an identity matrix.
- `IdentityMatrix[n]`: Gives the `n x n` identity matrix.
- `IdentityMatrix[{m, n}]`: Gives the `m x n` identity matrix.

**Features**:
- `Protected`.
- Generates exact integer outputs (`1` on main diagonal, `0` elsewhere).
- A number that is not a positive machine integer (`0`, `-1`, `2.5`, `5/2`,
  `2^70`, `1 + I`) emits `IdentityMatrix::dims` ("Dimension specification 0
  should be a positive machine integer or a pair of positive machine
  integers.") and stays unevaluated. A symbol, or a pair with a non-positive
  or non-integer entry (`{0, 3}`, `{2, -1}`), stays unevaluated silently.

```mathematica
In[1]:= IdentityMatrix[3]
Out[1]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}

In[2]:= IdentityMatrix[{2, 3}]
Out[2]= {{1, 0, 0}, {0, 1, 0}}
```

## Projection
The projection of one vector onto another (v0.6189, `src/linalg/projection.c`).
- `Projection[u, v]`: `((u . Conjugate[v]) / (v . Conjugate[v])) v`.
- `Projection[u, v, f]`: `f[u, v] / f[v, v] v` for an inner product `f`.

**Features**:
- `Protected`.
- A zero `v` gives the zero vector (`Projection[{1, 2}, {0, 0}]` is `{0, 0}`).
- Vectors of different lengths, or a vector beside a non-vector, give `Projection::prnv` and the call stays unevaluated. Two symbols give the formal `(y Conjugate[y] . x) / Conjugate[y] . y`, as in the reference.
- Exact and symbolic entries stay exact; MPFR entries keep their precision.
- Machine vectors take one kernel (`la_projection_machine`) on every surface: List, packed, NDArray (`ndla_projection`) and `Compile`, which lowers it at real and complex rank 1 (ND_FN2S row, `R2_VEC`). So `Compile[...][u, v] === Projection[u, v]`.
- SparseArray: a sparse `v` gives a SparseArray (the answer is a multiple of `v`); a sparse `u` beside a dense `v` gives a dense vector.

## TensorProduct
The tensor product of arrays (v0.6189, `src/linalg/tensorproduct.c`).
- `TensorProduct[t1, t2, ...]`: the outer product of all the elements, an array of rank `r1 + r2 + ...`.

**Features**:
- `Flat`, `OneIdentity`, `Protected`. `TensorProduct[t]` is `t`; `TensorProduct[]` is `1`.
- A number factor multiplies (`TensorProduct[2, {1, 2}]` is `{2, 4}`). A symbol stays a formal factor (`TensorProduct[x, {1, 2}]` and `TensorProduct[a, b]` are left alone, as in the reference).
- The product is computed by the buffer-aware, sparse-native `Outer[Times, ...]`. Any SparseArray factor gives a SparseArray.
- `Compile` lowers `TensorProduct[a, b]` for two machine arrays (ND_FN2S row with the new rank rule `R2_OUTER`, rank `ra + rb`; an integer pair stays integer).

## LowerTriangularMatrixQ
- `LowerTriangularMatrixQ[m]`: True if every entry above the main diagonal is zero.
- `LowerTriangularMatrixQ[m, k]`: True if every entry above the `k`-th diagonal is zero.

**Features** (v0.6189; this builtin was missing):
- `Protected`, with `Tolerance -> t` as for `UpperTriangularMatrixQ`, and rectangular matrices.
- Messages: `::nonopt` for a third positional argument, `::nmint` for a non-integer `k`. `UpperTriangularMatrixQ` now also says `::nmint` for a non-integer `k`, as the reference does; it used to say `::nonopt`.
- The two share one implementation (`triangular_matrix_q`, `src/list/matrixq.c`).

### Matrix predicates of a SparseArray (v0.6189)
`SymmetricMatrixQ`, `HermitianMatrixQ`, `SquareMatrixQ`, `DiagonalMatrixQ`, `UpperTriangularMatrixQ`, `LowerTriangularMatrixQ` and `PositiveDefiniteMatrixQ` of a SparseArray are answered in O(nnz) from the pattern and the stored values (`src/sparse/sparse_la_pred.c`). They used to be densified: O(n^2), and refused past the 2^27-element cap.
- The answers are the dense builtins' answers on `Normal[s]`.
- `PositiveDefiniteMatrixQ` of a machine SparseArray above order 128 tests the Hermitian part with a Cholesky attempt (CHOLMOD when linked).
- Options, a diagonal offset, a nonzero background, or a boxed `HermitianMatrixQ` go through the dense gate.
- `Cross` of SparseArray vectors answers on the dense form.

## KroneckerProduct
The Kronecker (direct) product of arrays.
- `KroneckerProduct[a, b]`: for matrices the block matrix whose `(i, j)` block
  is `a[[i, j]] b`, of dimensions `(m p) x (n q)`.
- `KroneckerProduct[a, b, c, ...]`: folds from the left.

**Features**:
- `Protected`.
- **Vectors.** Two vectors give their outer product (a column times a row):
  `KroneckerProduct[{1, 2}, {3, 4}]` is `{{3, 4}, {6, 8}}`. A leading vector
  beside a matrix is a column, a trailing one a row:
  `KroneckerProduct[{{1, 2}, {3, 4}}, {a, b}]` is
  `{{a, b, 2 a, 2 b}, {3 a, 3 b, 4 a, 4 b}}`.
- **Any rank.** For ranks `ra`, `rb` (`R = ra + rb`) the product regroups the
  tensor product's `R` indices: with `h = Floor[R/2]`, result axis `j <= h`
  joins index `j` (major) with index `h + j` (minor), and an odd `R` leaves the
  last index as a final axis. So equal ranks multiply dimension-wise
  (`{2,3,4} (x) {5,6,7}` is `{10, 18, 28}`), and `{2,3,4} (x) {5,6}` is
  `{8, 15, 6}`. A vector beside an array of even rank `>= 4` gives the outer
  product. With three or more arguments every argument but the last must have
  rank 1 or 2 (otherwise the call stays unevaluated).
- **Values.** Each entry is `a_i b_j`. If both operands consist of numbers and
  either holds a machine number the whole product is machine
  (`KroneckerProduct[{1/2, 2}, {3., 1}]` is `{{1.5, 0.5}, {6., 2.}}`); exact
  integers stay exact past 2^63; a symbolic operand gives `Times` products.
- **Buffers.** Packed / `NDArray` operands are read directly and the product is
  written into a fresh buffer (int64 for two integer operands, abandoning to the
  exact path on overflow; float64 or complex otherwise), inheriting the
  operand's presentation. Compiles: `Compile[{{m, _Real, 2}, {p, _Real, 2}},
  KroneckerProduct[m, p]]`.
- A `SparseArray` operand is taken by the native sparse kernel (see
  [data-structures.md](data-structures.md)).
- Messages: `KroneckerProduct::argm` / `::argmu` for fewer than two arguments;
  `KroneckerProduct::rect` ("Nonrectangular array encountered.") for a number
  or a ragged array. A symbolic operand leaves the call unevaluated silently.

```
In[1]:= KroneckerProduct[{{1, 2}, {3, 4}}, {{0, 5}, {6, 7}}]
Out[1]= {{0, 5, 0, 10}, {6, 7, 12, 14}, {0, 15, 0, 20}, {18, 21, 24, 28}}

In[2]:= KroneckerProduct[{1, 2}, {3, 4}, {5, 6}]
Out[2]= {{15, 18, 20, 24}, {30, 36, 40, 48}}
```

## UpperTriangularize, LowerTriangularize
The triangular parts of a matrix.
- `UpperTriangularize[m]`: `m` with every entry below the leading diagonal
  replaced by zero.
- `UpperTriangularize[m, k]`: keeps the entries on and above the `k`-th
  diagonal (column - row >= `k`; `k > 0` above, `k < 0` below the leading one).
- `LowerTriangularize[m]`, `LowerTriangularize[m, k]`: the entries on and
  below the (k-th) diagonal (column - row <= `k`).

**Features**:
- `Protected`. `m` may be rectangular; `k` any machine integer.
- **The zero.** A matrix holding a machine number anywhere gets the machine
  zero `0.` (`UpperTriangularize[{{1.5, I}, {2, 3}}]` is `{{1.5, I}, {0., 3}}`);
  any other matrix, an arbitrary-precision one included, gets the exact `0`.
- A packed / `NDArray` matrix is answered as a buffer of the same dtype
  (row runs copied, the rest zero); compiles (`Compile[{{m, _Real, 2}, {k,
  _Integer}}, UpperTriangularize[m, k]]`). A `SparseArray` is taken by the
  native sparse kernel.
- Messages: `::matrix` ("Argument ... at position 1 is not a non-empty
  rectangular matrix.") for a vector, a ragged matrix or a non-list; `::intm`
  for a `k` that is not a machine integer; `::nonopt` for a third non-option
  argument; `::argt` for no arguments. `{}` and an array of rank 3 or more stay
  unevaluated silently.

```
In[1]:= UpperTriangularize[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, 1]
Out[1]= {{0, 2, 3}, {0, 0, 6}, {0, 0, 0}}

In[2]:= LowerTriangularize[{{1, 2, 3, 4}, {5, 6, 7, 8}}]
Out[2]= {{1, 0, 0, 0}, {5, 6, 0, 0}}
```

## DiagonalMatrix
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
- **The zeros take the diagonal's exactness.** A machine `Real` anywhere on the
  diagonal makes the whole matrix machine-real — the invented zeros and the
  exact entries alike — so `DiagonalMatrix[{1, 2, 3.}]` is
  `{{1., 0., 0.}, {0., 2., 0.}, {0., 0., 3.}}` and `DiagonalMatrix[{1/2, 1.}]`
  is `{{0.5, 0.}, {0., 1.}}`. An exact diagonal keeps exact zeros. Only
  *machine* `Real` is contagious: an MPFR entry keeps them exact
  (``DiagonalMatrix[{1.`30, 2}]``), and a symbolic entry stays symbolic while
  the zeros around it still turn `Real` (`DiagonalMatrix[{a, 1.}]` is
  `{{a, 0.}, {0., 1.}}`).
- A uniformly exact-integer or uniformly machine-real result is a
  [packed list](packed-arrays.md) written directly, without building the
  elements — which is what the exactness rule above buys: a matrix of two heads
  fits no buffer, and the exact zeros used to cost a `Real` diagonal 320×
  against NumPy's `np.diag`. The single-vector form **compiles**
  (`Compile[{{v, _Real, 1}}, DiagonalMatrix[v]]`, rank 1 → rank 2). See
  [`packed-arrays.md`](packed-arrays.md#which-array-heads-compile).

```mathematica
In[1]:= DiagonalMatrix[{a, b, c}]
Out[1]= {{a, 0, 0}, {0, b, 0}, {0, 0, c}}

In[2]:= DiagonalMatrix[{a, b}, 1]
Out[2]= {{0, a, 0}, {0, 0, b}, {0, 0, 0}}

In[3]:= DiagonalMatrix[{1, 2, 3}, 0, {3, 5}]
Out[3]= {{1, 0, 0, 0, 0}, {0, 2, 0, 0, 0}, {0, 0, 3, 0, 0}}

In[4]:= DiagonalMatrix[{1., 2., 3.}]
Out[4]= {{1., 0., 0.}, {0., 2., 0.}, {0., 0., 3.}}
```

## HilbertMatrix
Generates a Hilbert matrix, whose `(i, j)` entry is `1/(i + j - 1)`.
- `HilbertMatrix[n]`: Gives the `n x n` Hilbert matrix.
- `HilbertMatrix[{m, n}]`: Gives the `m x n` Hilbert matrix.

**Features**:
- `Protected`.
- Entries are exact `Rational`s by default. The matrix is symmetric and
  notoriously ill-conditioned, making it a standard test case for numeric
  linear-algebra routines.
- The `WorkingPrecision` option chooses the entry representation:
  - `WorkingPrecision -> Infinity` (default): exact `Rational`s.
  - `WorkingPrecision -> MachinePrecision`: machine-precision `Real`s.
  - `WorkingPrecision -> d`: `d`-digit arbitrary-precision (MPFR) `Real`s.
    A digit count at or below machine precision (or a build without MPFR)
    degrades to machine `Real`s.

| Option | Default | Meaning |
|--------|---------|---------|
| `WorkingPrecision` | `Infinity` | precision at which to create entries |

**Diagnostics** (the call is returned unevaluated):
```
  HilbertMatrix::argx: HilbertMatrix called with 0 arguments; 1 argument is expected.
  HilbertMatrix::dims: Dimension specification <spec> should be a positive machine integer or a pair of positive machine integers.
  HilbertMatrix::nonopt: Options expected (instead of <expr>) beyond position 1 in HilbertMatrix[...]. An option must be a rule or a list of rules.
```

```mathematica
In[1]:= HilbertMatrix[3]
Out[1]= {{1, 1/2, 1/3}, {1/2, 1/3, 1/4}, {1/3, 1/4, 1/5}}

In[2]:= HilbertMatrix[{3, 5}]
Out[2]= {{1, 1/2, 1/3, 1/4, 1/5}, {1/2, 1/3, 1/4, 1/5, 1/6}, {1/3, 1/4, 1/5, 1/6, 1/7}}

In[3]:= HilbertMatrix[3, WorkingPrecision -> MachinePrecision]
Out[3]= {{1., 0.5, 0.333333}, {0.5, 0.333333, 0.25}, {0.333333, 0.25, 0.2}}

In[4]:= Det[HilbertMatrix[3]]
Out[4]= 1/2160

In[5]:= Inverse[HilbertMatrix[3]]
Out[5]= {{9, -36, 30}, {-36, 192, -180}, {30, -180, 180}}
```

## HankelMatrix
Generates a Hankel matrix — a matrix that is constant along its
antidiagonals. The `(i, j)` entry is `c_{i+j-1}` when `i + j - 1 <= m`, and
`r_{i+j-m}` otherwise (with `m` the number of rows). Hankel matrices arise in
approximation theory, functional analysis, numerical analysis, and signal
processing.
- `HankelMatrix[n]`: Gives the `n x n` Hankel matrix whose first row and
  first column are the successive integers `1..n`, with zeros below the main
  antidiagonal.
- `HankelMatrix[{c1, ..., cm}]`: Gives the `m x m` Hankel matrix whose first
  column is the given list, with zeros below the antidiagonal.
- `HankelMatrix[{c1, ..., cm}, {r1, ..., rn}]`: Gives the `m x n` Hankel
  matrix with the first list down the first column and the second list across
  the last row.

**Features**:
- `Protected`.
- Entries are copied verbatim, so symbolic, exact, complex, machine and
  arbitrary-precision entries all flow through unchanged; precision comes from
  the entries themselves (e.g. ``1`20``) or from wrapping the result in `N`.
  The antidiagonal fill is the exact integer `0`.
- For `m = n` the matrix is symmetric, and has real eigenvalues when the
  entries are real.
- The shared corner element `c_m` must equal `r_1`. If they differ, the
  column element `c_m` is used (the formula never reads `r_1`) and a
  `HankelMatrix::crs` warning is emitted; the matrix is still produced.
- An all-integer or all-machine-real argument writes the rank-2 result buffer
  directly (`ndbuild_open`) instead of building `m*n` `Expr` cells — packed out,
  bit-identical to the boxed path. Both the single-vector form
  (`Compile[{{v, _Real, 1}}, HankelMatrix[v]]`, rank 1 → rank 2) and the
  two-vector form (`Compile[{{c, _Real, 1}, {r, _Real, 1}}, HankelMatrix[c, r]]`,
  rank 1 × rank 1 → rank 2) **compile**, delegating to this builtin. See
  [`packed-arrays.md`](packed-arrays.md#which-array-heads-compile).

**Diagnostics**:
```
  HankelMatrix::argb: HankelMatrix called with 0 arguments; between 1 and 3 arguments are expected.
  HankelMatrix::crs: Warning: the column element <c_m> and row element <r_1> at positions <m> and 1 are not the same. Using column element.
```
A first argument that is neither a positive integer nor a list (and any
over-arity call) is returned unevaluated.

```mathematica
In[1]:= HankelMatrix[4]
Out[1]= {{1, 2, 3, 4}, {2, 3, 4, 0}, {3, 4, 0, 0}, {4, 0, 0, 0}}

In[2]:= HankelMatrix[{a, b, c, d}]
Out[2]= {{a, b, c, d}, {b, c, d, 0}, {c, d, 0, 0}, {d, 0, 0, 0}}

In[3]:= HankelMatrix[{x, y, z}, {z, a, b, c, d}]
Out[3]= {{x, y, z, a, b}, {y, z, a, b, c}, {z, a, b, c, d}}

In[4]:= HankelMatrix[{1, 1 + 2 I, 3 + 4 I}]
Out[4]= {{1, 1 + 2 I, 3 + 4 I}, {1 + 2 I, 3 + 4 I, 0}, {3 + 4 I, 0, 0}}

In[5]:= N[HankelMatrix[3]]
Out[5]= {{1., 2., 3.}, {2., 3., 0.}, {3., 0., 0.}}
```

## ToeplitzMatrix
Generates a Toeplitz matrix — a matrix that is constant along its diagonals.
The `(i, j)` entry is `c_{i-j+1}` when `i >= j`, and `r_{j-i+1}` otherwise.
Toeplitz matrices arise in approximation theory, signal processing, statistics
and time series.
- `ToeplitzMatrix[n]`: Gives the `n x n` Toeplitz matrix whose first row and
  first column are the successive integers `1..n` (entry `(i, j)` is
  `|i - j| + 1`, so the matrix is symmetric).
- `ToeplitzMatrix[{c1, ..., cn}]`: Gives the `n x n` symmetric Toeplitz matrix
  whose first column (and first row) is the given list.
- `ToeplitzMatrix[{c1, ..., cm}, {r1, ..., rn}]`: Gives the `m x n` Toeplitz
  matrix with the first list down the first column and the second list across
  the first row.

**Features**:
- `Protected`.
- Entries are copied verbatim, so symbolic, exact, complex, machine and
  arbitrary-precision entries all flow through unchanged; precision comes from
  the entries themselves (e.g. ``1`20``) or from wrapping the result in `N`.
  The single-list form is plain symmetric (no conjugation).
- For `m = n` the matrix is symmetric, and has real eigenvalues when the
  entries are real.
- The shared corner element `c_1` must equal `r_1`. If they differ, the
  column element `c_1` is used (the formula never reads `r_1`, which sits on
  the diagonal as `c_1`) and a `ToeplitzMatrix::crs` warning is emitted; the
  matrix is still produced.
- An all-integer or all-machine-real argument writes the rank-2 result buffer
  directly (`ndbuild_open`); both the single-vector form
  (`Compile[{{v, _Real, 1}}, ToeplitzMatrix[v]]`, rank 1 → rank 2) and the
  two-vector form (`ToeplitzMatrix[c, r]`, rank 1 × rank 1 → rank 2)
  **compile**. See [`packed-arrays.md`](packed-arrays.md#which-array-heads-compile).

**Diagnostics**:
```
  ToeplitzMatrix::argb: ToeplitzMatrix called with 0 arguments; between 1 and 3 arguments are expected.
  ToeplitzMatrix::crs: Warning: the column element <c_1> and row element <r_1> at positions 1 and 1 are not the same. Using column element.
```
A first argument that is neither a positive integer nor a list (and any
over-arity call) is returned unevaluated.

```mathematica
In[1]:= ToeplitzMatrix[4]
Out[1]= {{1, 2, 3, 4}, {2, 1, 2, 3}, {3, 2, 1, 2}, {4, 3, 2, 1}}

In[2]:= ToeplitzMatrix[{a, b, c, d}]
Out[2]= {{a, b, c, d}, {b, a, b, c}, {c, b, a, b}, {d, c, b, a}}

In[3]:= ToeplitzMatrix[{1, 2, 3, 4, 5}, {1, 6, 7}]
Out[3]= {{1, 6, 7}, {2, 1, 6}, {3, 2, 1}, {4, 3, 2}, {5, 4, 3}}

In[4]:= ToeplitzMatrix[{1, 2, 3}, {1, 4, 5, 6, 7}]
Out[4]= {{1, 4, 5, 6, 7}, {2, 1, 4, 5, 6}, {3, 2, 1, 4, 5}}

In[5]:= N[ToeplitzMatrix[3]]
Out[5]= {{1., 2., 3.}, {2., 1., 2.}, {3., 2., 1.}}
```

## VandermondeMatrix
Generates a Vandermonde matrix — a matrix whose rows are the successive powers
of a sequence of nodes. The `(i, j)` entry is `xi^(j-1)`, so the first column
is all ones, the second column is the nodes themselves, the third their
squares, and so on. Vandermonde matrices arise in polynomial interpolation and
in computing moments in the monomial basis.
- `VandermondeMatrix[{x1, ..., xn}]`: Gives the `n x n` Vandermonde matrix on
  the nodes `xi`.
- `VandermondeMatrix[{x1, ..., xn}, k]`: Gives the `n x k` Vandermonde matrix.

**Features**:
- `Protected`.
- The nodes need not be numerical and need not be distinct. Symbolic nodes stay
  as `Power` expressions; numeric powers fold to their value. Precision comes
  from the nodes themselves (e.g. ``2`20``) or from wrapping the result in `N`.
- The first column is the literal integer `1` (`xi^0`), so a zero node reads as
  `1` there rather than `Indeterminate` — matching the interpolation semantics.
- For distinct nodes the matrix is non-confluent and `LinearSolve[V, b]`
  recovers the coefficients `{a0, a1, ...}` of the polynomial
  `p(x) = a0 + a1 x + ...` through the points `{xi, bi}`. The determinant is the
  product of node differences `prod_{i<j} (xj - xi)`.
- An all-integer or all-machine-real node list writes the rank-2 result buffer
  directly. The integer path uses checked powers and falls back to the exact
  bignum path on `int64` overflow, so `VandermondeMatrix[{2, 100}, 20]` is still
  exact. (A *mixed* int+real list keeps its unpacked, two-head result, since the
  `Power` cells of an integer node stay exact beside a real one.) The
  single-vector form **compiles** (`Compile[{{v, _Real, 1}}, VandermondeMatrix[v]]`,
  rank 1 → rank 2). See [`packed-arrays.md`](packed-arrays.md#which-array-heads-compile).

**Diagnostics**:
```
  VandermondeMatrix::argt: VandermondeMatrix called with 0 arguments; 1 or 2 arguments are expected.
```
An empty node list, a non-list first argument, a non-positive integer `k`, the
single-argument structured-array conversion form `VandermondeMatrix[vmat]`
(unsupported — Mathilda has no structured-array representation), and any
over-arity call are returned unevaluated.

```mathematica
In[1]:= VandermondeMatrix[{x1, x2, x3, x4}]
Out[1]= {{1, x1, x1^2, x1^3}, {1, x2, x2^2, x2^3}, {1, x3, x3^2, x3^3}, {1, x4, x4^2, x4^3}}

In[2]:= VandermondeMatrix[{2, 3, 5}]
Out[2]= {{1, 2, 4}, {1, 3, 9}, {1, 5, 25}}

In[3]:= VandermondeMatrix[{a, b, c}, 2]
Out[3]= {{1, a}, {1, b}, {1, c}}

In[4]:= Factor[Det[VandermondeMatrix[{a, b, c}]]]
Out[4]= (-a + b) (-a + c) (-b + c)

In[5]:= LinearSolve[VandermondeMatrix[{1, 2, 3}], {6, 11, 18}]
Out[5]= {3, 2, 1}
```

## LinearSolve
Finds `x` that solves the matrix equation `m . x == b`.
- `LinearSolve[m, b]`
- `LinearSolve[m, b, Method -> "<name>"]`

**Features**:
- `Protected`.
- The matrix `m` may be square or rectangular.
- The argument `b` may be a vector (in which case the result is a
  vector) or a matrix (in which case the result is a matrix whose
  `k`-th column solves `m . x == b[[All, k]]`).
- **FLINT acceleration** (when built with FLINT): a square, nonsingular,
  all-integer/rational system is solved exactly via `fmpq_mat_solve` in
  polynomial time; the unique solution matches the classical division-free
  result. Non-square, singular, or symbolic systems fall through to the
  classical solver (which handles rectangular / underdetermined / inconsistent
  cases). Exposed directly as `` FLINT`LinearSolve `` (see *Structural
  Manipulation*).
- Higher-rank tensor inputs are supported. A rank-N `m` with
  dimensions `{d_1, ..., d_{N-1}, n}` is interpreted as a
  `(d_1 * ... * d_{N-1}) x n` linear system whose leading dimensions
  combine into rows; `b` then has dimensions
  `{d_1, ..., d_{N-1}, e_1, ..., e_p}` and the result has dimensions
  `{n, e_1, ..., e_p}`. When `p == 0` the result is a flat vector of
  length `n`.
- For under-determined systems LinearSolve returns one particular
  solution, with every free (non-pivot) variable set to 0; `Solve`
  returns the general solution.
- Issues `LinearSolve::nosol` and returns unevaluated when no solution
  exists.
- Issues `LinearSolve::matrix` / `::lvec` / `::lvec1` and returns
  unevaluated for shape errors.
- Default method is fraction-free Gauss-Jordan elimination on the
  augmented matrix `[m | b]` (the same Bareiss-like routine used by
  `RowReduce` and `Inverse`), so exact integer, rational, and
  symbolic inputs flow through with no spurious denominator blow-up.
- Lives in `src/linalg/linsolve.c`.
- Accepts an optional `Method -> "<name>"` argument:
  - `Method -> Automatic` or `Method -> "Automatic"` — default (alias for `"DivisionFreeRowReduction"`).
  - `Method -> "DivisionFreeRowReduction"` — Bareiss-like fraction-free Gauss-Jordan on `[m | b]`. Recommended for exact / symbolic inputs.
  - `Method -> "OneStepRowReduction"` — classical Gauss-Jordan with one division per pivot per element on `[m | b]`. Each per-cell update is canonicalised via `Together` so symbolic rationals reduce. Fast on numeric matrices.
  - `Method -> "CofactorExpansion"` — Cramer's rule. Requires square non-singular `m`; emits `LinearSolve::cofnsq` on a non-square `m`, and `LinearSolve::cofsng` on a structurally singular `m`. For matrix `b` the rule is applied column-by-column.
- Unknown method names emit `LinearSolve::method` and the call remains unevaluated.
- A singular machine-real system that has solutions sets each free variable to
  `0.`, not the exact `0`: `LinearSolve[{{1., 2.}, {2., 4.}}, {3., 6.}]` is
  `{3., 0.}`. One with no solution emits `LinearSolve::nosol` once; neither
  recurses to `$RecursionLimit` any more.

```mathematica
In[1]:= LinearSolve[{{r, s}, {t, u}}, {y, z}]
Out[1]= {(u y - s z)/(r u - s t), (r z - t y)/(r u - s t)}

In[2]:= LinearSolve[{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}]
Out[2]= {{-3, -4}, {4, 5}}

In[3]:= LinearSolve[{{1, 5}, {2, 6}, {3, 7}, {4, 8}}, {9, 10, 11, 12}]
Out[3]= {-1, 2}

In[4]:= LinearSolve[{{1, 2, 3}, {4, 5, 6}}, {6, 15}]
Out[4]= {0, 3, 0}

In[5]:= LinearSolve[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, {1, 1, 2}]
LinearSolve::nosol: Linear equation encountered that has no solution.
Out[5]= LinearSolve[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, {1, 1, 2}]

In[6]:= LinearSolve[{{1, 2}, {3, 4}}, {5, 6}, Method -> "CofactorExpansion"]
Out[6]= {-4, 9/2}

In[7]:= LinearSolve[{{1, 5}, {2, 6}, {3, 7}, {4, 8}}, {9, 10, 11, 12}, Method -> "OneStepRowReduction"]
Out[7]= {-1, 2}

(* Rank-3 m and rank-4 b: the leading {2, 3} dimensions of m combine
 * into 6 rows of a 6x6 system; the trailing {4, 5} dimensions of b
 * pass through to the result. *)
In[8]:= a = RandomInteger[{-3, 3}, {2, 3, 6}];
        b = RandomInteger[{-3, 3}, {2, 3, 4, 5}];
        Dimensions[LinearSolve[a, b]]
Out[8]= {6, 4, 5}
```


### LinearSolve: sparse matrices, machine methods, LinearSolveFunction (v0.6180)

- **`LinearSolve[m]`** and `LinearSolve[m, opts]` give a
  `LinearSolveFunction[{rows, cols}, {key, m, {opts}}]` (Length 2, as the
  reference's). `f[b]`, `f[B]` (several right-hand sides),
  `f[b, "T"]` (transpose), `f[b, "C"]` (conjugate transpose), `f[b, "J"]`
  (conjugate), `f[b, "N"]`; `f["Properties"]`, `f["MatrixDimensions"]`,
  `f["ConditionNumber"]` (1-norm estimate). A machine matrix, dense or
  SparseArray, is factored once (UMFPACK / CHOLMOD / LAPACK, see below) and the
  factor kept in a small LRU cache keyed by `key`; because the object carries
  `m`, an evicted factor is recomputed transparently. An **arbitrary-precision**
  square matrix (MPFR entries, no machine number) is likewise factored once,
  by an MPFR LU at its own precision, and each application is two triangular
  solves at that precision (`O(n^2)`; a right-hand side carrying a machine
  number or fewer digits solves afresh). Its singularity test is a pivot at
  or below `64 n 2^(1 - bits) max|m_ij|`, relative to the working precision:
  `LinearSolve[N[HilbertMatrix[256], 500]]` (condition number `1.8*10^389`) is
  regular. (Through v0.6225 an MPFR matrix went through a machine-precision
  `MatrixRank`, was reported `LinearSolve::sing1`, and every application
  re-eliminated.) A square arbitrary-precision `LinearSolve[m, b]` takes the
  same LU and two triangular solves (`n^3/3` fused multiply-adds instead of
  Gauss-Jordan's `n^3`); a matrix singular to working precision falls back to
  the tolerance-aware elimination. Exact and symbolic matrices, the `"Krylov"` method, and
  singular (`LinearSolve::sing1`) or rectangular (`LinearSolve::sqmat1`)
  matrices make an object that solves afresh on every application. Errors:
  `LinearSolveFunction::lslc` (wrong length), `LinearSolve::notrhsm` (bad mode).
- **Machine methods**, for a dense or SparseArray machine matrix:
  `Method -> "Multifrontal"` (sparse LU), `"Cholesky"` (`LinearSolve::npdef`
  when not positive definite), `"Krylov"` or
  `{"Krylov", Method -> "ConjugateGradient" | "BiCGSTAB" | "GMRES",
  "Preconditioner" -> "ILU0" | "ILUT" | "ILUTP" | None | Automatic | f,
  "MaxIterations" -> k, Tolerance -> t, "StartingVector" -> v}`. The Krylov
  default is BiCGSTAB with ILU0 (the reference's). Messages: `::kryme` (unknown
  method), `::krypre` (unknown preconditioner), `::krymit` (iteration limit:
  the last iterate is returned). Exact input with `"Krylov"` says `::kryinfp`
  and is solved exactly.
- **SparseArray coefficient matrix.** A machine SparseArray (real or complex,
  zero background) is solved by a native kernel (`src/sparse/sparse_la_solve.c`)
  without densifying: up to order 128 by a dense LAPACK LU (bit-identical to
  the dense matrix), above that by UMFPACK, or CHOLMOD first for a Hermitian
  matrix with a positive diagonal, when SuiteSparse is linked; otherwise by
  the dense LU within the 2^27-element densify cap and by Krylov beyond it.
  The answer is a packed List (a SparseArray right-hand side gives a dense
  answer too). As in the reference, and unlike the dense matrix, a **singular**
  machine SparseArray gives `LinearSolve::sing` and stays unevaluated (even for
  a consistent b), and a **rectangular** one gives `LinearSolve::sqmat`. Exact,
  MPFR and symbolic SparseArrays behave exactly like the dense matrix
  (particular solutions, `::nosol`). An **exact rational** SparseArray of
  order >= 64 that is sparse (density < 10%), or that exceeds the densify cap,
  is solved by an exact sparse LU over FLINT `fmpq` with Markowitz-style
  pivoting (`src/sparse/sparse_la_exact.c`, v0.6191): fill-free systems cost
  O(nnz) rational operations (the exact 20000 x 20000 tridiagonal solves in
  0.04 s; it was refused at the cap), a singular system gets the reference's
  particular solution (free unknowns 0) or `::nosol`; `Det` and `MatrixRank`
  of such matrices come from the same elimination. `Det` and `Inverse` of a machine
  SparseArray above order 128 use the same factorisation (`Inverse::sing`
  for a singular one); `Inverse` answers a List, as in the reference.

### LinearSolve: singular machine systems, conditioning, arbitrary precision
- A square machine system is solved by LU (`?gesv`); `LinearSolve::luc` when
  the infinity-norm rcond is below `10^6 eps`.
- An exactly singular square or a rectangular machine system answers the
  minimum-norm solution, as the reference does:
  `LinearSolve[{{1., 2.}, {2., 4.}}, {1., 2.}]` is `{0.2, 0.4}`. It is
  inconsistent (`LinearSolve::nosol`, unevaluated) when the residual exceeds
  `1e-12 |b|` (probed: `{1, 2 + 3*^-12}` solves, `{1, 2 + 5*^-12}` does not),
  and for a zero matrix even with a zero right-hand side.
- At arbitrary precision the tolerance-aware Gauss-Jordan on `[m | b]` answers
  with every free variable `0` (`LinearSolve[N[{{1, 2, 3}, {4, 5, 6}}, 30], N[{1, 2}, 30]]`
  is `{-1/3, 2/3, 0}` at 30 digits); a 54 x 54 complex system takes
  milliseconds (it took over 10 s).
- **TimeConstrained** is honoured inside the arbitrary-precision LU of
  `LinearSolve`, `LUDecomposition` and `Det`: the factorisation (worker
  threads, no evaluator step) and the loading of the matrix into MPFR arrays
  poll the deadline and unwind cleanly. A 10000 x 10000 factor at 60 digits
  ran more than 12 minutes past a 60 s cap; the overrun is now the time to
  free what was allocated (seconds at that size).

## LeastSquares
Finds `x` that solves the linear least-squares problem for `m . x == b`,
i.e. an `x` minimising `Norm[m . x - b]`.
- `LeastSquares[m, b]`
- `LeastSquares[m, b, Method -> "<name>"]`
- `LeastSquares[m, b, Tolerance -> t]`

**Features**:
- `Protected`.
- The matrix `m` may be square or rectangular and of any rank. When `m`
  has full column rank the minimiser is unique; when `m` is
  rank-deficient the result is the minimum-norm minimiser
  (`PseudoInverse[m] . b`).
- The right-hand side `b` may be a vector or a matrix. When `b` is a
  matrix, the result has one column per RHS — column `j` is the
  least-squares solution for `b[[All, j]]`, i.e. the `x` minimising
  `Norm[m . x - b, "Frobenius"]` over the multi-RHS system.
- Works on every input family supported by `PseudoInverse`:
  integer / rational, symbolic, machine-precision real / MPFR,
  exact and inexact complex.
- Method and Tolerance options may appear in any order; each may
  appear at most once. Duplicates or unknown option names leave the
  call unevaluated.
- Accepted Method names:
  - `Method -> Automatic` or `Method -> "Automatic"` — alias for `"Direct"` (default).
  - `Method -> "Direct"` — Moore-Penrose solve `PseudoInverse[m] . b`. Works on every input family (dense or sparse, exact or numeric, real or complex). The workhorse method. A **machine-real** coefficient matrix takes the same thin SVD `PseudoInverse` does, applying `V`, `S^+` and `U^T` to `b` in turn rather than forming the pseudo-inverse and multiplying: `LeastSquares[A500, b500]` did not finish in 180 s on the exact pipeline and costs 77 ms here. The iterative methods below are deliberately left on the exact path — they are chosen for their iteration, not their answer.  An **arbitrary-precision** real matrix with `rows >= cols` and the Automatic tolerance is solved by Householder QR, `x = LinearSolve[r, q . b]`, when `r` is square with every `|r_ii|` above `2^(-bits/2) max |r_jj|` (full column rank, condition below `10^(digits/2)`, where the SVD would truncate nothing); otherwise, and for an explicit `Tolerance`, the SVD `PseudoInverse`.  A degree-20 fit to 101 points at 60 digits: 45 ms (SVD: 180 ms), within `10^-54` of the exact solution.  Exact matrices take the exact `PseudoInverse`, whose rational products are FLINT matrix products (see `Dot`).
  - `Method -> "IterativeRefinement"` — residual-correction loop on top of Direct: `x <- PseudoInverse[m] . b`, then repeatedly `r = b - m . x`, `dx = PseudoInverse[m] . r`, `x <- x + dx`, capped at 50 iterations and terminated when `Total[Flatten[dx]^2] <= Tolerance^2`. For exact inputs the first correction is exactly zero (Moore-Penrose identity) so the loop converges in one pass; for inexact inputs the loop drives round-off down to Tolerance.
  - `Method -> "Krylov"` — Conjugate-Gradient-on-Least-Squares (Hestenes-Stiefel CG applied to the normal equations). Iterates `q = A p`, `alpha = |s|^2 / |q|^2`, `x <- x + alpha p`, `r <- r - alpha q`, `s = A^T r`, `beta = |s_new|^2 / |s|^2`, `p <- s_new + beta p` from `x_0 = 0` (so the iterate stays in `range(A^T)` and converges to the minimum-norm LS solution for rank-deficient `m`). Capped at `2 cols(m) + 10` iterations. Matrix RHS are solved column-by-column and recombined via `Transpose`. Symbolic inputs (the convergence test is undecidable for them) fall back to Direct.
  - `Method -> "LSQR"` — Paige-Saunders LSQR (ACM TOMS 1982): Lanczos bidiagonalisation of `m` with a Givens-rotation update of the resulting upper triangular factor. Dispatches by input grammar: free-symbol inputs go to Direct; exact (Integer / Rational) and Complex inputs go to Krylov / CGLS (mathematically equivalent without the square-root growth in exact arithmetic); pure-real inputs with at least one Real entry run the canonical double-precision algorithm. Uses the Paige-Saunders estimate `|phi_bar * alpha_{k+1}|` of `||A^T r||` for the convergence test, scaled against the initial gradient `||A^T b||`. Cap is `2 cols(m) + 10` iterations. Detects rank deficiency by monitoring `alpha_new / max(|A_ij|)` and `beta_new / max(|A_ij|)`, avoiding the catastrophic blowup of dividing by a near-zero `alpha`.
- `Tolerance -> Automatic` (default), or a non-negative integer / real /
  `Rational`. Forwarded verbatim as the Tolerance option of the
  underlying `PseudoInverse` call so any future singular-value
  truncation pass in `PseudoInverse` is picked up automatically.
- When `m . x == b` is consistent, `LeastSquares[m, b]` coincides with
  `LinearSolve[m, b]`.
- Satisfies the identity `LeastSquares[m, b] == PseudoInverse[m] . b`.
- Issues `LeastSquares::matrix` / `::lvec` / `::lvec1` and returns
  unevaluated for shape errors.
- Lives in `src/linalg/lstsq.c`.

```mathematica
In[1]:= LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}]
Out[1]= {19/3, 1/2}

In[2]:= LeastSquares[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, {2, -4, 2}]
Out[2]= {0, 0, 0}

In[3]:= LeastSquares[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}}, {1, 2, 4, 8}]
Out[3]= {157/180, 23/90, -13/36}

In[4]:= LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {{7, 1}, {7, 2}, {8, 3}}]
Out[4]= {{19/3, 0}, {1/2, 1}}

In[5]:= LeastSquares[IdentityMatrix[4], {1, 2, 3, 4}]
Out[5]= {1, 2, 3, 4}

In[6]:= LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 7}] ==
        LinearSolve[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 7}]
Out[6]= True

In[7]:= LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8},
                     Method -> "IterativeRefinement", Tolerance -> 1/100]
Out[7]= {19/3, 1/2}

In[8]:= LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, Method -> "Krylov"]
Out[8]= {19/3, 1/2}

In[9]:= LeastSquares[{{3.2, 2.2, 1.2}, {2.1, 7.1, 8.5}, {9.5, 6.7, 3.7}},
                     {7., 8., 9.}, Method -> "LSQR"]
Out[9]= {73.9499, -174.379, 128.329}

In[10]:= LeastSquares[{{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}},
                      {2., -4., 2.}, Method -> "LSQR"]
Out[10]= {0.0, 0.0, 0.0}
```

## Fit

`Fit[data, {f1, ..., fn}, vars]` finds a linear least-squares fit
`a1 f1 + ... + an fn` to `data` for basis functions `fi` of the variables
`vars` (a single symbol `x` or a list `{x, y, ...}`). `Fit[{m, v}]` returns the
coefficient vector `a` minimising `Norm[m . a - v]` for a design matrix `m` and
response vector `v`.

- Data shapes: a flat list `{v1, ..., vn}` (coordinates assumed to be
  `1, 2, ..., n`), univariate pairs `{{x1, v1}, ...}`, or multivariate rows
  `{{x1, ..., xk, v1}, ...}` whose last entry is the response.
- The design matrix has entries `m_ij = f_i(coords_j)`; the fit minimises the
  sum of squared residuals `Sum[(a . f(coords_j) - v_j)^2]`.
- **WorkingPrecision** `-> Automatic` (default) converts exact input to machine
  reals while preserving the precision of already-approximate (`N[..., p]`)
  input; `-> n` uses `n`-digit MPFR arithmetic; `-> Infinity` keeps exact
  rational arithmetic.
- **FitRegularization** `-> {"Tikhonov", lambda}` (aliases `"L2"`,
  `"RidgeRegression"`) minimises `Norm[m.a-v]^2 + lambda Norm[a]^2` (ridge),
  reduced to ordinary least squares on the augmented system
  `[m; Sqrt[lambda] I]` / `[v; 0]`. `-> {"LASSO", lambda}` (alias `"L1"`)
  minimises `Norm[m.a-v]^2 + lambda Norm[a, 1]`, solved by cyclic coordinate
  descent with soft-thresholding (machine precision; selects basis functions by
  driving small coefficients to exactly zero).
- **NormFunction** `-> Function[Norm[#, 1]]` fits in the 1-norm (least absolute
  deviations) via iteratively reweighted least squares; `Norm[#, 2]` is the
  default 2-norm; any other norm (or a norm combined with regularisation) is
  minimised by `FindMinimum`, warm-started from the L2 solution.
- L2 and ridge solves route through `LeastSquares`, so they inherit exact,
  machine and MPFR arithmetic. Malformed data (ragged rows, a variable count
  that does not match the data coordinates, empty data) issues a `Fit::...`
  message and returns unevaluated.
- **Machine fast path.** When the data is fully numeric and `WorkingPrecision`
  resolves to machine precision, the design matrix is built column-wise — each
  basis function is evaluated once with the variables bound to the packed
  coordinate *vectors*, reusing the packed elementwise kernels — and assembled
  as a packed `float64` array, so the solve reaches LAPACK. This accepts a
  packed array or a visible `NDArray[...]` as `data` (previously an `NDArray`
  errored) and takes a 100k-point quadratic fit from ~3.6 s to ~7 ms. Exact
  (`WorkingPrecision -> Infinity`), MPFR, and symbolic-basis fits fall back to
  the general path unchanged.
- **Conditioning.** The plain-L2 path column-scales the design matrix
  (numpy.polyfit-style equilibration) before the solve and unscales the
  coefficients afterwards, which keeps a high-degree monomial fit accurate where
  the raw Vandermonde is ill-conditioned. A basis that is still nearly
  rank-deficient after scaling issues `Fit::cond`; a basis function that is
  numerically zero over the data gets coefficient exactly `0`.
- Attribute `Protected` (and packed-array aware). Lives in `src/fit.c`.

```mathematica
In[1]:= Fit[{{0,1},{1,0},{3,2},{5,4}}, {1, x}, x]
Out[1]= 0.186441 + 0.694915 x

In[2]:= Fit[{{0,1},{1,0},{3,2},{5,4}}, {1, x, x^2}, x, WorkingPrecision -> Infinity]
Out[2]= 135/199 - 53/199 x + 38/199 x^2

In[3]:= Fit[{N[HilbertMatrix[4]], Range[4]}]
Out[3]= {-64.0, 900.0, -2520.0, 1820.0}

In[4]:= Fit[{{0,0,0},{1,0,1},{0,1,2},{1,1,0},{1/2,1/2,1}}, {1, x, y}, {x, y}]
Out[4]= 0.8 - 0.5 x + 0.5 y

In[5]:= Fit[{{0.,0.},{0.001,1},{0.01,1}}, {1, x, x^2}, x,
            FitRegularization -> {"Tikhonov", 1}]
Out[5]= 0.499985 + 0.00549961 x + 5.0496e-05 x^2

In[6]:= Fit[{{0,1},{1,0},{3,2},{5,4}}, {1, x}, x,
            NormFunction -> Function[Norm[#, 1]]]
Out[6]= -1.0 + 1.0 x
```

## DesignMatrix

`DesignMatrix[data, {f1, ..., fn}, vars]` gives the design matrix with entries
`f_i` evaluated at the data coordinates — the matrix `Fit` forms internally. The
data shapes are identical to `Fit`, and `data` may be a nested `List`, a packed
array, or a visible `NDArray[...]`. By default the entries are kept exact; the
`WorkingPrecision` option (`MachinePrecision` or a digit count) converts them to
approximate numbers. Attribute `Protected` (and packed-array aware). Lives in
`src/fit.c`.

```mathematica
In[1]:= DesignMatrix[{{0,1},{1,0},{3,2},{5,4}}, {1, x}, x]
Out[1]= {{1, 0}, {1, 1}, {1, 3}, {1, 5}}

In[2]:= DesignMatrix[{{0,0,0},{1,0,1},{0,1,2}}, {1, x, y}, {x, y}]
Out[2]= {{1, 0, 0}, {1, 1, 0}, {1, 0, 1}}
```

## CharacteristicPolynomial
Gives the characteristic polynomial of a square matrix.
- `CharacteristicPolynomial[m, x]`: the polynomial `Det[m - x I]` in `x`.
- `CharacteristicPolynomial[{m, a}, x]`: the generalised characteristic
  polynomial `Det[m - x a]`.

**Features**:
- `Protected`.
- Entries may be integer, rational, machine- or arbitrary-precision real,
  complex, or symbolic. The result is an expanded polynomial in `x`.
- An all-rational matrix (integers and rationals, List or packed) goes to
  FLINT's multimodular `fmpq_mat_charpoly`: a 40 x 40 integer matrix takes
  0.002 s (was 1.25 s through Faddeev–LeVerrier over `Expr` arithmetic; sympy
  0.22 s).
- Any other ordinary case reuses the eigen module's Faddeev–LeVerrier–Souriau fast
  path (`O(n^4)` matrix multiplications), so the characteristic polynomial of a
  large numeric matrix is computed in polynomial time — the naïve
  `Det[m - x I]` would face an `O(n!)` Laplace expansion of a symbolic-in-`x`
  matrix. The leading `(-1)^n` sign of `Det[m - x I]` is applied so odd-`n`
  polynomials are monic-negative (`CharacteristicPolynomial[IdentityMatrix[3],
  x]` → `1 - 3 x + 3 x^2 - x^3`).
- The generalised case is `Det[m - x a]` via Laplace expansion. A shared null
  space of `m, a` drops the leading term(s): an infinite generalised eigenvalue
  shows as a degree deficit (`CharacteristicPolynomial[{a, b}, x]` of degree
  `< n`).
- The second argument may be a symbol, a number, or any expression (the value
  is substituted for the polynomial variable).
- Called with other than two arguments, emits `CharacteristicPolynomial::argrx`
  and stays unevaluated; a non-square matrix stays unevaluated.

```mathematica
In[1]:= CharacteristicPolynomial[{{1, 2}, {3, 4}}, x]
Out[1]= -2 - 5 x + x^2

In[2]:= CharacteristicPolynomial[{{a, b}, {c, d}}, x]
Out[2]= -b c + a d - a x - d x + x^2

(* Generalised: Det[m - x a]; a shared null space drops the leading term. *)
In[3]:= CharacteristicPolynomial[{{{1, 1, 1}, {1, 0, 1}, {0, 0, 1}},
          {{0, 1, 1}, {0, 1, 1}, {1, 0, 0}}}, x]
Out[3]= -1 - x + x^2
```

## Eigenvalues
Gives a list of the eigenvalues of a square matrix.
- `Eigenvalues[m]`: eigenvalues of the n×n matrix `m`.
- `Eigenvalues[{m, a}]`: generalised eigenvalues of `m` with respect to `a`
  (the values `lambda` such that `m.v == lambda a.v` for some non-zero `v`).
- `Eigenvalues[m, k]`: first `k` eigenvalues largest in absolute value.
- `Eigenvalues[m, -k]`: `k` eigenvalues smallest in absolute value.
- `Eigenvalues[m, UpTo[k]]`: up to `k` eigenvalues.
- A machine `SparseArray` with a count, and `Method -> "Arnoldi"`: see
  "Eigenvalues / Eigenvectors / Eigensystem / SingularValueList /
  SingularValueDecomposition: SparseArray and Method -> \"Arnoldi\"" after `Eigensystem`.

**Features**:
- `Protected`.
- Implemented via the characteristic polynomial `Det[m - lambda I]`
  (or `Det[m - lambda a]`) followed by `Solve`. The ordinary case uses a
  Faddeev–Leverrier–Souriau fast path in `O(n^4)` matrix multiplications,
  so eigenvalues of large rational / diagonal matrices return instantly.
- Approximate (`Real` / MPFR) input flows through Solve's rationalise →
  solve → numericalize pipeline; tiny imaginary noise introduced by the
  Cardano formula on real cubics is chopped automatically.
- Repeated eigenvalues appear with their algebraic multiplicity.
- Eigenvalues are sorted in order of decreasing absolute value. On an exact
  modulus tie — the case that matters is a complex-conjugate pair `a ± b I` —
  the `+imag` member is listed first, matching the reference
  (`Eigenvalues[{{0,-1},{1,0}}]` → `{I, -I}`). This holds for both the numeric
  Direct kernel and the symbolic (`Root[]` / radical) path, so
  `N[Eigenvalues[m]]` and `Eigenvalues[N[m]]` agree on the sign order
  (e.g. `Eigenvalues[{{0,1,0},{0,0,1},{1,1,0}}]` returns the characteristic
  roots as `{Root[…,1], Root[…,3], Root[…,2]}`, whose `N` is
  `{1.32472, -0.662+0.562 I, -0.662-0.562 I}`). Symbolic eigenvalues that carry
  free variables (no concrete modulus) retain Solve's natural order.
- When `m, a` have a shared null space, `Eigenvalues[{m, a}]` returns
  `Infinity` for each degree drop in the characteristic polynomial.
- A numeric pencil `{m, a}` is solved at the precision of its entries: machine
  input through LAPACK (`dsygv` for a symmetric `m` with symmetric positive
  definite `a`, `dggev` otherwise), arbitrary-precision real input in MPFR —
  the Cholesky reduction `L^-1 m L^-T` (`a = L L^T`) and the symmetric kernel
  when `a` is symmetric positive definite, so the spectrum is real by
  construction, else `a^-1 m` by partial-pivoting elimination and the general
  kernel.  A singular `a` (an infinite eigenvalue) and a complex
  arbitrary-precision pencil take the machine / symbolic path.
- Options: `Cubics -> False`, `Quartics -> False` (defaults), `Method`.  With
  the defaults a general irreducible cubic/quartic characteristic polynomial is
  returned as held `Root[]` objects (matching the reference); `Cubics -> True` /
  `Quartics -> True` force explicit radicals for the general case.  The special
  always-solvable families — binomials (`a x^n + b`), quadratic-in-`x^m`, and
  biquadratic-after-depression quartics — are always returned in radical form
  regardless of these options, so a matrix whose eigenvalues are compact nested
  radicals keeps its closed form.  `Root[]` objects numericalize (including
  after numeric substitution), so `N` and the numeric optimizers work either
  way.
- **Arbitrary-precision nonsymmetric matrices** (the Direct kernel):
  Householder reduction to Hessenberg form with every inner product fused
  (`arb_approx_dot`; `v = A u`, `w = A^T u` as fused dots and each entry one
  two-term fused update `A - u y^T - x u^T` per step, rows threaded when
  large), then Francis double-shift QR. For eigenvalues alone the reflectors
  reach only the active block (rows above it and columns past it are never
  read again). A subdiagonal also deflates once it is below `u ||A||_F`
  (`u = 2^-bits`, the input's own rounding level), so the noise eigenvalues
  of a numerically low-rank matrix need no iteration: the 64 x 64 Nystrom
  matrix of a smooth kernel at 20 digits took 83 sweeps (70 of them on ~50
  eigenvalues of size `u ||A||`), 13 now; 0.29 s -> 0.04 s.
- For approximate-numeric (`Real` / MPFR) matrices `Method` selects the
  numerical kernel.  Accepted values: `Automatic` (default; routes
  among Direct / Arnoldi / Banded based on shape and `k`), `"Direct"`,
  `"Arnoldi"`, `"Banded"`, and `"FEAST"`.  Each accepts a sub-option
  list form (`Method -> {"Name", "Key" -> value, ...}`); see the
  per-method sections in the `?Eigenvalues` REPL docstring and the
  weekly changelog entries under `docs/spec/changelog/`.
- Machine real symmetric eigenvalues (no vectors) with LAPACK: a narrow band
  (half-bandwidth `<= n/64`, `n >= 128`) goes to the band solver `dsbevd`
  (`O(n^2 kd)`); otherwise, at `n >= 800` with Accelerate, the two-stage
  `dsyev_2stage`, whose dense-to-band reduction is BLAS-3 (1.8x on a
  1600 x 1600 Laplacian); else `dsyev`.  With LAPACK present `Automatic`
  no longer prefers the in-file Givens `"Banded"` kernel (it is `O(n^3)` with
  full-row updates); `Method -> "Banded"` runs `dsbevd`.

### `Method -> "FEAST"` (interval slice)

For Hermitian (real symmetric or complex Hermitian) inexact numeric
input, `Method -> "FEAST"` returns only the eigenpairs whose
eigenvalues lie in a user-supplied real interval `[a, b]` — a
spectral-slice query rather than a full decomposition.  Implements
the contour-integral spectral-projector algorithm of Polizzi (2009):
`Ne`-point Gauss-Legendre quadrature on the upper half of the
elliptic contour through `(a, 0)` and `(b, 0)` approximates
`P_[a, b](A) Y` (Schwarz symmetry halves the number of complex
linear solves), then a Rayleigh–Ritz reduction with Cholesky
`B_q = L L^*` extracts the in-interval eigenpairs.

Sub-options:

| Key | Default | Meaning |
|-----|---------|---------|
| `"Interval" -> {a, b}` | required | Real interval to slice.  Auto-swapped if `a > b`. |
| `"ContourPoints" -> Ne` | `8` | Gauss-Legendre order.  Supported: `2`, `4`, `8`, `16`. |
| `"SubspaceSize" -> m0` | `Max[20, n/4]` (capped at `n`) | Working subspace dimension; must be `>=` spectral count in `[a, b]`. |
| `"MaxIterations" -> k` | `20` | Outer iteration cap. |
| `"Tolerance" -> t` | precision-aware | Residual stopping criterion. |

Automatic never routes to FEAST; explicit `Method -> "FEAST"` is the
only way in.  The output is sorted by `|lambda|` descending so the
optional `k_spec` (`Eigenvalues[m, k]`, `Eigenvalues[m, -k]`,
`Eigenvalues[m, UpTo[k]]`) composes naturally with the in-interval
filter.  MPFR inputs run a parallel kernel at the input's combined
precision.

**Fail-soft cascade**: FEAST returns `NULL` (and the call falls
through to Direct) on any of: non-Hermitian input, missing
`"Interval"`, `interval_high <= interval_low` (catches degenerate
`{c, c}` and NaN coercion failures), generalised eigenproblem
(`Eigenvalues[{m, a}]`), Cholesky failure on `B_q` (subspace too
small for the spectral count in the interval), LU singular at any
quadrature node, or non-convergence within `MaxIterations`.  The
first such fall-back per process emits a single
`Eigenvalues::feast: ... falling back to the Direct method.` stderr
line tagged with the reason, so an explicit FEAST call always
returns *some* sensible answer — at worst the full Direct spectrum.

```mathematica
In[6]:= Eigenvalues[N[{{2, -1, 0, 0, 0}, {-1, 2, -1, 0, 0},
                       {0, -1, 2, -1, 0}, {0, 0, -1, 2, -1},
                       {0, 0, 0, -1, 2}}],
                    Method -> {"FEAST", "Interval" -> {2.5, 4}}]
Out[6]= {3.7320508075688776, 3.0}

In[7]:= Eigenvalues[N[{{2, -1, 0, 0, 0}, {-1, 2, -1, 0, 0},
                       {0, -1, 2, -1, 0}, {0, 0, -1, 2, -1},
                       {0, 0, 0, -1, 2}}], 1,
                    Method -> {"FEAST", "Interval" -> {0, 4}}]
Out[7]= {3.7320508075688776}
```

```mathematica
In[1]:= Eigenvalues[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= {1/2 (15 + 3 Sqrt[33]), 1/2 (15 - 3 Sqrt[33]), 0}

In[2]:= Eigenvalues[{{a, b}, {c, d}}]
Out[2]= {1/2 (a + d + Sqrt[(-a - d)^2 - 4 (-b c + a d)]),
         1/2 (a + d - Sqrt[(-a - d)^2 - 4 (-b c + a d)])}

In[3]:= Eigenvalues[{{7/2, 0, 1/2, 0}, {0, 3, 0, 1},
                     {1/2, 0, 7/2, 0}, {0, 1, 0, 3}}]
Out[3]= {4, 4, 3, 2}

In[4]:= Eigenvalues[{{{1, 1, 1}, {1, 0, 1}, {0, 0, 1}},
                     {{0, 1, 1}, {0, 1, 1}, {1, 0, 0}}}]
Out[4]= {Infinity, 1/2 (1 + Sqrt[5]), 1/2 (1 - Sqrt[5])}

In[5]:= Eigenvalues[IdentityMatrix[12]]
Out[5]= {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
```

## Eigenvectors
Gives a list of the eigenvectors of a square matrix.
- `Eigenvectors[m]`: eigenvectors of the n×n matrix `m`.
- `Eigenvectors[{m, a}]`: generalised eigenvectors of `m` with respect to `a`.
- `Eigenvectors[m, k]`: first `k` eigenvectors.
- `Eigenvectors[m, UpTo[k]]`: up to `k` eigenvectors.

**Features**:
- `Protected`.
- For each eigenvalue `lambda_i` (with multiplicity μ), Eigenvectors
  computes the null space of `m - lambda_i I` via `RowReduce` and emits
  up to `μ` basis vectors. When the matrix is defective for that
  eigenvalue, the shortfall is padded in-line with zero vectors so the
  `i`-th eigenvector still corresponds positionally to the `i`-th
  eigenvalue.
- The returned list always has length `n` for an `n×n` matrix.
- **Irrational algebraic eigenvalues** bypass `RowReduce` entirely. Its
  pivot test (`is_zero_poly`) is a polynomial-*identity* test that treats
  each distinct radical as an independent generator, so for a *casus
  irreducibilis* eigenvalue — three real roots expressible in radicals
  only through complex cube roots — it cannot prove a genuinely-zero
  entry zero, and reports `m - lambda I` as full rank. When the null
  space comes back **empty** (never a legal answer: every eigenvalue has
  an eigenvector), the eigenvector is instead read off a column of
  `adj(m - x I)` reduced modulo the minimal polynomial `q` of `lambda`.
  The adjugate is built from cofactor determinants — no division, no
  pivoting — so nothing must be decided zero while it is computed, and
  the surviving "is this column zero?" test runs on a univariate
  polynomial over `Q`, where it is exact and complete. `q` is the
  characteristic polynomial itself when that is irreducible (checked
  with `IrreduciblePolynomialQ`), else `MinimalPolynomial[lambda, x]`.
  A *repeated* irrational eigenvalue drops the rank to `n-2` or below,
  vanishing the whole adjugate; the routine then declines rather than
  guessing, and the zero-vector padding above applies.
- For approximate matrices the result is computed in the rationalised
  domain, then numericalized and normalised to unit `Norm`. Exact /
  symbolic matrices return un-normalised eigenvectors.
- Generalised case: vectors that fall in the shared null space of `m`
  and `a` are returned as zero vectors (matching the reference's
  `Eigenvectors::geinsl1` warning behaviour).
- Options: same as Eigenvalues, including the numerical `Method`
  dispatch (`Automatic`, `"Direct"`, `"Arnoldi"`, `"Banded"`,
  `"FEAST"`).  `Method -> "FEAST"` returns the eigenvectors whose
  eigenvalues lie in the supplied `"Interval"` — orthonormal for
  real symmetric input, unitary for complex Hermitian input — and
  shares the same fail-soft cascade as `Eigenvalues`.

```mathematica
In[1]:= Eigenvectors[{{2, 1, 0}, {0, 2, 0}, {0, 0, 1}}]
Out[1]= {{1, 0, 0}, {0, 0, 0}, {0, 0, 1}}

In[2]:= Eigenvectors[{{1, 0, 1}, {0, 1, 0}, {0, 0, 1}}]
Out[2]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 0}}

In[3]:= Norm /@ Eigenvectors[{{1., 2.}, {2., 1.}}]
Out[3]= {1.0, 1.0}

In[4]:= Norm /@ Eigenvectors[{{1, 2}, {2, 1}}]
Out[4]= {Sqrt[2], Sqrt[2]}

(* Irreducible cubic characteristic polynomial: the eigenvector entries
   are polynomials in the eigenvalue.  Returned {{0,0,0},{0,0,0},{0,0,0}}
   before 2026-08-03. *)
In[5]:= FreeQ[Eigenvectors[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}], {0, 0, 0}]
Out[5]= True
```


## Eigensystem
Gives the eigenvalues and the matching eigenvectors of a square matrix.
- `Eigensystem[m]` — `{Eigenvalues[m], Eigenvectors[m]}`, sorted by decreasing
  absolute value.
- `Eigensystem[m, k]` — the first `k` (`k < 0`: the last `|k|`); also
  `UpTo[k]`. `Eigensystem[m, 0]` is `{{}, {}}`.
- `Eigensystem[{m, a}]` — the generalized eigensystem of `m` with respect to `a`.

**Features**:
- `Protected`. Every argument form and option of `Eigenvalues` /
  `Eigenvectors`: `k`, `-k`, `UpTo[k]`, `{m, a}`, `Method -> Automatic |
  "Direct" | "Arnoldi" | "Banded" | "FEAST"` (with sub-options), `Cubics`,
  `Quartics`.
- **One decomposition.** The numerical kernels (Direct, Arnoldi, Banded, FEAST;
  machine and MPFR, real, complex, symmetric/Hermitian, generalized) run ONCE
  with `MATEIGEN_WANT_BOTH`: they sort the values, apply the one permutation to
  the vectors, and k-select both parts together (`direct_finish` in
  `eigen_direct.c`). A visible `NDArray` / packed matrix takes one LAPACK call
  (`dsyev` / `zheevd` / `dgeev` / `zgeev` with vectors) in `ndla_eigensystem`.
  A complex Hermitian matrix (relative tolerance `1e-12`) takes `zheevd`, on the
  packed path and the boxed path alike, and its eigenvalues are Reals, as in the
  reference (general `zgeev` left roundoff imaginary parts). Every `zheevd`
  answer is checked in `O(n^2)` (`A (V x)` against `V (diag(w) x)`), as the
  `zgees` / `zgeev` answers are. The symbolic
  path (exact or symbolic input, or an inexact matrix no numerical kernel
  accepts) computes the eigenvalues from the characteristic polynomial once and
  returns them together with the eigenvectors it found for them. So the result
  equals `{Eigenvalues[...], Eigenvectors[...]}` in order and pairing — up to the
  last bits of a machine eigenvalue, since a values-only LAPACK call
  (`Eigenvalues`) and a values-and-vectors call (`Eigensystem`) may round
  differently.
- A defective matrix pads with zero vectors in place, as `Eigenvectors` does:
  `Eigensystem[{{1, 1}, {0, 1}}]` is `{{1, 1}, {{1, 0}, {0, 0}}}`.
- Messages (reference text): `::matsq`, `::seqs` (a bad sequence
  specification), `::take` (`|k|` larger than the order — `Eigenvalues` itself
  clips silently).
- `SparseArray` input with a count: the Krylov path below; the whole spectrum
  is densified by the transparency gate.
- **Compile**: not lowered — the result is the ragged pair `{values, vectors}`
  (a vector and a matrix, real or complex by data), which is no single machine
  array (`EXEMPT` in `tools/compile_coverage.py`). The head is on `pack.c`'s
  AWARE list, so a packed argument reaches `ndla_eigensystem` without being
  unpacked.

```mathematica
In[1]:= Eigensystem[{{2, 1}, {1, 2}}]
Out[1]= {{3, 1}, {{1, 1}, {-1, 1}}}

In[2]:= Eigensystem[{{1., 2., 3.}, {4., 5., 6.}, {7., 8., 10.}}, -1]
Out[2]= {{0.198247}, {{0.278296, -0.831847, 0.48019}}}
```

### Eigenvalues / Eigenvectors / Eigensystem / SingularValueList / SingularValueDecomposition: SparseArray and Method -> "Arnoldi"

A few eigenpairs or singular values of a machine `SparseArray` (real or complex
entries, zero background) come from a restarted Krylov method, without
densifying (`src/sparse/sparse_la_eigen.c`, `src/linalg/krylov_schur.c`):

- **Forms**: `Eigenvalues[s, k]`, `[s, -k]`, `[s, UpTo[k]]`, the same for
  `Eigenvectors` and `Eigensystem`, the pencil `{a, b}` (either may be sparse),
  `SingularValueList[s, k]` / `[s, -k]`, `SingularValueDecomposition[s, k]` /
  `[s, -k]`. Sizes far above the densify cap work: the 6 smallest eigenvalues
  of the 90000 x 90000 2-D Laplacian take ~2 s.
- **When the reference densifies instead** (and so does Mathilda, with the same
  answer as the dense matrix): `n <= 14`; more than a fifth of the matrix
  stored (`5 nnz > n^2`); a count over a fifth of the dimension (`5 |k| > n`;
  `min(m, n)` for the singular-value heads); the whole spectrum; the pencil
  under the densify cap unless `Method -> "Arnoldi"`. For a densified request
  on a large, genuinely sparse matrix (`max(m, n) > 100`, `5 nnz <= m n`) the
  reference says `::arh` ("Because finding k out of the n eigenvalues and/or
  eigenvectors is likely to be faster with dense matrix methods ...") for
  machine input and `::arhm` for exact input; so does Mathilda. `|k|` beyond
  the dimension is `::take`, unevaluated.
- **Algorithms**: `k > 0` takes the `k` largest in magnitude of `A` (of
  `b^-1 a` for a pencil); `-k` the smallest, by shift-and-invert about 0
  through a sparse factor (CHOLMOD for a Hermitian matrix with a positive
  diagonal, UMFPACK, or a dense LAPACK LU up to order 128). An exactly
  singular matrix at the implicit shift 0 is shifted by `1e-7 phi ||A||`
  instead -- the eigenvalues are still `sigma + 1/theta` exactly, and the
  reference answers there too. The Krylov layer is ARPACK (`dsaupd` /
  `dnaupd` / `znaupd`, mode 1 on a callback) when linked, else -- and when
  `MATHILDA_NO_ARPACK=1` is set, or a run of ARPACK does not converge -- a
  native Krylov-Schur: full re-orthogonalisation (CGS2), Ritz pairs of the
  projected matrix by LAPACK (`dsyev`/`zheev` for a Hermitian operator,
  `dgeev`/`zgeev` otherwise), thick restart onto the wanted invariant subspace
  (thick-restart Lanczos for a Hermitian operator), real arithmetic for a real
  operator, and every pair confirmed by its explicit residual before the run
  stops. Singular values use the Hermitian eigenproblem of `A^H A` (or `A A^H`,
  the smaller), `(A^H A)^-1` through a factor of `A` for `-k` on a square
  matrix (the normal matrix itself, factored, for a rectangular one); each
  `sigma` is refined as `||A v||`.
- **`Method -> {"Arnoldi", opts}`** (the reference's sub-options):
  `"Shift" -> s` (the eigenvalues nearest `s`, for `k` and `-k` alike;
  `::ssing`, unevaluated, when `a - s b` is exactly singular), `"Criteria" ->
  "Magnitude" | "RealPart" | "ImaginaryPart" | "BothEnds"` (applied to
  `1/(lambda - s)` when a shift or `-k` is in force -- so `-k` with
  `"RealPart"` gives the eigenvalues just above 0; `"BothEnds"` is `k - k/2`
  from the top and `k/2` from the bottom, for a real symmetric matrix only,
  else `::arbsym`; `::eslct` for another name), `"BasisSize"` (an integer
  over 3 and at most `n`, else `::arbss`), `"MaxIterations"` (restarts; when
  reached, `::maxit2` and the current values are returned), `"Tolerance"`,
  `"StartingVector"` (`::svecnm` for a wrong length); `::moptx` for another
  name and `::wrgopt` for a value of the wrong kind. Explicit `"Arnoldi"` always
  takes the Krylov path.
- **Ordering and normalisation**, as the dense answer: decreasing `|lambda|`
  (ties: the more negative real part, then `+imag`, first), a conjugate pair
  split by `k` exactly as the dense answer splits it (the solver computes one
  extra value to see the pair); eigenvalues of a Hermitian matrix (and of a
  real definite pencil) are `Real`. Eigenvectors have unit 2-norm and their
  largest component real and positive -- the Krylov sign / phase is otherwise
  arbitrary (the reference's own sparse and dense signs differ). Singular
  values decrease; a value under the tolerance (`Tolerance`, default
  `max(m, n) eps sigma_max`) is reported as `0.`.
- **Result heads** (the reference's): Lists; `SingularValueDecomposition` of a
  machine SparseArray, with or without `k`, is `{List, SparseArray, List}` -- the
  singular-value matrix a SparseArray with background `0.`.
- **Accuracy**: every returned pair has `||A v - lambda v|| <= ~1e-12 ||A||`.
  Where the reference's own Arnoldi stops short (`::maxit2` on the
  1000 x 1000 Laplacian with `k = 6`, `SingularValueList[s, -3]` ~1e-3 off)
  Mathilda converges.
- **Dense `Method -> "Arnoldi"`** (machine) runs the same driver on the CSR form
  of the matrix, with all the sub-options above, `-k`, and the pencil. It
  replaced a single, un-restarted Arnoldi pass (eigenvector residuals up to
  `3e-2`). Without a count, or for `k >= n - 1`, Direct answers. A dense
  matrix with `Method -> Automatic` and a count is Direct (LAPACK) trimmed --
  `Eigenvalues[m, k] === Take[Eigenvalues[m], k]`, as in the reference (the
  former automatic switch to the single-pass Arnoldi for small `k` was ~1e-3
  off on medium matrices). MPFR input keeps its own Arnoldi.
- A complex general machine matrix (dense, or a SparseArray densified) is now
  solved by LAPACK `zgeev`; the real `2n x 2n` embedding it replaces could pick
  `conj(lambda)` for one member of a clustered pair (`1162.53 + 0.089 I` for
  `1162.53 - 0.089 I`, vector residuals `5e-2`). It remains the fallback
  without LAPACK.

```mathematica
In[1]:= L = KroneckerProduct[#, IdentityMatrix[300, SparseArray]] +
            KroneckerProduct[IdentityMatrix[300, SparseArray], #] &[
          SparseArray[{Band[{1, 1}] -> 2., Band[{1, 2}] -> -1., Band[{2, 1}] -> -1.}, {300, 300}]];
In[2]:= Eigenvalues[L, -3]
Out[2]= {0.000544657, 0.000544657, 0.000217868}

In[3]:= Eigenvalues[SparseArray[{Band[{1, 1}] -> 2., Band[{1, 2}] -> -1., Band[{2, 1}] -> -1.}, {60, 60}], 3,
          Method -> {"Arnoldi", "Shift" -> 1.05}]
Out[3]= {1.15217, 1.06005, 0.970414}
```

## HessenbergDecomposition
Gives the Hessenberg decomposition of a numerical square matrix.
- `HessenbergDecomposition[m]` — `{p, h}` with `p` unitary, `h` upper
  Hessenberg and `m == p . h . ConjugateTranspose[p]`.

**Features**:
- `Protected` (the reference also sets `ReadProtected`; Mathilda never does).
  Options `Pivoting -> False | True` (`True` adds the balancing matrix `d` from
  LAPACK `gebal`, with `m . d == d . p . h . ConjugateTranspose[p]`) and
  `TargetStructure` (dense either way); any other option gives `::optx`.
- Numerical, as in the reference: an exact matrix gives `::schurf` ("received a
  matrix with infinite precision") and a symbolic one `::schurn`, both left
  unevaluated; a non-square argument gives `::matsq`.
- **Machine** (List, packed, `NDArray`): LAPACK `dgehrd` + `dorghr` (real) or
  `zgehrd` + `zunghr` (complex), new `mat_lapack_*` wrappers in
  `src/linalg/lapack.c`. The factors match the reference digit for digit (same
  reflector convention). A 2×2 matrix is already Hessenberg: `p` is the identity.
- **Arbitrary precision**: a Householder reduction at the input precision using
  LAPACK's reflector convention (`zlarfg`: a real `beta` of sign opposite to
  `Re x0`, and the last reflector that only makes a complex subdiagonal entry
  real), so the MPFR factors agree with the machine ones to machine precision.
  Entries the reduction never touches — `p`'s first row and column and `h`
  below the subdiagonal — are exact `0`, as in the reference. Pivoting at
  arbitrary precision is computed at machine precision.
- **Compile**: lowers for `_Real` / `_Complex` rank-2 operands. `{p, h}` is two
  `n × n` factors of one element type, i.e. a rank-3 `2 × n × n` machine array:
  the `ND_FNS` row uses the new `rank_rule` 6 (rank 2 → rank 3) and a delegate
  that stacks the interpreter's own `gehrd`/`orghr` factors, so the compiled
  answer equals the interpreted one element for element. Pivoting is not lowered.

```mathematica
In[1]:= {p, h} = HessenbergDecomposition[{{1., 2., 3.}, {4., 5., 6.}, {7., 8., 10.}}];
        h
Out[1]= {{1., -3.59701, -0.248069}, {-8.06226, 14.8, 2.4}, {0., 0.4, 0.2}}

In[2]:= HessenbergDecomposition[{{1, 2}, {3, 4}}]
HessenbergDecomposition::schurf: HessenbergDecomposition has received a matrix with infinite precision.
Out[2]= HessenbergDecomposition[{{1, 2}, {3, 4}}]
```
