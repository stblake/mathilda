# Linear Algebra

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
- Zero positional arguments emits a Mathematica-compatible

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
- Zero positional arguments emits a Mathematica-compatible

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
  Mathematica-compatible

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
  Mathematica-compatible

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
  Mathematica-compatible

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

## Outer
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
## Det
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

## Cross
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

## Inverse
Gives the inverse of a square matrix.
- `Inverse[m]`
- `Inverse[m, Method -> "<name>"]`

**Features**:
- `Protected`.
- Works on both symbolic and numerical matrices.
- For matrices with approximate real or complex numbers, the inverse is generated to the maximum possible precision given the input.
- Issues `Inverse::sing` warning and returns unevaluated if the matrix is singular.
- Issues `Inverse::matsq` warning and returns unevaluated if the argument is not a non-empty square matrix.
- Satisfies the relation `a . Inverse[a] == Inverse[a] . a == IdentityMatrix[n]`.
- Satisfies the relation `Inverse[a . b] == Inverse[b] . Inverse[a]`.
- Accepts an optional `Method -> "<name>"` argument that selects the inversion algorithm. Shares the same method-name grammar as `RowReduce` and `LinearSolve`.
  - `Method -> Automatic` or `Method -> "Automatic"` (default) — alias for `"DivisionFreeRowReduction"`.
  - `Method -> "DivisionFreeRowReduction"` — Bareiss-like fraction-free Gauss-Jordan elimination on the augmented matrix `[A | I]`. Best choice for exact integer / rational / symbolic input — never produces a denominator larger than necessary.
  - `Method -> "OneStepRowReduction"` — classical Gauss-Jordan on `[A | I]` with one division per pivot per row entry. Each entry is canonicalised via `Together` so symbolic cancellations are still detected. Fast on numeric matrices.
  - `Method -> "CofactorExpansion"` — adjugate / determinant formula `A^-1[i,j] = (-1)^(i+j) Det[M[j,i]] / Det[A]`, with `Det` computed via Laplace cofactor expansion. Time complexity is `O(n! n^2)`; intended for small `n` or closed-form symbolic inverses of small matrices. Emits `Inverse::sing` if `Det[A]` is structurally zero.
- Unknown method names emit `Inverse::method` and the call remains unevaluated.

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
- Inexact (`Real` / MPFR) matrices are rationalised at the input
  precision, computed exactly to preserve rank, and numericalised back.
- The `Tolerance` option accepts `Automatic` (default), a non-negative
  number, or a non-negative `Rational`.
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

## RowReduce
Gives the row-reduced form of the matrix `m`.
- `RowReduce[m]`
- `RowReduce[m, Method -> "<name>"]`

**Features**:
- `Protected`.
- Uses fraction-free division logic to perform exact algorithmic reduction across numerical, rational, and symbolics expressions natively avoiding division errors.
- Lives in `src/linalg/linsolve.c`; the helper primitives (Laplace cofactor determinant, exact polynomial division, tensor flatten / dimensions) are exposed from `src/linalg/util.c` and `src/linalg/det.c`.
- Accepts an optional `Method -> "<name>"` argument:
  - `Method -> Automatic` or `Method -> "Automatic"` (default) — alias for `"DivisionFreeRowReduction"`.
  - `Method -> "DivisionFreeRowReduction"` — Bareiss-like fraction-free Gauss-Jordan. Best for exact integer / rational / symbolic input — never produces a denominator larger than necessary.
  - `Method -> "OneStepRowReduction"` — classical Gauss-Jordan with one division per pivot per element. Each entry is canonicalised via `Together` so symbolic cancellations are still detected. Fast on numeric matrices.
  - `Method -> "CofactorExpansion"` — for a non-singular square matrix, returns the identity (verified via `Det[m] != 0` computed by Laplace cofactor expansion). On singular or rectangular input, falls back to `"DivisionFreeRowReduction"` and emits `RowReduce::cofnsq`.
- Unknown method names emit `RowReduce::method` and the call remains unevaluated.

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

## NullSpace
Gives a basis for the null space of a matrix.
- `NullSpace[m]`
- `NullSpace[m, Method -> "<name>"]`

**Features**:
- `Protected`.
- Returns a list of linearly-independent vectors whose span equals
  `{v : m . v == 0}`. If `m` has full column rank the result is the
  empty list `{}`.
- Works on numerical (Integer / Rational / Real / MPFR / Complex),
  big-integer, and symbolic matrices.
- The matrix `m` may be square or rectangular.
- Basis vectors are returned with the **rightmost free column first**,
  matching Mathematica's ordering.
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
  - *Exact path* (every leaf is exact, no `Tolerance`): routes
    through `RowReduce[m, Method -> "<name>"]` and counts the
    non-zero rows of the RREF, using `is_zero_poly` for structural
    zero. Honors the same `Method` grammar as NullSpace / RowReduce /
    LinearSolve / Inverse (`Automatic` / `"DivisionFreeRowReduction"`
    / `"OneStepRowReduction"` / `"CofactorExpansion"`).
  - *Numerical path* (any inexact leaf, or any explicit `Tolerance`):
    runs partial-pivot Gaussian forward-elimination over a portable
    `double`-complex kernel with tolerance-aware pivot selection: a
    column is skipped when its largest sub-pivot `|entry|` is `<= t`.
    `Method` does not affect this path.
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
  - **MPFR inputs (`min_bits > 53`)** use a hand-rolled Doolittle
    kernel over row-major MPFR arrays (paired re/im for complex).
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
    Mathematica (e.g. `LUDecomposition[{{1/2, 1/3}, {1/5, 1/7}}]`
    picks the `1/5` pivot, keeping intermediate `L` entries small).
    For any column containing a free symbol, `Sqrt`, or other
    non-exact-numeric leaf, the rule falls back to "first non-zero"
    — matching the Mathematica example
    `LUDecomposition[{{a, b}, {c, d}}]` returning `p = {1, 2}`.
  - When BLAS/LAPACK is unavailable at build time (`USE_LAPACK=0`)
    or MPFR is unavailable (`USE_MPFR=0`) the corresponding fast
    path transparently routes to the symbolic kernel, which itself
    understands inexact input via the standard rationalise → exact
    → numericalise round-trip.
- Singular inputs emit `LUDecomposition::sing` and the factorisation
  completes with a zero on U's diagonal at the singular step
  (matching Mathematica's behaviour).
- Ill-conditioned numerical inputs emit `LUDecomposition::luc`: the
  factorisation completes but the reported `c` exceeds the
  precision-loss threshold (`1 / $MachineEpsilon` for machine input,
  `2^bits` for MPFR input), matching Mathematica's
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

## IdentityMatrix
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

```mathematica
In[1]:= DiagonalMatrix[{a, b, c}]
Out[1]= {{a, 0, 0}, {0, b, 0}, {0, 0, c}}

In[2]:= DiagonalMatrix[{a, b}, 1]
Out[2]= {{0, a, 0}, {0, 0, b}, {0, 0, 0}}

In[3]:= DiagonalMatrix[{1, 2, 3}, 0, {3, 5}]
Out[3]= {{1, 0, 0, 0, 0}, {0, 2, 0, 0, 0}, {0, 0, 3, 0, 0}}
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
  - `Method -> "Direct"` — Moore-Penrose solve `PseudoInverse[m] . b`. Works on every input family (dense or sparse, exact or numeric, real or complex). The workhorse method.
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

## Eigenvalues
Gives a list of the eigenvalues of a square matrix.
- `Eigenvalues[m]`: eigenvalues of the n×n matrix `m`.
- `Eigenvalues[{m, a}]`: generalised eigenvalues of `m` with respect to `a`
  (the values `lambda` such that `m.v == lambda a.v` for some non-zero `v`).
- `Eigenvalues[m, k]`: first `k` eigenvalues largest in absolute value.
- `Eigenvalues[m, -k]`: `k` eigenvalues smallest in absolute value.
- `Eigenvalues[m, UpTo[k]]`: up to `k` eigenvalues.

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
- Numeric eigenvalues are sorted in order of decreasing absolute value;
  symbolic eigenvalues retain Solve's natural order.
- When `m, a` have a shared null space, `Eigenvalues[{m, a}]` returns
  `Infinity` for each degree drop in the characteristic polynomial.
- Options: `Cubics -> True` (use radicals for cubics; default true so the
  closed-form pipeline can numericalize), `Quartics -> True`, `Method`.
- For approximate-numeric (`Real` / MPFR) matrices `Method` selects the
  numerical kernel.  Accepted values: `Automatic` (default; routes
  among Direct / Arnoldi / Banded based on shape and `k`), `"Direct"`,
  `"Arnoldi"`, `"Banded"`, and `"FEAST"`.  Each accepts a sub-option
  list form (`Method -> {"Name", "Key" -> value, ...}`); see the
  per-method sections in the `?Eigenvalues` REPL docstring and the
  monthly changelog entries under `docs/spec/changelog/`.

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
- For approximate matrices the result is computed in the rationalised
  domain, then numericalized and normalised to unit `Norm`. Exact /
  symbolic matrices return un-normalised eigenvectors.
- Generalised case: vectors that fall in the shared null space of `m`
  and `a` are returned as zero vectors (matching Mathematica's
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
```

