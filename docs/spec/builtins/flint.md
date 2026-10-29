# FLINT context

Mathilda transparently accelerates a number of core algebraic and numeric
operations with [FLINT](https://flintlib.org/) (Fast Library for Number Theory)
and its Arb/`acb` ball-arithmetic layer. Those kernels normally run *behind* the
ordinary `System`` builtins — `Factor`, `Cancel`, `Det`, `Zeta`, and friends —
choosing the FLINT path automatically when the input is in scope and falling
back to the native implementation otherwise.

The ``FLINT` `` context exposes those same kernels **directly**, with no
dispatch heuristics, coercion, or symbolic post-processing in the way. Each
routine either returns the exact/rigorous FLINT result or leaves the call
unevaluated when its argument is out of FLINT's scope (a non-polynomial, a
matrix with a symbolic entry, a symbolic numeric argument, a pole, …). They are
useful for benchmarking the accelerated path, for forcing the exact kernel, and
for understanding exactly what FLINT contributes to the system.

All ``FLINT` `` routines carry the `Protected` attribute. They are only present
when Mathilda is built with FLINT support (`USE_FLINT`; the numeric routines
additionally require `USE_MPFR`); otherwise the symbols are absent and the calls
stay unevaluated.

Polynomial algebra works over the rationals **Q** (`fmpq_mpoly`), the linear
algebra over rational matrices (`fmpq_mat`), and the numerics over the complex
plane with rigorous `acb` ball arithmetic (`acb_dirichlet`).

---

## Polynomial algebra

Exact, multivariate polynomial operations over Q, computed directly by FLINT's
`fmpq_mpoly` layer.

## FLINT`PolynomialGCD

`` FLINT`PolynomialGCD[a, b] `` gives the monic greatest common divisor of the
polynomials `a` and `b` over the rationals, computed directly via FLINT
(`fmpq_mpoly_gcd`). It is multivariate and returns unevaluated if either
argument is not a polynomial over Q.

```mathematica
In[1]:= FLINT`PolynomialGCD[x^2 - 1, x^2 - x]
Out[1]= -1 + x

In[2]:= FLINT`PolynomialGCD[x^2 - y^2, x^2 + 2 x y + y^2]
Out[2]= x + y
```

## FLINT`Resultant

`` FLINT`Resultant[a, b, x] `` gives the resultant of the polynomials `a` and
`b` eliminating the variable `x`, over the rationals, computed directly via
FLINT (`fmpq_mpoly_resultant`). Other variables are treated as coefficients.
Returns unevaluated if out of scope.

```mathematica
In[1]:= FLINT`Resultant[x^2 - 1, x - 2, x]
Out[1]= 3
```

## FLINT`Factor

`` FLINT`Factor[p] `` gives the irreducible factorisation of the polynomial `p`
over the rationals, computed directly via FLINT (`fmpq_mpoly_factor`), as
`Times[const, factor^exp, …]`. It is multivariate and returns unevaluated if `p`
is not a polynomial over Q.

```mathematica
In[1]:= FLINT`Factor[x^4 - 1]
Out[1]= (-1 + x) (1 + x) (1 + x^2)
```

## FLINT`FactorSquareFree

`` FLINT`FactorSquareFree[p] `` gives the squarefree factorisation of the
polynomial `p` over the rationals, computed directly via FLINT
(`fmpq_mpoly_factor_squarefree`). Returns unevaluated if out of scope.

```mathematica
In[1]:= FLINT`FactorSquareFree[(x - 1)^2 (x + 1)]
Out[1]= (1 + x) (-1 + x)^2
```

---

## Exact linear algebra

Exact rational-matrix operations, computed directly by FLINT's `fmpq_mat`
layer. Every entry must be an integer or a rational; a matrix with any
non-rational entry leaves the call unevaluated.

## FLINT`Det

`` FLINT`Det[m] `` gives the determinant of the square matrix `m` when every
entry is an integer or rational, computed exactly and directly via FLINT
(`fmpq_mat_det`). Returns unevaluated for a matrix with any non-rational entry.

```mathematica
In[1]:= FLINT`Det[{{1, 2}, {3, 4}}]
Out[1]= -2
```

## FLINT`Inverse

`` FLINT`Inverse[m] `` gives the inverse of the square matrix `m` when every
entry is an integer or rational, computed exactly via FLINT (`fmpq_mat_inv`).
Returns unevaluated for a singular or non-rational matrix.

```mathematica
In[1]:= FLINT`Inverse[{{1, 2}, {3, 4}}]
Out[1]= {{-2, 1}, {3/2, -1/2}}
```

## FLINT`LinearSolve

`` FLINT`LinearSolve[m, b] `` solves the square system `m.x == b` exactly via
FLINT (`fmpq_mat_solve`) when `m` is a nonsingular rational matrix and `b` a
rational vector or matrix. Returns unevaluated for a non-square, singular, or
non-rational system.

```mathematica
In[1]:= FLINT`LinearSolve[{{1, 2}, {3, 4}}, {5, 6}]
Out[1]= {-4, 9/2}
```

## FLINT`RowReduce

`` FLINT`RowReduce[m] `` gives the reduced row echelon form of the matrix `m`
when every entry is an integer or rational, computed exactly via FLINT
(`fmpq_mat_rref`). Returns unevaluated for a matrix with any non-rational entry.

```mathematica
In[1]:= FLINT`RowReduce[{{1, 2, 3}, {4, 5, 6}}]
Out[1]= {{1, 0, -1}, {0, 1, 2}}
```

## FLINT`MatrixRank

`` FLINT`MatrixRank[m] `` gives the rank of the matrix `m` when every entry is
an integer or rational, computed exactly via FLINT (`fmpq_mat_rref`). Returns
unevaluated for a matrix with any non-rational entry.

```mathematica
In[1]:= FLINT`MatrixRank[{{1, 2}, {2, 4}}]
Out[1]= 1
```

---

## Polynomial matrices over K[x]

Exact operations on matrices of univariate polynomials over a field K, on FLINT's
generic rings (`gr_poly`, `gr_mat`) with one implementation for three choices of K:

| K | option | coefficients |
|---|---|---|
| Q | none | rationals |
| GF(p), p a word-size prime | `Modulus -> p` | integers in `[0, p)` (a rational input `a/b` is `a b^-1 mod p`) |
| Q(θ) = Q[t]/(μ), μ irreducible over Q | `Extension -> {t, μ}` | polynomials in `t` of degree < deg μ, reduced on input |

Results are expanded polynomials in canonical form. A call stays unevaluated when
it is out of scope: an entry that is not a polynomial in `x` (and `t`), a ragged or
empty matrix, a rank-deficient module, a non-prime or oversized modulus, a reducible
μ, or an unknown option. These routines are the kernels of the divisor-class
arithmetic of `Integrate[f, x, Method -> "Trager"]` (`src/internal/trager/`). The
package's own routines stay as the reference and fallback, selected by
`tragerOption["NativeKernels"]`.

## FLINT`PolynomialHermiteForm

`` FLINT`PolynomialHermiteForm[m, x, opts] `` gives the Hermite normal form of the
K[x]-module spanned by the rows of `m`. The form is upper triangular with monic
pivots, every entry above a pivot has lower degree than the pivot, and zero rows are
dropped. The module must have full column rank. The form is unique, so it identifies
the module.

| option | effect |
|---|---|
| `"DeterminantMultiple" -> D` | `D` is a multiple of the determinant of the module, so the module contains `D K[x]^n`; the rows `D e_k` are adjoined and entries are reduced modulo `D` during the elimination, which bounds their degrees. A constant `D` is ignored. |
| `"Products" -> {gens, table}` | the module spanned by the products `r g` (`r` a row of `m`, `g` in `gens`) in the K[x]-algebra with basis `w_0..w_(n-1)` and `w_i w_j = f w_k`, `table[[i+1, j+1]] = {k, f}`: an ideal product formed inside the kernel. |
| `"Denominator" -> den` | the fractional module `m/den` in lowest terms: returns `{H, den'}` with every factor of `den` *over the prime field* that divides all entries removed, and `den'` monic. Over Q(θ) a factor of `den` that splits over Q(θ) is not removed piecemeal. |

```mathematica
In[1]:= FLINT`PolynomialHermiteForm[{{x^2 + 1, x}, {x, 1 + x}}, x]
Out[1]= {{1, -x^2}, {0, 1 + x + x^3}}

In[2]:= FLINT`PolynomialHermiteForm[{{x^2 + t, t x}, {x, 1 + t}}, x, Extension -> {t, t^2 - 2}]
Out[2]= {{1, -1/2 t x}, {0, 2 + t + x^2}}

In[3]:= FLINT`PolynomialHermiteForm[{{x (x - 1), 0}, {0, x (x + 1)}}, x, "Denominator" -> x^2]
Out[3]= {{{-1 + x, 0}, {0, 1 + x}}, x}
```

## FLINT`PolynomialDetAdjugate

`` FLINT`PolynomialDetAdjugate[m, x, opts] `` gives `{det, adj}` for a square
polynomial matrix over K[x], with `adj . m == det IdentityMatrix`. It takes the field
options and `"Products"` as above; with `"Products" -> {{h}, table}` applied to an
identity matrix, it gives the determinant (the norm) and adjugate of the
multiplication matrix of `h`.

```mathematica
In[1]:= FLINT`PolynomialDetAdjugate[{{x^2 + 1, x}, {x, 1 + x}}, x]
Out[1]= {1 + x + x^3, {{1 + x, -x}, {-x, 1 + x^2}}}
```

## FLINT`PolynomialWeakPopovForm

`` FLINT`PolynomialWeakPopovForm[m, x, shifts, opts] `` reduces the rows of `m`
without changing the module they span, until their leading columns are distinct.
The leading column of a row is the one that maximises
`k Exponent[e_j, x] + shifts[[j]]`, with `k` set by `"DegreeScale" -> k`
(default 1); on a tie the first such column wins. While two rows share a leading
column, the first such pair is taken and the row with the greater key loses its
leading term to `x^s` times the other. The result therefore matches step for step
the package's `reducedRows`, which uses it for the Riemann–Roch reduction of divisor
classes.

```mathematica
In[1]:= FLINT`PolynomialWeakPopovForm[{{x, 1}, {1, x}}, x, {0, 5}, "DegreeScale" -> 2]
Out[1]= {{x, 1}, {1 - x^2, 0}}
```

## FLINT`PolynomialQuotientRemainder

`` FLINT`PolynomialQuotientRemainder[a, b, x, opts] `` gives `{q, r}` with
`a == q b + r` and `Exponent[r, x] < Exponent[b, x]` over K[x].

```mathematica
In[1]:= FLINT`PolynomialQuotientRemainder[t^2 x, 1, x, Extension -> {t, t^2 - 2}]
Out[1]= {2 x, 0}
```

## FLINT`PolynomialExtendedGCD

`` FLINT`PolynomialExtendedGCD[a, b, x, opts] `` gives `{g, {s, t}}` with `g` the
monic gcd over K[x] and `s a + t b == g`.

```mathematica
In[1]:= FLINT`PolynomialExtendedGCD[x^2 - 2, x - t, x, Extension -> {t, t^2 - 2}]
Out[1]= {-t + x, {0, 1}}
```

## FLINT`PolynomialFieldGCD

`` FLINT`PolynomialFieldGCD[a, b, x, opts] `` gives the monic gcd over K[x], without
cofactors. Over Q(θ) it uses the certified modular algorithm of the field-gcd engine.
θ is first scaled to an algebraic integer when μ is not integral. The Euclidean gcd
over a number field is fraction-based, and its remainders' coefficients grow with
the field degree. With this kernel, the splitting field of degree 36 built for a
Trager residue polynomial took 14 s where the Euclidean route took 74 s.

```mathematica
In[1]:= FLINT`PolynomialFieldGCD[x^2 - 2, x^2 - t x, x, Extension -> {t, t^2 - 2}]
Out[1]= -t + x
```

---

## Rigorous numerics

Numeric evaluation of the zeta family, computed to the precision of the argument
(machine precision for exact input) with FLINT's rigorous `acb` ball arithmetic
over the full complex plane and standard branch cuts.

## FLINT`Zeta

`` FLINT`Zeta[s] `` gives the numeric value of the Riemann zeta function at the
numeric argument `s` (real or complex), computed to the precision of `s`
(machine precision for exact `s`) via FLINT's rigorous `acb` arithmetic
(`acb_dirichlet_zeta`). Unevaluated for symbolic `s` or at the pole `s = 1`.

```mathematica
In[1]:= FLINT`Zeta[2]
Out[1]= 1.644934066848226

In[2]:= FLINT`Zeta[4]
Out[2]= 1.082323233711138
```

## FLINT`HurwitzZeta

`` FLINT`HurwitzZeta[s, a] `` gives the numeric value of the Hurwitz zeta
function via FLINT (`acb_dirichlet_hurwitz`), to the precision of the arguments.
Unevaluated for symbolic arguments or at a pole.

```mathematica
In[1]:= FLINT`HurwitzZeta[2, 1/2]
Out[1]= 4.934802200544679
```

## FLINT`PolyGamma

`` FLINT`PolyGamma[n, z] `` gives the numeric value of the `n`-th derivative of
the digamma function (`n = 0` is digamma) via FLINT (`acb_polygamma`), to the
precision of the arguments. Unevaluated for symbolic arguments or at a pole.

```mathematica
In[1]:= FLINT`PolyGamma[1, 2]
Out[1]= 0.6449340668482264
```

## FLINT`StieltjesGamma

`` FLINT`StieltjesGamma[n] `` and `` FLINT`StieltjesGamma[n, a] `` give the
numeric value of the `n`-th Stieltjes constant (generalized, at `a`) for a
non-negative integer `n`, via FLINT (`acb_dirichlet_stieltjes`). Unevaluated for
negative or non-integer `n`.

```mathematica
In[1]:= FLINT`StieltjesGamma[0]
Out[1]= 0.5772156649015329
```
