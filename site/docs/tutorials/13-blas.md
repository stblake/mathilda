# BLAS kernels

This tutorial is a hands-on tour of the ``BLAS` `` context — direct access to
the machine-precision **BLAS** (Basic Linear Algebra Subprograms) kernels that
power Mathilda's numerical linear algebra. If you know the netlib BLAS names
(`ddot`, `dgemv`, `dgemm`, …) you can call them by hand.

Every transcript below was produced by the actual Mathilda binary. Type the
`In[...]` lines yourself (without the prompt) and you will see the same
`Out[...]`.

## When to reach for BLAS

Mathilda's high-level operators (`Dot`, `Norm`, …) are the right tool for
symbolic or exact work. The ``BLAS` `` context is for when you want a *specific*
machine-precision kernel — to benchmark the numerical path, to teach BLAS, or
to transcribe a numerical recipe written in BLAS terms.

The context is available whenever Mathilda is built with BLAS/LAPACK
(`USE_LAPACK`, on by default). If it was built without, the ``BLAS` `` calls
simply stay unevaluated.

## Passing arrays

Array arguments can be ordinary lists or `NDArray[...]` values. Vectors are
flat lists; matrices are lists of rows. You never pass dimensions — they are
read from the shapes.

Real results come back as an `NDArray`; wrap them in `Normal` to see a plain
list.

## Level 1: vectors

The simplest kernels act on one or two vectors:

```mathematica
In[1]:= BLAS`ddot[{1, 2, 3}, {4, 5, 6}]
Out[1]= 32.0

In[2]:= BLAS`dnrm2[{3, 4}]
Out[2]= 5.0

In[3]:= BLAS`idamax[{1, -9, 3}]
Out[3]= 2
```

`ddot` is the dot product, `dnrm2` the Euclidean norm, and `idamax` the
**1-based** index of the largest-magnitude entry (here `-9`, at position 2).

`daxpy` computes `alpha x + y` and `dscal` computes `alpha x`:

```mathematica
In[4]:= BLAS`daxpy[2, {1, 1, 1}, {3, 4, 5}]
Out[4]= NDArray[{5.0, 6.0, 7.0}]

In[5]:= BLAS`dscal[3, {1, 2, 3}]
Out[5]= NDArray[{3.0, 6.0, 9.0}]
```

## Level 2: matrix times vector

`dgemv` forms `alpha A.x + beta y`. Here `alpha = 1`, `beta = 0`, so it is just
`A.x`:

```mathematica
In[6]:= BLAS`dgemv[1, {{1, 2}, {3, 4}}, {1, 1}, 0, {0, 0}]
Out[6]= NDArray[{3.0, 7.0}]
```

The `beta y` term lets you accumulate into an existing vector — pass the vector
you want to add as the last argument and a nonzero `beta`.

## Level 3: matrix times matrix

`dgemm` is the workhorse: `alpha A.B + beta C`, with `A` of shape `m*k`, `B` of
shape `k*n`, and `C` of shape `m*n`.

```mathematica
In[7]:= BLAS`dgemm[1, {{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}, 0, {{0, 0}, {0, 0}}]
Out[7]= NDArray[{{19.0, 22.0}, {43.0, 50.0}}]
```

That matches `{{1,2},{3,4}} . {{5,6},{7,8}}`. The `beta C` term accumulates:
`dgemm[2, A, B, 1, C]` gives `2 A.B + C`.

The symmetric and triangular Level-3 kernels take advantage of matrix
structure. `dsyrk` forms the symmetric product `alpha A.Transpose[A] + beta C`
(for `A` here a `1*2` row, `A.Transpose[A] = 1 + 4 = 5`):

```mathematica
In[8]:= BLAS`dsyrk[1, {{1, 2}}, 0, {{0}}]
Out[8]= NDArray[{{5.0}}]
```

`dtrsm` solves a triangular system `A.X = alpha B` in place. With the
upper-triangular `A = {{2,1},{0,3}}`:

```mathematica
In[9]:= Normal[BLAS`dtrsm[1, {{2, 1}, {0, 3}}, {{2, 0}, {0, 3}}]]
Out[9]= {{1.0, -0.5}, {0.0, 1.0}}
```

## Complex kernels

The `z*` routines take complex data. Complex results come back as lists of
`Complex[...]` (an `NDArray` holds reals only). `zdotc` conjugates its first
argument, `zdotu` does not:

```mathematica
In[10]:= BLAS`zdotc[{Complex[1, 1]}, {Complex[1, 1]}]
Out[10]= 2.0

In[11]:= BLAS`zdotu[{I}, {I}]
Out[11]= -1.0

In[12]:= BLAS`dznrm2[{Complex[3, 4]}]
Out[12]= 5.0
```

`zgemm` is the complex matrix product:

```mathematica
In[13]:= BLAS`zgemm[1, {{Complex[1, 0], Complex[0, 1]}}, {{I}, {1}}, 0, {{0}}]
Out[13]= {{0.0 + 2.0*I}}
```

## When a call stays put

If an entry is symbolic, or the shapes don't line up, the call is left
unevaluated and flows through the evaluator untouched — so you can pattern-match
or inspect it:

```mathematica
In[14]:= BLAS`dgemm[1, {{a, b}}, {{c}, {d}}, 0, {{0}}]
Out[14]= BLAS`dgemm[1, {{a, b}}, {{c}, {d}}, 0, {{0}}]
```

## Where to go next

The [LAPACK tutorial](14-lapack.md) covers the higher-level drivers — solving
systems, factorizations, SVD, and eigenproblems. The full reference is in
[`builtins/blas.md`](../../../docs/spec/builtins/blas.md).
