# Norm

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Norm[expr]
    gives the 2-norm of a number, vector, or matrix (Frobenius norm for
    matrices).
Norm[expr, p]
    gives the p-norm: Abs[expr] for scalars; (Sum |xi|^p)^(1/p) for
    vectors with 1 <= p < Infinity; Max[Abs[expr]] for p == Infinity;
    induced operator norm for matrices when p is 1, 2, or Infinity.
```

## Examples

All examples below are verified against the current Mathilda build.

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

## Implementation notes

**Algorithm.** `builtin_norm` uses `get_tensor_dims` to classify the argument and then builds a symbolic norm expression that it evaluates. A rank-0 scalar reduces to `Abs[x]` (the `p` argument is rejected for scalars). For a rank-1 vector (or any tensor when `p` is the string `"Frobenius"`) it flattens to `Expr**` and assembles: `Norm[v, Infinity]` → `Max[Abs[v_i]]`; otherwise (default `p = 2`, or `"Frobenius"` → 2, or a user `p`) → `Power[Sum[Abs[v_i]^p], 1/p]`. Every `Abs`, `Power`, `Plus`, `Max` is created and run through the evaluator (`eval_and_free`), so exact/symbolic/complex entries produce exact symbolic norms.

**Limits.** Genuine matrix p-norms (e.g. the spectral 2-norm via the largest singular value) are not implemented — a rank-≥2 argument with a non-`"Frobenius"` `p` falls through to NULL and the call stays symbolic. Jagged arrays (`rank < 0`) also return NULL.

- `Protected`.
- For scalars, `Norm[z]` is `Abs[z]`.
- For vectors, `Norm[v]` defaults to the 2-norm: `Sqrt[v . Conjugate[v]]`.
- For vectors, `Norm[v, p]` is `Total[Abs[v]^p]^(1/p)`.
- For vectors, `Norm[v, Infinity]` is the $\infty$-norm given by `Max[Abs[v]]`.
- `Norm[m, "Frobenius"]` gives the Frobenius norm of a matrix `m`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/linalg/norm.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/norm.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Norm[{3, 4}]
Out[1]= 5
```

A symbolic 2-norm is kept exact, written with `Abs` so it is valid for complex
components too:

```mathematica
In[1]:= Norm[{a, b, c}]
Out[1]= Sqrt[Abs[a]^2 + Abs[b]^2 + Abs[c]^2]
```

The same vector under different p-norms — the 1-norm, the max (Infinity) norm,
and the exact 3-norm:

```mathematica
In[1]:= Norm[{1, 2, 3, 4}, 1]
Out[1]= 10

In[2]:= Norm[{1, 2, 3, 4}, Infinity]
Out[2]= 4

In[3]:= Norm[{1, 2, 3, 4}, 3]
Out[3]= 10^(2/3)
```

The norm of a complex vector, taken to 40 significant digits:

```mathematica
In[1]:= N[Norm[{1, 2, 3}], 40]
Out[1]= 3.7416573867739413855837487323165493017559
```

### Notes

`Norm[expr]` gives the 2-norm of a scalar or vector; `Norm[expr, p]` gives the
p-norm (`Abs` for scalars; `(Sum |xi|^p)^(1/p)` for vectors with `1 <= p <
Infinity`; `Max[Abs[expr]]` for `p == Infinity`). Results stay exact for exact
input, and symbolic components are written with `Abs` so the formula is valid for
complex entries.
