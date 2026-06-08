# Coefficient

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Coefficient[expr, form]
    gives the coefficient of form^1 in expr.  form is matched
    structurally against the bases of products in the expanded form
    of expr.
Coefficient[expr, form, n]
    gives the coefficient of form^n.  n may be a non-negative integer
    or (for Laurent / Puiseux expressions) a rational.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Coefficient[(x+1)^3, x, 2]
Out[1]= 3

In[2]:= Coefficient[(x+y)^4, x y^3]
Out[2]= 4

In[3]:= Coefficient[x^s x, x^s]
Out[3]= x
```

## Implementation notes

**Algorithm.** `builtin_coefficient` (in `src/poly/poly.c`) extracts the coefficient of `form^n` (default n = 1) from `expr`. It first runs `expr_expand` to a flat sum of monomials, then decomposes both the target `form` and each summand into base–exponent pairs via `decompose_to_bp`. The helper `get_k` computes the integer power at which the target's base(s) divide each term; terms with `k == n` contribute, with the matching base factors stripped out (`get_k` handles multi-factor monomial forms like `x y`). For `n == 0` the whole term is taken. Surviving residual factors are reassembled with `internal_times`/`internal_plus`.

**Data structures.** `BPList` — a list of `{base, exp}` pairs (initialised with `bp_init`, freed with `bp_free`) — is the core representation for both the target form and each term.

- `Protected`, `Listable`.
- `Coefficient[expr, form, 0]` picks out terms that do NOT contain `form`.
- Works whether or not `expr` is explicitly given in expanded form (it automatically expands internally).
- Treats distinct transcendental powers as algebraically unrelated (e.g., `x^s` is treated as a separate base from `x`).

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), Ch. 3 (monomial extraction from polynomial normal forms).
- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Coefficient[x^2 + 3 x + 2, x]
Out[1]= 3
```

```mathematica
In[1]:= Coefficient[x^2 + 3 x + 2, x, 2]
Out[1]= 1
```

```mathematica
In[1]:= Coefficient[a x^2 + b x + c, x, 0]
Out[1]= c
```

```mathematica
In[1]:= Coefficient[3 x^2 y + 2 x y, x, 2]
Out[1]= 3 y
```

### Notes

`Coefficient[expr, form]` returns the coefficient of `form^1` after expanding
`expr`, while the three-argument `Coefficient[expr, form, n]` extracts the
coefficient of `form^n`. Using `n = 0` recovers the term free of `form`, e.g.
the constant `c` in `a x^2 + b x + c`. The extracted coefficient retains any
other variables, so the `x^2` coefficient of `3 x^2 y + 2 x y` is `3 y`. The
`form` is matched structurally against the bases of products, and `n` may be a
non-negative integer (or a rational for Laurent/Puiseux expressions).
