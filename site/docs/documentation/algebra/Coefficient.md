# Coefficient

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Coefficient[expr, form]`**

gives the coefficient of form^1 in expr.  form is matched structurally against the bases of products in the expanded form of expr.

**`Coefficient[expr, form, n]`**

gives the coefficient of form^n.  n may be a non-negative integer or (for Laurent / Puiseux expressions) a rational.

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Coefficient[(x+1)^3, x, 2]
Out[1]= 3

In[2]:= Coefficient[(x+y)^4, x y^3]
Out[2]= 4

In[3]:= Coefficient[x^s x, x^s]
Out[3]= x
```

### Applications (6)

```mathematica
In[4]:= Coefficient[x^2 + 3 x + 2, x]
Out[4]= 3

In[5]:= Coefficient[x^2 + 3 x + 2, x, 2]
Out[5]= 1

In[6]:= Coefficient[a x^2 + b x + c, x, 0]
Out[6]= c

In[7]:= Coefficient[3 x^2 y + 2 x y, x, 2]
Out[7]= 3 y

In[8]:= Coefficient[(1 + x)^50, x, 25]
Out[8]= 126410606437752

In[9]:= Coefficient[(a + b)^4, a^2 b^2]
Out[9]= 6
```

## Implementation notes

**Algorithm.** `builtin_coefficient` (in `src/poly/poly.c`) extracts the coefficient of `form^n` (default n = 1) from `expr`. It first runs `expr_expand` to a flat sum of monomials, then decomposes both the target `form` and each summand into base–exponent pairs via `decompose_to_bp`. The helper `get_k` computes the integer power at which the target's base(s) divide each term; terms with `k == n` contribute, with the matching base factors stripped out (`get_k` handles multi-factor monomial forms like `x y`). For `n == 0` the whole term is taken. Surviving residual factors are reassembled with `internal_times`/`internal_plus`.

**Data structures.** `BPList` — a list of `{base, exp}` pairs (initialised with `bp_init`, freed with `bp_free`) — is the core representation for both the target form and each term.

- `Protected`, `Listable`.
- `Coefficient[expr, form, 0]` picks out terms that do NOT contain `form`.
- Works whether or not `expr` is explicitly given in expanded form (it automatically expands internally).
- Treats distinct transcendental powers as algebraically unrelated (e.g., `x^s` is treated as a separate base from `x`).

**Attributes:** `Listable`, `Protected`.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), Ch. 3 (monomial extraction from polynomial normal forms).
- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_expand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expand.c)
- Tests: [`tests/test_poly.c`](https://github.com/stblake/mathilda/blob/main/tests/test_poly.c)
- Tests: [`tests/test_risch_field.c`](https://github.com/stblake/mathilda/blob/main/tests/test_risch_field.c)
- Tests: [`tests/test_risch_hypertangent.c`](https://github.com/stblake/mathilda/blob/main/tests/test_risch_hypertangent.c)

## Notes & additional examples

### Notes

`Coefficient[expr, form]` returns the coefficient of `form^1` after expanding
`expr`, while the three-argument `Coefficient[expr, form, n]` extracts the
coefficient of `form^n`. Using `n = 0` recovers the term free of `form`, e.g.
the constant `c` in `a x^2 + b x + c`. The extracted coefficient retains any
other variables, so the `x^2` coefficient of `3 x^2 y + 2 x y` is `3 y`. The
`form` is matched structurally against the bases of products, and `n` may be a
non-negative integer (or a rational for Laurent/Puiseux expressions). It scales
to high degrees — extracting the central binomial coefficient
`C(50, 25) = 126410606437752` from `(1 + x)^50` directly — and `form` need not be
a single variable: passing the monomial `a^2 b^2` picks out its coefficient `6`
in `(a + b)^4`.
