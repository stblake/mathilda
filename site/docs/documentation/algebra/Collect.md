# Collect

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Collect[expr, x]
    expands expr and gathers terms with the same power of x, returning
    a sum of the form Sum[c_k x^k] with each c_k free of x.
Collect[expr, {x1, x2, ...}]
    collects with respect to each xi in turn (nested grouping).
Collect[expr, x, f]
    applies f to each coefficient before re-assembling the sum, useful
    for f = Simplify or f = Factor.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Collect[a x + b y + c x + d y, x]
Out[1]= (a + c) x + b y + d y

In[2]:= Collect[a x^4 + b x^4 + 2 a^2 x - 3 b x + x - 7, x]
Out[2]= -7 + (1 + 2 a^2 - 3 b) x + (a + b) x^4

In[3]:= Collect[a Sqrt[x] + Sqrt[x] + x^(2/3) - c x + 3x - 2b x^(2/3) + 5, x]
Out[3]= 5 + (1 + a) Sqrt[x] + (1 - 2 b) x^(2/3) + (3 - c) x
```

## Implementation notes

**Algorithm.** `builtin_collect` (in `src/poly/poly.c`) groups an expression by powers of one or more keys, delegating to the recursive worker `collect_internal`. For each key it expands the expression with respect to that key (`expr_expand_patt`), except when the key is itself a `Plus` — there expansion would distribute the subterm and is skipped so e.g. `Collect[a(c+s)+b(c+s), c+s]` stays grouped. Each summand is decomposed into base–power form; terms are bucketed by the exponent at which the key appears (single-base keys group by the rational/symbolic exponent ratio, multi-factor monomial keys by the integer multiplicity from `get_k`). The collected coefficients are summed and, if a third head argument `h` is given, wrapped with `h`. It threads over `List`, equations and inequalities (skipping operator slots in `Inequality`), and recurses across multiple keys.

**Data structures.** Base–power lists (as in `Coefficient`) for term decomposition; results are rebuilt with `internal_times`/`Plus` and re-evaluated through `eval_and_free`.

- `Protected`.
- Automatically threads over lists, equations, inequalities, and logic functions.
- Effectively writes `expr` as a polynomial in `x` or a fractional power of `x`.
- `Collect[expr, var, h]` applies `h` to the expression that forms the coefficient of each term obtained.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Collect[a x^2 + b x^2 + c x + x, x]
Out[1]= (1 + c) x + (a + b) x^2
```

```mathematica
In[1]:= Collect[(1 + x + y)^3, x]
Out[1]= 1 + x^3 + 3 y + 3 y^2 + y^3 + x (3 + 6 y + 3 y^2) + x^2 (3 + 3 y)
```

```mathematica
In[1]:= Collect[Expand[(1 + x)^4 (1 + y)^2], x, Factor]
Out[1]= (1 + y)^2 + 4 x (1 + y)^2 + 6 x^2 (1 + y)^2 + 4 x^3 (1 + y)^2 + x^4 (1 + y)^2
```

```mathematica
In[1]:= Collect[(x + y + z)^2, {x, y}]
Out[1]= x^2 + y^2 + 2 y z + z^2 + x (2 y + 2 z)
```

### Notes

`Collect[expr, x]` expands `expr` and gathers terms by power of `x`, leaving each coefficient free of `x`. A third argument applies a function (e.g. `Factor` or `Simplify`) to each coefficient before reassembly; a list of variables groups by each in turn.
