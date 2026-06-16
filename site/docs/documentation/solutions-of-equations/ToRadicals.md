# ToRadicals

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ToRadicals[expr]
    attempts to express all Root objects in expr in terms of radicals.

ToRadicals can always give expressions in terms of radicals when the
    highest degree of the polynomial that appears in any Root object is
    four.  Binomial Root objects of the form Root[Function[a #^n + b], k]
    are also reduced to radicals for any degree n.  Other Root objects
    of degree five or higher are returned unchanged.
If Root objects in expr contain parameters, ToRadicals[expr] may yield
    a result that is not equal to expr for all values of the parameters.
ToRadicals automatically threads over lists, equations, inequalities,
    and logic functions.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_to_radicals` (`builtin_to_radicals`) converts held `Root[Function[poly], k]` objects into closed-form radical expressions. The top-level walker is a structural recurrence that rebuilds every `EXPR_FUNCTION` node, so a `Root` buried inside `List`/`Equal`/`Less`/`And`/`Or`/... is handled identically.

Per `Root` node: (1) extract the polynomial body, accepting both the `Slot[1]` form `Function[expr]` and the bound-variable form `Function[t, expr]`; (2) substitute the slot/variable with a fresh symbol `x$` so the standard `get_coeff`/`get_degree_poly` univariate machinery applies; (3) dispatch on degree d — `d=1` linear (`-c0/c1`), `d=2` quadratic formula, `d=3` Cardano, `d=4` Ferrari (depressed quartic + resolvent cubic), `d≥5` only the binomial fast-path `a·x^n + b`, otherwise the `Root` is left untouched — each path producing all d radical roots as a fresh `Expr**`; (4) select the k-th root in Mathilda's canonical `Root` ordering by computing `N[Root[poly, k]]` at machine precision (`root_numericalize`) and picking the radical root closest in the complex plane, falling back to the natural per-formula order (index `k-1`) when coefficients are parametric and numeric evaluation is unavailable.

**Data structures.** `Expr*` trees; degree dispatch reuses the polynomial coefficient extractors from `src/poly/poly.c`. Intermediate radical-expression bookkeeping rides on `eval_and_free` and the Plus/Times/Power normalisation. Inputs are borrowed and deep-copied into the output.

**Complexity / limits.** Closed forms exist only up to degree 4 (Abel–Ruffini); degree ≥ 5 is supported solely for binomials. Root selection costs one numeric `Root` evaluation per node.

- `Protected`.
- Closed-form radicals are always returned when the polynomial has degree

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- G. Cardano, *Ars Magna*, 1545 (cubic); L. Ferrari (quartic resolvent, via Cardano).
- Source: [`src/radicals.c`](https://github.com/stblake/mathilda/blob/main/src/radicals.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= ToRadicals[Root[#^2 - 2 &, 1]]
Out[1]= -Sqrt[2]
```

The smaller root of `x^2 + x - 1` is the negative reciprocal of the golden ratio,
recovered exactly in radicals:

```mathematica
In[1]:= ToRadicals[Root[#^2 + # - 1 &, 1]]
Out[1]= 1/2 (-1 - Sqrt[5])
```

Cardano's formula appears automatically for the real root of a depressed cubic:

```mathematica
In[1]:= ToRadicals[Root[1 + #1 + #1^3 &, 1]]
Out[1]= -1/3 ((1/2 (27 + 3 Sqrt[93]))^(1/3) - 3/(1/2 (27 + 3 Sqrt[93]))^(1/3))
```

It threads through the implicit `Root` objects produced by `Solve`, here giving
the three cube roots of two:

```mathematica
In[1]:= ToRadicals[Solve[x^3 - 2 == 0, x]]
Out[1]= {{x -> 2^(1/3)}, {x -> -(-1)^(1/3) 2^(1/3)}, {x -> (-1)^(2/3) 2^(1/3)}}
```

### Notes

`ToRadicals[expr]` rewrites `Root` objects in `expr` using radicals. It always
succeeds when the underlying polynomial has degree at most four (and for
binomial `Root[a #^n + b &, k]` of any degree); degree-five-and-higher Root
objects are returned unchanged. It threads automatically over lists, equations,
and the results of `Solve`.
