# PolynomialSqrt

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PolynomialSqrt[p] gives a polynomial s with s^2 == p when p is a perfect square (every non-constant irreducible factor has even multiplicity; the numeric content is carried through Sqrt), and $Failed otherwise. PolynomialSqrt[p, x] treats p as a polynomial in x.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

(1 + x)^2

```mathematica
In[1]:= PolynomialSqrt[x^2 + 2 x + 1]
Out[1]= 1 + x
```

Returned in factored form

```mathematica
In[2]:= PolynomialSqrt[(x^2 - 1)^2]
Out[2]= (-1 + x) (1 + x)
```

Not a perfect square

```mathematica
In[3]:= PolynomialSqrt[x^2 + 1]
Out[3]= $Failed
```

## Implementation notes

**Algorithm.** `builtin_polynomialsqrt` returns a polynomial `s` with `s^2 == p` when
`p` is a perfect square, and `$Failed` otherwise. The input is `PolynomialSqrt[p]` or
`PolynomialSqrt[p, x]` (treating `p` as a polynomial in `x`). After `Expand`, a zero
input maps to `0`. Otherwise:

1. `Factor` `p` into `constant * prod(base_i^e_i)`.
2. Walk the factors. A numeric or (in the 2-arg form) `x`-free factor is **content** and
   is multiplied into a running constant `konst`. A non-constant base needs **even**
   multiplicity `e_i`: take `base_i^{e_i/2}`. A non-constant factor of odd multiplicity
   (including a bare multiplicity-1 irreducible) means `p` is not a perfect square, so
   the whole thing fails.
3. Assemble `s = Sqrt[konst] * prod(base_i^{e_i/2})` and emit it only after an **exact
   certificate**: `Expand[s^2 - p]` must be the zero polynomial. If the check fails
   (or any factor was odd-multiplicity), return `$Failed`.

The 2-arg form carries a symbolic/algebraic leading coefficient through `Sqrt` (e.g.
`Sqrt[-c]` for `-c x^2`) rather than rejecting it as a bare irreducible — needed by
Cherry's Erf-argument construction over an algebraically closed constant field.

**Data structures.** `Expr` trees built with the module's `internal_factor` and the
`expr_expand` / `eval_and_free` helpers; the half-power factors are collected in a small
`Expr**` array before being multiplied together.

**Complexity / limits.** Dominated by `Factor[p]`. Attributes: `Protected`. The square
root is returned in **factored** form (e.g. `(-1 + x)(1 + x)`, not the expanded
`x^2 - 1`), and the numeric content is carried as an exact `Sqrt` when it is not a
perfect square itself. Anything whose factorisation has an odd-multiplicity
non-constant factor is `$Failed`.

**Attributes:** `Protected`.

## References

- Source: [`src/poly/facpoly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/facpoly.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_polynomialsqrt.c`](https://github.com/stblake/mathilda/blob/main/tests/test_polynomialsqrt.c)

## Notes & additional examples

### Notes

`PolynomialSqrt[p]` returns a polynomial `s` with `s^2 == p` when `p` is a perfect
square — every non-constant irreducible factor must have even multiplicity, and the
numeric content is carried through `Sqrt` — and `$Failed` otherwise. `PolynomialSqrt[p,
x]` treats `p` as a polynomial in `x`, so any factor free of `x` counts as constant
content.

The result is given in **factored** form: `PolynomialSqrt[(x^2 - 1)^2]` is
`(-1 + x)(1 + x)`, not the expanded `x^2 - 1`. Every success carries an exact
certificate — `Expand[s^2 - p]` must be zero before `s` is returned — so a near-miss is
reported as `$Failed` rather than an approximate root. `PolynomialSqrt` is `Protected`.
