# E

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`E`**

is the exponential constant e (base of natural logarithms), with numerical value ~= 2.71828.

**`NumericQ[E] is True, and D[E, x] is 0. N[E, prec] evaluates it to any`**

<details>
<summary>Notes</summary>

E is a mathematical constant: it has attributes Constant and Protected, precision.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= Log[E^3]
Out[1]= 3

In[2]:= N[E, 40]
Out[2]= 2.7182818284590452353602874713526624977572

In[3]:= Sum[1/n!, {n, 0, Infinity}]
Out[3]= E

In[4]:= Limit[(1 + 1/n)^n, n -> Infinity]
Out[4]= E
```

## Implementation notes

**Definition.** `E` is Euler's number *e* = 2.71828…, the base of the natural
logarithm and the unique real with `D[E^x, x] == E^x`. It is a constant *symbol*:
row `{ "E", M_E, fill_mpfr_e }` of the `kConstants[]` table in `src/numeric.c` supplies
its numeric value, and the `numericalize_init` loop stamps it `Constant | Protected`
(so `Attributes[E]` is `{Constant, Protected}`; the base `Protected` flag is also set
in the default-attributes table in `src/attr.c`). Interned name `SYM_E`
(`src/sym_names.c`), docstring in `src/info.c`. `E` has no DownValues; the identities
that fold back to it — `Log[E] -> 1`, `Log[E^3] -> 3`, `Exp[1] -> E`, and the series /
limit characterisations `Sum[1/n!, {n,0,Infinity}] -> E`,
`Limit[(1+1/n)^n, n->Infinity] -> E` — are produced by the log-exp, sum and limit
subsystems. `D[E, x] -> 0` follows from `Constant`.

**Representation & numeric value.** `E` is a bare `EXPR_SYMBOL` that prints as `E` and
stays exact under evaluation. `N[E]` returns the machine double `M_E`; `N[E, k]` runs
`fill_mpfr_e`, which evaluates `mpfr_exp(1)` at the requested precision (MPFR ships no
dedicated `const_e`) and returns an `EXPR_MPFR`. Machine value 2.71828;
`N[E, 50] = 2.71828182845904523536028747135266249775724709369996`.

**Usage & limits.** `E` is the fixed point of the exponential closed forms and is
recognised as numeric (`NumericQ[E]` is `True`, via `is_numeric_quantity` in
`src/core.c`). Being transcendental it is never rationalised; its digits are produced
on demand by `N`.

- Attributes `Constant`, `Protected`. `Attributes[E] = {Constant, Protected}`;
  the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[E]` is `True` and
  `D[E, x] = 0`.
- `N[E]` gives the machine value `2.71828`; `N[E, prec]` gives any precision
  (MPFR `mpfr_exp` of 1), e.g.
  `N[E, 50] = 2.71828182845904523536028747135266249775724709369996`.
- Participates in exact numeric work, e.g.
  `Round[E^100] = 26881171418161354484126255515800135873611119`.

**Attributes:** `Constant`, `Protected`.

## References

- Source: [`src/numeric.c`](https://github.com/stblake/mathilda/blob/main/src/numeric.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)
- Tests: [`tests/test_dsolve_m19_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m19_stress.c)

## Notes & additional examples

### Notes

`E` is the exponential constant *e*, the base of the natural logarithm. It is a
protected `Constant`, so `D[E, x]` is `0` and it survives evaluation symbolically
until `N` is applied — `N[E, prec]` returns it to any requested precision via the
MPFR backend. The constant is recognised by the rest of the system, so the
classic limit and series characterisations of *e* both fold back to `E`, and
`Log[E^3]` simplifies to its exponent.
