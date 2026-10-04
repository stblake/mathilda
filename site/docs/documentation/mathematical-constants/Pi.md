# Pi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Pi`**

is pi, with numerical value ~= 3.14159.

**`NumericQ[Pi] is True, and D[Pi, x] is 0. N[Pi, prec] evaluates it to any`**

<details>
<summary>Notes</summary>

Pi is a mathematical constant: it has attributes Constant and Protected, precision.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (8)

```mathematica
In[1]:= N[Pi]
Out[1]= 3.14159

In[2]:= N[Pi, 40]
Out[2]= 3.1415926535897932384626433832795028841971

In[3]:= Sin[Pi/6]
Out[3]= 1/2

In[4]:= Cos[Pi/3] + Sin[Pi/4]^2
Out[4]= 1

In[5]:= ArcTan[1]
Out[5]= 1/4 Pi

In[6]:= Cos[Pi/5]
Out[6]= 1/4 (1 + Sqrt[5])

In[7]:= N[Pi^2/6, 40]
Out[7]= 1.644934066848226436472415166646025189219

In[8]:= D[Pi, x]
Out[8]= 0
```

## Implementation notes

**Definition.** `Pi` is the circle constant π = 3.14159…, the ratio of a circle's
circumference to its diameter. It is a first-class constant *symbol*, not a function:
the `kConstants[]` table in `src/numeric.c` registers it with a machine value and an
MPFR filler, and the init loop at the bottom of `numericalize_init` stamps every row
of that table `Constant | Protected`, so `Attributes[Pi]` is `{Constant, Protected}`.
The interned name is `SYM_Pi` (`src/sym_names.c`); the docstring is in `src/info.c`.
`Pi` carries no DownValues of its own. The exact closed forms that *produce* it —
`Sin[Pi] -> 0`, `Cos[Pi] -> -1`, `ArcTan[1] -> Pi/4`, `Exp[I Pi] -> -1` — live in the
trigonometric and log-exp heads, and `D[Pi, x] -> 0` follows from the `Constant`
attribute.

**Representation & numeric value.** Internally `Pi` is a bare `EXPR_SYMBOL`; it
survives evaluation unevaluated and `FullForm[Pi]` is `Pi`. `N[Pi]` returns the
machine double `M_PI` (via `leaf_from_double`). `N[Pi, k]` goes through
`numericalize_symbol`, which calls the filler `fill_mpfr_pi` → MPFR's `mpfr_const_pi`
at the requested bit precision and returns an `EXPR_MPFR`. Machine value 3.14159;
`N[Pi, 50] = 3.1415926535897932384626433832795028841971693993751`.

**Usage & limits.** `Pi` drives the special-value reductions of the elementary and
inverse-trigonometric functions and appears throughout the Zeta / Gamma / series
results; `NumericQ[Pi]` is `True` (the symbol is whitelisted in `is_numeric_quantity`,
`src/core.c`). It is never reduced to a rational — being irrational, its only
"evaluation" is that its digits are materialised on demand by `N`.

- Attributes `Constant`, `Protected`. `Attributes[Pi] = {Constant, Protected}`;
  the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[Pi]` is `True` and
  `D[Pi, x] = 0`.
- `N[Pi]` gives the machine value `3.14159`; `N[Pi, prec]` gives any precision
  (MPFR `mpfr_const_pi`), e.g.
  `N[Pi, 50] = 3.1415926535897932384626433832795028841971693993751`.
- Participates in exact numeric work, e.g.
  `Round[Pi^100] = 51878483143196131920862615246303013562686760680406`.

**Attributes:** `Constant`, `Protected`.

## References

- Source: [`src/numeric.c`](https://github.com/stblake/mathilda/blob/main/src/numeric.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)
- Tests: [`tests/test_elliptic.c`](https://github.com/stblake/mathilda/blob/main/tests/test_elliptic.c)

## Notes & additional examples

### Notes

`Pi` is the mathematical constant π. It carries the `Constant` and `Protected`
attributes, `NumericQ[Pi]` is `True`, and `D[Pi, x]` is `0`. It remains an exact
symbol through symbolic computation — driving the special-value reductions of the
trigonometric and inverse-trigonometric functions — and is evaluated to arbitrary
precision only on demand via `N[Pi, prec]`, which uses the MPFR numeric core.
