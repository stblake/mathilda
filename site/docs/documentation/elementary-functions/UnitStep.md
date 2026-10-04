# UnitStep

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`UnitStep[x]`**

gives 0 for x \< 0 and 1 for x \>= 0 (the value at 0 is 1).

**`UnitStep[x1, x2, ...]`**

gives 1 only when none of the xi are negative, otherwise 0.

**`UnitStep[] is 1. The result is always exact. Exact symbolic real`**

<details>
<summary>Notes</summary>

arguments are resolved by numerical certification; non-real or unresolved arguments are left unevaluated. UnitStep is Listable and Orderless.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= UnitStep[0]
Out[1]= 1

In[2]:= UnitStep[1, Pi, 5.3]
Out[2]= 1

In[3]:= UnitStep[{-1.6, 3.200000000000}]
Out[3]= {0, 1}

In[4]:= UnitStep[Sqrt[2] - 99/70]
Out[4]= 0

In[5]:= D[UnitStep[x], x]
Out[5]= Piecewise[{{Indeterminate, x == 0}}, 0]

In[6]:= D[UnitStep[x, y, z], z]
Out[6]= UnitStep[x, y] Piecewise[{{Indeterminate, z == 0}}, 0]
```

### Applications (6)

A nonnegative argument gives 1

```mathematica
In[7]:= UnitStep[3]
Out[7]= 1
```

A negative argument gives 0

```mathematica
In[8]:= UnitStep[-2]
Out[8]= 0
```

The step is closed at zero: UnitStep[0] is 1

```mathematica
In[9]:= UnitStep[0]
Out[9]= 1
```

Listable, so it threads over a vector

```mathematica
In[10]:= UnitStep[{-2, 0, 3}]
Out[10]= {0, 1, 1}
```

Several arguments: 1 only when none is negative

```mathematica
In[11]:= UnitStep[2, 3, -1]
Out[11]= 0
```

An undecidable sign is left unevaluated

```mathematica
In[12]:= UnitStep[x]
Out[12]= UnitStep[x]
```

## Options & behaviour

**Derivative** -- via the product rule, each argument contributes
`Piecewise[{{Indeterminate, xi == 0}}, 0]`:

## Implementation notes

**Algorithm.** `builtin_unitstep` classifies each argument's sign with
`ustep_class` (a numerical-certification test that decides `< 0`, `>= 0`, or
unresolved). `UnitStep[]` is 1; any argument certified negative makes the whole
call 0; arguments certified non-negative contribute a factor of 1 and are
dropped; if every argument is non-negative the result is 1, and if some remain
unresolved the call returns `UnitStep` over just those (returning `NULL`,
unevaluated, when nothing could be resolved). The result is always the exact
integer 0 or 1 when fully determined.

**Data structures.** Plain `Expr` arguments plus a small `int` class array; the
reduced call is rebuilt with `expr_new_function`. The ND kernel is **narrowing**:
`ndk_UnitStep_i` takes a `double` to an `int64` 0/1 and `ndk_UnitStep_ii` takes an
`int64` to an `int64`, with no real-closed or complex arm on purpose (the answer
is always an integer). It is on `pack.c`'s AWARE + `INT64_OK` list, so a packed or
visible numeric `NDArray` narrows to an integer buffer rather than materialising
boxed reals.

**Complexity / limits.** `O(1)` per element; `Compile[]` lowers it to `CT_INT` at
scalar and rank-1 shapes (a `double` input still yields an integer: `UnitStep[0.5]`
is `1`, not `1.`). The sign test is `UnitStep[0] = 1` (the step is closed at zero,
following `Sign[-0.] = 0`). A non-real or sign-undecidable argument stays symbolic.

- `Listable`, `NumericFunction`, `Orderless`, `Protected`.
- The result is **always exact** -- an integer `0` or `1` -- for real numeric
  input, including `Real`/`MPFR` arguments (e.g. `UnitStep[{-1.6, 3.2}]` gives
  `{0, 1}`).
- **Exact symbolic real arguments** (`Pi`, `Sqrt[2]`, `E - 3`, ...) are
  resolved by numerical certification: the argument is numericalized to MPFR at
  increasing precision and the sign is accepted only once two successive
  precisions agree on the same non-zero sign. This separates tight cases such
  as `Sqrt[2] - 99/70` ($\approx -6.4\times10^{-5}$) from zero without ever
  guessing; an argument whose sign cannot be certified is left unevaluated.
- Non-real arguments (a `Complex` with non-zero imaginary part) and unresolved
  symbolic arguments are left unevaluated. In a multidimensional call the
  proven-non-negative arguments are dropped (they contribute a factor of `1`),
  so e.g. `UnitStep[1, x]` reduces to `UnitStep[x]`.

**Attributes:** `Listable`, `NumericFunction`, `Orderless`, `Protected`.

## References

**See also:** [Orderless](../../expression-information/Orderless/), [Real](../../other-advanced/Real/), [Pi](../../mathematical-constants/Pi/), [Complex](../../arithmetic/Complex/)

- Source: [`src/piecewise.c`](https://github.com/stblake/mathilda/blob/main/src/piecewise.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_deriv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_deriv.c)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)

## Notes & additional examples

### Notes

`UnitStep` is the Heaviside step, closed at the origin (`UnitStep[0] = 1`), and
its result is always the exact integer 0 or 1 once the sign is settled. Signs are
decided by numerical certification, so an exact symbolic real like
`UnitStep[Sqrt[2] - 1]` resolves while a genuinely unknown `UnitStep[x]` stays
symbolic.

The multi-argument form is the indicator of the non-negative orthant: it drops
each argument it can prove non-negative and keeps `UnitStep` over the rest. The
`NDArray` kernel is **narrowing** — a real buffer answers with an `int64` buffer,
not boxed reals — and the same narrowing holds under `Compile[]`, where
`UnitStep[0.5]` is `1`, not `1.`.
