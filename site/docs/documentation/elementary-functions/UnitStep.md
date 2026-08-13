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

```mathematica
In[7]:= UnitStep[-2]
Out[7]= 0

In[8]:= UnitStep[3]
Out[8]= 1

In[9]:= Table[UnitStep[k - 3], {k, 0, 6}]
Out[9]= {0, 0, 0, 1, 1, 1, 1}

In[10]:= UnitStep[Pi - 3]
Out[10]= 1

In[11]:= UnitStep[Log[2] - Log[3]]
Out[11]= 0

In[12]:= UnitStep[1, -1, 2]
Out[12]= 0
```

## Options & behaviour

**Derivative** -- via the product rule, each argument contributes
`Piecewise[{{Indeterminate, xi == 0}}, 0]`:

## Implementation notes

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

**See also:** [Orderless](../../expression-information/Orderless/), [Pi](../../mathematical-constants/Pi/), [Complex](../../arithmetic/Complex/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_deriv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_deriv.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

`UnitStep[x]` is `0` for `x < 0` and `1` for `x >= 0` (the value at `0` is `1`). The result is always exact: certifiable real arguments resolve numerically, while non-real or unresolved arguments are left unevaluated. UnitStep is `Listable` and `Orderless`.
