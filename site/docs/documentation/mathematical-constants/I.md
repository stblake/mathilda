# I

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`I`**

is the imaginary unit Sqrt\[-1\].

<details>
<summary>Notes</summary>

I represents the imaginary unit; I^2 evaluates to -1 and complex numbers are written a + b I. It has attribute Protected, and N\[I\] is 0. + 1. I.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (10)

The defining relation: the square of the imaginary unit is -1

```mathematica
In[1]:= I^2
Out[1]= -1
```

The principal square root of -1 folds back to I

```mathematica
In[2]:= Sqrt[-1]
Out[2]= I
```

The real part of the square

```mathematica
In[3]:= Re[I^2]
Out[3]= -1
```

Euler's identity

```mathematica
In[4]:= Exp[I Pi]
Out[4]= -1
```

A quarter turn around the unit circle

```mathematica
In[5]:= Exp[I Pi/2]
Out[5]= I
```

Z times its conjugate is the squared modulus

```mathematica
In[6]:= (1 + I) (1 - I)
Out[6]= 2
```

Conjugation flips the sign of the imaginary part

```mathematica
In[7]:= Conjugate[I]
Out[7]= -I
```

Modulus one, argument a right angle

```mathematica
In[8]:= {Abs[I], Arg[I]}
Out[8]= {1, 1/2 Pi}
```

Expanding a complex power

```mathematica
In[9]:= ComplexExpand[(1 + I)^2]
Out[9]= 2*I
```

Internally the imaginary unit is a Complex atom

```mathematica
In[10]:= FullForm[I]
Out[10]= Complex[0, 1]
```

## Implementation notes

**Definition.** `I` is the imaginary unit √(−1), the number with `I^2 == -1`. Unlike
the real constants, `I` is **not** a row of the `kConstants[]` table: `core_init`
(`src/core.c`) gives it an *OwnValue* mapping the symbol `I` to `Complex[0, 1]`
(`symtab_add_own_value("I", sym_I, make_complex(0,1))`), so the evaluator rewrites `I`
to that `Complex` atom on sight. Its only attribute is `Protected`, set in the
default-attributes table in `src/attr.c` (`{"I", ATTR_PROTECTED}`); it is *not*
`Constant`, so `Attributes[I]` is `{Protected}`. Interned name `SYM_I`
(`src/sym_names.c`); docstring in `src/info.c`. Because it resolves to a number,
`I^2 -> -1`, `Sqrt[-1] -> I`, `Conjugate[I] -> -I`, `Abs[I] -> 1` and `Arg[I] -> Pi/2`
all fall out of the generic complex arithmetic (`src/complex.c`, `src/times.c`,
`src/power.c`), and `D[I, x] -> 0` because the result is constant.

**Representation & numeric value.** The surface symbol `I` evaluates to the two-field
`Complex[0, 1]` atom — `FullForm[I]` is `Complex[0, 1]` and `Head[I]` is `Complex`.
`N[I]` numericalises the two exact integer parts to machine reals, giving
`0. + 1. I`. There is **no** arbitrary-precision MPFR filler for `I` and none is
needed: because the real and imaginary parts are the exact integers `0` and `1`,
`N[I, k]` leaves them machine-valued (`0.0 + 1.0 I`) rather than inflating them to
`k` digits — precision enters only when `I` is combined with an inexact quantity.

**Usage & limits.** `I` is how all complex literals are written (`a + b I` is
`Complex` once the parts are numeric); `NumericQ[I]` is `True` (whitelisted in
`is_numeric_quantity`, `src/core.c`) and it threads through every complex-aware head —
`ComplexExpand`, `Re`/`Im`/`Conjugate`/`Abs`/`Arg`, and the Euler identities
`Exp[I Pi] -> -1`, `Exp[I Pi/2] -> I`. The only subtlety is the one above: `I` is a
fixed exact number, so asking for extra precision on `I` alone changes nothing.

- Attribute `Protected`. `Attributes[I] = {Protected}`; the symbol cannot be
  reassigned.
- Carries the OwnValue `Complex[0, 1]`, so `I` evaluates to the imaginary unit;
  `I^2 = -1` and complex numbers are written `a + b I`.
- `N[I] = 0. + 1. I`.

**Attributes:** `Protected`.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)

## Notes & additional examples

### Notes

`I` is the imaginary unit √(−1). It is not stored as a mathematical constant in its
own right: the symbol `I` carries an OwnValue rewriting it to `Complex[0, 1]`, so
`FullForm[I]` is `Complex[0, 1]` and `Head[I]` is `Complex`. Its only attribute is
`Protected` (it is *not* `Constant`), yet `NumericQ[I]` is `True` and `D[I, x]` is `0`
because it resolves to a number. Every complex identity above — `I^2 == -1`,
`Sqrt[-1] == I`, `Conjugate[I] == -I`, `Abs[I] == 1`, `Arg[I] == Pi/2`, and the Euler
relations `Exp[I Pi] == -1`, `Exp[I Pi/2] == I` — falls out of the generic complex
arithmetic. Since the real and imaginary parts are the exact integers `0` and `1`,
`N[I]` is simply `0. + 1. I`: there is no arbitrary-precision form of `I` to request.
