# Power

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
x ^ y or Power[x, y] represents x to the power y.
Power is Listable, NumericFunction, and OneIdentity. Integer exponents
are reduced exactly (repeated squaring on GMP); Rational and Real
exponents evaluate numerically when the base is numeric; Power[0, 0]
stays Indeterminate; Power[x, 1/2] is canonicalised to Sqrt[x].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Sqrt[45]
Out[1]= 3 Sqrt[5]

In[2]:= (a * b)^2
Out[2]= a^2 b^2

In[3]:= (-1)^(3/2)
Out[3]= -I

In[4]:= (-1)^(7/4)
Out[4]= -(-1)^(3/4)

In[5]:= 18^(1/3)
Out[5]= 2^(1/3) 3^(2/3)

In[6]:= 12^(1/3)
Out[6]= 2^(2/3) 3^(1/3)

In[7]:= 60^(1/3)            (* 3 and 5 share eff 1/3 -> grouped *)
Out[7]= 2^(2/3) 15^(1/3)

In[8]:= 6^(1/3)             (* uniform exps -> stays *)
Out[8]= 6^(1/3)
```

## Implementation notes

**Algorithm.** `builtin_power` evaluates `Power[base, exp]`. `Power[x]` is x; `Power[b, e1, e2, ...]` is right-associated into `Power[b, Power[e1, e2, ...]]` (right-associative grouping). The two-argument core handles, in order: infinity/`Indeterminate` algebra (`0^Infinity -> 0`, `1^Infinity -> Indeterminate` with message, `Infinity^n` by sign of n, etc.); numeric exact folding (integer/rational/bigint powers via GMP, e.g. exact `2^10`, `(1/2)^3`); inexact Real/MPFR exponentiation; partial radical simplification of `integer^(p/q)` (pulling out perfect-power factors so `Sqrt[8] -> 2 Sqrt[2]`); `(b^m)^n -> b^(m·n)` and product/zero/one identities; and `Sqrt`-style rational-exponent canonicalisation. `Sqrt[x]` is a thin wrapper (`builtin_sqrt`) that rewrites to `Power[x, 1/2]`. Symbolic cases that cannot be reduced return `NULL`, leaving the call unevaluated. `Power` is `ONEIDENTITY | LISTABLE | NUMERICFUNCTION | PROTECTED` (note: not Flat/Orderless — exponentiation is neither associative nor commutative).

**Data structures.** `Expr*` trees; exact integer/bigint exponentiation uses GMP `mpz`, rationals via `make_rational`, and MPFR for high-precision reals. Radical factor extraction works on integer factorisation of the base.

**Complexity / limits.** Integer powers are `O(log exp)` GMP multiplies; radical canonicalisation costs a factorisation of the integer base.

- `Listable`.
- Simplifies integer powers of integers.
- Returns `Overflow[]` if the result exceeds 64-bit integer limits.
- Reduces radicals (e.g., `8^(1/2)` becomes `2*Sqrt[2]`).
- Supports complex results for negative bases (e.g., `(-1)^(1/2)` becomes `I`).

**Attributes:** `Listable`, `NumericFunction`, `OneIdentity`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on binary exponentiation.
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on simplification of radical powers.
- Source: [`src/power.c`](https://github.com/stblake/mathilda/blob/main/src/power.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= 2^200
Out[1]= 1606938044258990275541962092341162602522202993782792835301376
```

```mathematica
In[1]:= (1/2)^-5
Out[1]= 32
```

```mathematica
In[1]:= 27^(2/3)
Out[1]= 9
```

```mathematica
In[1]:= 0^0
Out[1]= Indeterminate
```

A Gaussian-integer base is raised exactly, keeping the real and imaginary parts as integers:

```mathematica
In[1]:= (3 + 4 I)^10
Out[1]= -9653287 + 1476984*I
```

A negative radicand is split into its real radical and the principal complex unit:

```mathematica
In[1]:= Sqrt[-12]
Out[1]= (2*I) Sqrt[3]
```

Irrational powers numericalise to arbitrary precision on request:

```mathematica
In[1]:= N[2^(1/2), 40]
Out[1]= 1.4142135623730950488016887242096980785697
```

### Notes

Integer powers use binary exponentiation and promote to GMP bigints, so `2^200`
is exact. A rational base with a negative integer exponent inverts and raises,
giving `(1/2)^-5 = 32`. Rational exponents trigger perfect-power extraction:
`27^(2/3)` reduces to `9`, while non-extractable cases such as `8^(1/3)` of a
non-cube stay symbolic. The indeterminate form `0^0` evaluates to
`Indeterminate` rather than `1`. Complex bases (Gaussian integers, negative
radicands) are handled in closed form, and irrational powers of numeric bases
evaluate to the requested precision under `N[...]`.
