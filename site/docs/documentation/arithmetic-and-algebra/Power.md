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
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)

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

### Notes

Integer powers use binary exponentiation and promote to GMP bigints, so `2^200`
is exact. A rational base with a negative integer exponent inverts and raises,
giving `(1/2)^-5 = 32`. Rational exponents trigger perfect-power extraction:
`27^(2/3)` reduces to `9`, while non-extractable cases such as `8^(1/3)` of a
non-cube stay symbolic. The indeterminate form `0^0` evaluates to
`Indeterminate` rather than `1`.
