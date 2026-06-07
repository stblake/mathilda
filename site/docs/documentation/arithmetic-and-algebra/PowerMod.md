# PowerMod

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PowerMod[a, b, m] gives a^b mod m.
PowerMod[a, -1, m] finds the modular inverse of a modulo m.
PowerMod[a, 1/r, m] finds a modular r-th root of a.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PowerMod[2, 10, 3]
Out[1]= 1

In[2]:= PowerMod[3, -2, 7]
Out[2]= 4

In[3]:= PowerMod[3, 1/2, 2]
Out[3]= 1

In[4]:= PowerMod[2, {10, 11, 12, 13, 14}, 5]
Out[4]= {4, 3, 1, 2, 4}

In[5]:= PowerMod[100, 1/2, 17 * 19 * 23]
Out[5]= 10

In[6]:= PowerMod[2, 1/2, 10^18 + 9]
Out[6]= 742174169206529574
```

## Implementation notes

- `Protected`, `Listable`.
- Evaluates much more efficiently than `Mod[a^b, m]`.
- Integer-exponent path uses GMP `mpz_powm` / `mpz_invert`; `a`, `b`, `m`

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on modular exponentiation and modular inverses.
- von zur Gathen & Gerhard, "Modern Computer Algebra", on the extended Euclidean algorithm.
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= PowerMod[3, 1000000, 7]
Out[1]= 4
```

```mathematica
In[1]:= PowerMod[2, 100, 101]
Out[1]= 1
```

```mathematica
In[1]:= PowerMod[7, -1, 11]
Out[1]= 8
```

```mathematica
In[1]:= PowerMod[2, -1, 4]
Out[1]= PowerMod[2, -1, 4]
```

### Notes

`PowerMod[a, b, m]` computes `a^b mod m` using binary modular exponentiation, so a
huge exponent like `PowerMod[3, 1000000, 7] = 4` is evaluated without forming the
full power. `PowerMod[2, 100, 101] = 1` illustrates Fermat's little theorem, since
`101` is prime and `2^100 ≡ 1 (mod 101)`. A negative exponent computes the modular
inverse via the extended Euclidean algorithm: `PowerMod[7, -1, 11] = 8` because
`7*8 = 56 ≡ 1 (mod 11)`. When the inverse does not exist (the base shares a factor
with the modulus, as with `2` and `4`), the call returns unevaluated rather than
producing a spurious result.
