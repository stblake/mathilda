# Mod

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Mod[m, n]`**

gives the remainder of m on division by n, carrying the sign of the divisor n (a floored modulus, so Mod\[-17, 5\] is 3, not C's -2).

**`Mod[m, n, d]`**

gives the representative congruent to m modulo n in the half-open range d \<= r \< d + n.

<details>
<summary>Notes</summary>

Reduction is exact on Integer / BigInt / Rational inputs (no overflow), and works for Real / MPFR arguments at the input precision. Mod is Listable; non-numeric arguments are left unevaluated.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= Mod[-17, 5]
Out[1]= 3

In[2]:= Mod[17, -5]
Out[2]= -3

In[3]:= Mod[2^100, 7]
Out[3]= 2

In[4]:= Mod[10, 3, 1]
Out[4]= 1
```

## Implementation notes

`builtin_mod` handles the two- and three-argument forms. For Integer and BigInt operands it reduces with GMP's `mpz_fdiv_r` — a *floored* remainder, so the result carries the sign of the divisor `n`; a rational path computes `m - n*Floor[m/n]` exactly in `mpq`; MPFR operands compute `m - n*floor(m/n)` at the larger of the two input precisions; and a machine-`double` fallback covers mixed Integer/Real arguments. The three-argument `Mod[m, n, d]` reduces `m - d` and adds `d` back, landing the representative in the half-open range `[d, d+n)`. `Mod` is registered `PROTECTED | NUMERICFUNCTION | LISTABLE` and runs element-wise on a packed machine buffer. Non-numeric arguments return NULL (left symbolic).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Quotient](../../arithmetic/Quotient/), [QuotientRemainder](../../arithmetic/QuotientRemainder/), [Union](../../structural-manipulation/Union/), [Tally](../../data-structures/Tally/), [DeleteDuplicates](../../data-structures/DeleteDuplicates/)

- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on the division algorithm and modular reduction.
- von zur Gathen & Gerhard, "Modern Computer Algebra", on residues and modular arithmetic.
- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_bigint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bigint.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)

## Notes & additional examples

### Notes

`Mod[m, n]` returns a residue that takes the sign of the divisor `n`, so
`Mod[-17, 5] = 3` (non-negative) while `Mod[17, -5] = -3` (non-positive). This is
the mathematician's floored modulus, not C's truncated `%`. Reduction is exact on
bigints, so `Mod[2^100, 7] = 2` without overflow. The three-argument form
`Mod[m, n, d]` returns the representative in the offset range `[d, d+n)`, e.g.
`Mod[10, 3, 1] = 1`.
