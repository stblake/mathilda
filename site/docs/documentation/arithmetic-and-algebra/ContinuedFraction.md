# ContinuedFraction

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ContinuedFraction[x, n]
    gives a list of the first n terms in the continued-fraction
    representation of x.
ContinuedFraction[x]
    gives all terms determinable from the precision of x.
The list {a1, a2, a3, ...} corresponds to a1 + 1/(a2 + 1/(a3 + ...)).
Exact rationals give a finite (canonical, last term >= 2) expansion.
For Sqrt[d] with d a non-square integer the no-count form returns
{a1, ..., {b1, ...}}, the bracketed block repeating cyclically. Inexact
Real / MPFR inputs yield terms only as far as the precision determines
them. ContinuedFraction is Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ContinuedFraction[47/17]
Out[1]= {2, 1, 3, 4}

In[2]:= ContinuedFraction[Sqrt[13]]
Out[2]= {3, {1, 1, 1, 1, 6}}

In[3]:= ContinuedFraction[Pi, 20]
Out[3]= {3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2}

In[4]:= ContinuedFraction[N[Pi]]
Out[4]= {3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14}

In[5]:= ContinuedFraction[Exp[Pi Sqrt[163]], 10]
Out[5]= {262537412640768743, 1, 1333462407511, 1, 8, 1, 1, 5, 1, 4}
```

## Implementation notes

**Algorithm.** `builtin_continued_fraction` computes the *simple* continued-fraction expansion, dispatching on input regime: (1) **exact rationals** use the Euclidean algorithm with floor division (`cf_rational`), producing the canonical terminating form (last term `>= 2`); (2) **`Sqrt[D]`** for a non-square positive integer `D` uses the classic periodic-surd recurrence `m'=a q-m, q'=(D-m'^2)/q, a'=floor((a0+m')/q')`, detecting the period when the `(m, q)` state first repeats (`cf_sqrt_period`), and returns `{a0, {period...}}` (or unrolls to `n` terms); (3) **inexact reals** (machine or MPFR) extract terms by repeated reciprocation while tracking absolute uncertainty (`|x| 2^-prec` plus 64 guard bits), stopping when the integer part is no longer determined by available precision (`cf_inexact`); (4) **exact symbolic reals with explicit `n`** (Pi, etc.) are numericised via `N[x, digits]` at doubling precision until `n` terms stabilise (`cf_exact_numeric`).

**Data structures.** Terms accumulate in a growable `TermVec` of GMP `mpz_t`; the MPFR path uses `mpfr_t` working registers; output is a `List` of integers (with a nested `List` for the repeating block in the unbounded surd case).

**Complexity / limits.** Rationals terminate in `O(log)` Euclidean steps. General quadratic irrationals beyond bare `Sqrt[D]` (e.g. `(1+Sqrt[5])/2`) are not recognised symbolically — supply an explicit `n` to use the numeric path. Inexact inputs stop at the precision floor.

- `Protected`, `Listable`.
- **Exact rationals** (Integer / BigInt / Rational) use the Euclidean

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/contfrac.c`](https://github.com/stblake/mathilda/blob/main/src/contfrac.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
