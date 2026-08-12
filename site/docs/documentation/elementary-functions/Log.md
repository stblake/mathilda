# Log

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Log[z]`**

gives the principal natural logarithm of z, with branch cut along the negative real axis.

**`Log[b, z]`**

gives the logarithm to base b, i.e. Log\[z\] / Log\[b\].

<details>
<summary>Notes</summary>

Log is Listable. Log\[1\] = 0, Log\[E\] = 1, Log\[E^n\] = n for symbolic n. Numeric inputs route to libm / MPFR; negative reals yield I Pi + Log\[|z|\].

</details>

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (9)

```mathematica
In[1]:= Log[E]
Out[1]= 1

In[2]:= Log[E^2]
Out[2]= 2

In[3]:= Log[2, 8]
Out[3]= 3

In[4]:= Log[{1, E, E^2}]
Out[4]= {0, 1, 2}
```

```mathematica
In[1]:= Log[-1]
Out[1]= I Pi
```

```mathematica
In[1]:= D[Log[Sin[x]], x]
Out[1]= Cot[x]
```

```mathematica
In[1]:= Series[Log[1 + x], {x, 0, 6}]
Out[1]= x - 1/2 x^2 + 1/3 x^3 - 1/4 x^4 + 1/5 x^5 - 1/6 x^6 + O[x]^7
```

```mathematica
In[1]:= N[Log[2], 40]
Out[1]= 0.69314718055994530941723212145817656807549
```

```mathematica
In[1]:= N[Log[1 + I], 40]
Out[1]= 0.34657359027997265470861606072908828403779 + 0.78539816339744830961566084581987572104928*I
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| v^2.5 over 4x10^6 | 6.66 s | 3.79 s | 16.3 s |
| integrate Sin[x] Exp[x] | 6.62 s | 0.726 s | 5.47 s |
| Sin[Exp[Log[v]]] fused? | 6.52 s | 5.61 s | 42.5 s |
| integrate Exp[x]/x (non-elementary) | 5.38 s | 0.318 s | 33.5 s |
| Sqrt over 4x10^6 | 5.35 s | 0.444 s | 0.925 s |
| integrate Tan[x]^3 | 4.07 s | 0.196 s | 3.62 s |

## Implementation notes

**Algorithm.** `builtin_log` accepts 1 arg (natural log) or 2 (`Log[b, z]`, base-b), emitting `Log::argt` otherwise. For `Log[z]`: exact special values `Log[0] = -Infinity` (exact), `Log[0.] = Indeterminate`, `Log[1] = 0`, `Log[E] = 1`, `Log[Infinity] = Infinity`; negative integers fold to `I Pi + Log[-n]` (the principal-branch shift). The key simplification is `Log[E^k] -> k` — but **only when `k` is a real numeric** (`is_real_numeric_expr`), because `E > 0` keeps us on the principal branch and an unrestricted fold would cross the branch cut for complex `k`. Numeric inputs go through MPFR (`mpfr_log` for positive reals, a complex `log|z| + I arg(z)` helper otherwise) when `USE_MPFR`, else `clog` on a `double complex`, returning a real result only when the input is real-positive. For `Log[b, z]`: `Log[b,b] = 1`; `Log[b, b^k] -> k` under the same positive-base, real-`k` branch-cut guard (`is_positive_known`); exact integer powers (`Log[2, 8] = 3`) are found by repeated division; integer-zero argument resolves the directed `±Infinity` from the sign of `Log[b]`. Everything else falls back to the rewrite `Log[b, z] -> Log[z] / Log[b]`. Unmatched inputs return NULL (left symbolic).

**Data structures.** Operates on `Expr*` trees; numeric paths use `double complex` (machine) or `mpfr_t`/the complex MPFR helpers (arbitrary precision). Negative-integer handling routes through GMP `mpz_t` to negate `INT64_MIN`/bigints safely.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/logexp.c`](https://github.com/stblake/mathilda/blob/main/src/logexp.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_besselk.c`](https://github.com/stblake/mathilda/blob/main/tests/test_besselk.c)
- Tests: [`tests/test_bessely.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bessely.c)
- Tests: [`tests/test_cherry_dilog.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_dilog.c)

## Notes & additional examples

### Notes

`Log[z]` is the principal natural logarithm; `Log[b, z]` gives the base-`b` logarithm `Log[z]/Log[b]`. Log is Listable.
