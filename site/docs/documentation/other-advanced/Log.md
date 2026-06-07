# Log

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Log[z]
    gives the principal natural logarithm of z, with branch cut along
    the negative real axis.
Log[b, z]
    gives the logarithm to base b, i.e. Log[z] / Log[b].
Log is Listable. Log[1] = 0, Log[E] = 1, Log[E^n] = n for symbolic n.
Numeric inputs route to libm / MPFR; negative reals yield I Pi + Log[|z|].
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_log` accepts 1 arg (natural log) or 2 (`Log[b, z]`, base-b), emitting `Log::argt` otherwise. For `Log[z]`: exact special values `Log[0] = -Infinity` (exact), `Log[0.] = Indeterminate`, `Log[1] = 0`, `Log[E] = 1`, `Log[Infinity] = Infinity`; negative integers fold to `I Pi + Log[-n]` (the principal-branch shift). The key simplification is `Log[E^k] -> k` — but **only when `k` is a real numeric** (`is_real_numeric_expr`), because `E > 0` keeps us on the principal branch and an unrestricted fold would cross the branch cut for complex `k`. Numeric inputs go through MPFR (`mpfr_log` for positive reals, a complex `log|z| + I arg(z)` helper otherwise) when `USE_MPFR`, else `clog` on a `double complex`, returning a real result only when the input is real-positive. For `Log[b, z]`: `Log[b,b] = 1`; `Log[b, b^k] -> k` under the same positive-base, real-`k` branch-cut guard (`is_positive_known`); exact integer powers (`Log[2, 8] = 3`) are found by repeated division; integer-zero argument resolves the directed `±Infinity` from the sign of `Log[b]`. Everything else falls back to the rewrite `Log[b, z] -> Log[z] / Log[b]`. Unmatched inputs return NULL (left symbolic).

**Data structures.** Operates on `Expr*` trees; numeric paths use `double complex` (machine) or `mpfr_t`/the complex MPFR helpers (arbitrary precision). Negative-integer handling routes through GMP `mpz_t` to negate `INT64_MIN`/bigints safely.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/logexp.c`](https://github.com/stblake/mathilda/blob/main/src/logexp.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

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

### Notes

`Log[z]` is the principal natural logarithm; `Log[b, z]` gives the base-`b` logarithm `Log[z]/Log[b]`. Log is Listable.
