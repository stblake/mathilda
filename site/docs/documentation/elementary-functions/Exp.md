# Exp

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Exp[z]
    gives the exponential E^z.
Exp is Listable. Exp[0] = 1, Exp[Log[x]] -> x, Exp[I Pi] = -1.
Numeric inputs route to libm / MPFR.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_exp` (1-arg). Exact special values: `Exp[0] = 1`, `Exp[-Infinity] = 0`, `Exp[Infinity] = Infinity`. The notable closed-form path is **Euler's formula** for `Exp[I q Pi]`: when the argument is `Times[Complex[0,q], Pi]` with `q` an integer or rational, it rewrites to `Cos[q Pi] + I Sin[q Pi]` (which the trig kernels then evaluate to exact algebraic values). Numeric arguments go through MPFR when `USE_MPFR` — `mpfr_exp` for pure reals, the `exp(a)(cos b + i sin b)` complex helper otherwise — or `cexp` on a `double complex` for machine inputs, collapsing to a real result when the imaginary part is zero. Anything that matches none of these is returned as `Power[E, z]` (the canonical internal form for an unevaluated exponential), so `Exp` is effectively a thin front-end that produces `E^z` and lets the `Power` machinery carry the symbolic case.

**Data structures.** `Expr*` trees throughout; numeric evaluation uses `double complex` (machine) or `mpfr_t`/complex-MPFR helpers (arbitrary precision). The Euler path scans the `Times` argument list for a `Pi` factor and a single pure-imaginary `Complex[0, q]` coefficient.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/logexp.c`](https://github.com/stblake/mathilda/blob/main/src/logexp.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Exp[0]
Out[1]= 1

In[2]:= Exp[I Pi]
Out[2]= -1

In[3]:= Exp[Log[x]]
Out[3]= x
```

```mathematica
In[1]:= Exp[2 Log[x]]
Out[1]= x^2
```

```mathematica
In[1]:= D[Exp[Sin[x]], x]
Out[1]= Cos[x] E^Sin[x]
```

```mathematica
In[1]:= Series[Exp[x], {x, 0, 6}]
Out[1]= 1 + x + 1/2 x^2 + 1/6 x^3 + 1/24 x^4 + 1/120 x^5 + 1/720 x^6 + O[x]^7
```

```mathematica
In[1]:= N[Exp[1], 50]
Out[1]= 2.71828182845904523536028747135266249775724709369996
```

### Notes

`Exp[z]` is `E^z`; it inverts `Log` and evaluates Euler's identity exactly.
Logarithmic arguments collapse (`Exp[2 Log[x]] -> x^2`), it differentiates and
series-expands symbolically, and numeric inputs route to libm / MPFR for
arbitrary precision. Exp is Listable.
