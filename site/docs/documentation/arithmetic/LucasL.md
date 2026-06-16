# LucasL

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
LucasL[n]
    gives the nth Lucas number L_n.
LucasL[n, x]
    gives the nth Lucas polynomial L_n(x).
Exact integer orders are computed via GMP fast doubling (numbers, using
L_m = 2 F_{m+1} - F_m) or the recurrence L_k = x L_{k-1} + L_{k-2} with
L_0 = 2, L_1 = x (polynomials); negative orders use L_{-n} = (-1)^n L_n.
Inexact or complex orders evaluate the generalized closed form
L_n = phi^n + Cos[Pi n] phi^-n (phi = GoldenRatio) numerically.
Listable; symbolic orders stay unevaluated.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Table[LucasL[n], {n, 10}]
Out[1]= {1, 3, 4, 7, 11, 18, 29, 47, 76, 123}

In[2]:= LucasL[7, x]
Out[2]= 7 x + 14 x^3 + 7 x^5 + x^7

In[3]:= LucasL[-11.]
Out[3]= -199.0

In[4]:= N[LucasL[11/3], 50]
Out[4]= 5.92396265296195541013569786219401262875198554223617
```

## Implementation notes

**Algorithm.** `builtin_lucasl` mirrors `Fibonacci`. For exact integer order it fast-doubles the Fibonacci pair `(F_m, F_{m+1})` in GMP and derives `L_m = 2 F_{m+1} - F_m` in `O(log n)`, with negative orders via `L_{-m} = (-1)^m L_m`. `LucasL[n, x]` (Lucas polynomial) iterates `L_k = x L_{k-1} + L_{k-2}` from `L_0 = 2, L_1 = x`, Expand-ing each step. For inexact/non-integer order it builds the closed form `phi^n + Cos[Pi n] phi^-n` (`phi = GoldenRatio`, or `beta = (x+Sqrt[x^2+4])/2`) and hands it to `numericalize`. Purely symbolic order returns `NULL`.

**Data structures.** GMP `mpz_t` integer pair; `Expr` trees through `eval_and_free` for the polynomial and closed-form branches.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/lucas.c`](https://github.com/stblake/mathilda/blob/main/src/lucas.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= LucasL[10]
Out[1]= 123
```

```mathematica
In[1]:= LucasL[5, x]
Out[1]= 5 x + 5 x^3 + x^5
```

```mathematica
In[1]:= LucasL[100]
Out[1]= 792070839848372253127
```

```mathematica
In[1]:= LucasL[-7]
Out[1]= -29
```

```mathematica
In[1]:= N[LucasL[20]/LucasL[19], 20]
Out[1]= 1.6180339887498948482
```

### Notes

`LucasL[n]` is the `n`th Lucas number `L_n` (`L_0 = 2`, `L_1 = 1`, `L_k = L_{k-1} + L_{k-2}`); `LucasL[n, x]` is the Lucas polynomial `L_n(x)`. Integer orders use GMP fast doubling for arbitrary size (so `LucasL[100]` is exact), and negative orders follow `L_{-n} = (-1)^n L_n`. Consecutive ratios `L_{n+1}/L_n` converge to `GoldenRatio` (last example). Inexact or complex orders use the closed form `phi^n + Cos[Pi n] phi^-n`. Listable; symbolic orders stay unevaluated.
