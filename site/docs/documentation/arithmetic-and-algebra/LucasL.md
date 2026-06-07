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

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
