# Fibonacci

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Fibonacci[n]
    gives the nth Fibonacci number F_n.
Fibonacci[n, x]
    gives the nth Fibonacci polynomial F_n(x).
Exact integer orders are computed via GMP fast doubling (numbers) or the
recurrence F_k = x F_{k-1} + F_{k-2} (polynomials); negative orders use
F_{-n} = (-1)^(n+1) F_n. Inexact or complex orders evaluate the
generalized closed form numerically. Listable; symbolic orders stay
unevaluated.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Table[Fibonacci[n], {n, 10}]
Out[1]= {1, 1, 2, 3, 5, 8, 13, 21, 34, 55}

In[2]:= Fibonacci[7, x]
Out[2]= 1 + 6 x^2 + 5 x^4 + x^6

In[3]:= Fibonacci[5.8, 3]
Out[3]= 283.483

In[4]:= N[Fibonacci[15/17], 50]
Out[4]= 0.956519913924311225085822634276922986486069690120617
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
