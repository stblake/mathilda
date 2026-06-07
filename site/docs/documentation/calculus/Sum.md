# Sum

!!! warning "Status: Partial"
    implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## Description

```text
Sum[f, {i, imax}] gives the sum of f for i from 1 to imax. Sum[f, {i, imin, imax}], Sum[f, {i, imin, imax, di}] and Sum[f, {i, {i1, i2, ...}}] use the standard iterator forms; multiple iterators give nested sums. Sum[f, i] gives the indefinite sum (antidifference). Symbolic and infinite sums are evaluated in closed form via Method -> "Polynomial" | "Geometric" | "Gosper".
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Sum[i^2, {i, 1, 100}]
Out[1]= 338350

In[2]:= Sum[i^2, {i, 1, n}]
Out[2]= 1/6 n (1 + n) (1 + 2 n)

In[3]:= Sum[f[i, j], {i, 1, 3}, {j, 1, i}]
Out[3]= f[1, 1] + f[2, 1] + f[2, 2] + f[3, 1] + f[3, 2] + f[3, 3]
```

```mathematica
In[1]:= Sum[i^3, {i, 1, n}]
Out[1]= 1/4 n^2 (1 + n)^2

In[2]:= Sum[i^2, i]
Out[2]= 1/6 i (-1 + i) (-1 + 2 i)
```

```mathematica
In[1]:= Sum[a^i, i]
Out[1]= a^i/(-1 + a)

In[2]:= Sum[q1^i q2^i, i]
Out[2]= (q1 q2)^i/(-1 + q1 q2)
```

```mathematica
In[1]:= Sum[k k!, k]
Out[1]= Factorial[k]
```

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Partial** — implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## References

- Petkovšek, Wilf & Zeilberger, "A=B" (A K Peters, 1996).
- Graham, Knuth & Patashnik, "Concrete Mathematics", 2nd ed. (Addison-Wesley, 1994), ch. 2 & 6.
- Source: [`src/sum/sum.c`](https://github.com/stblake/mathilda/blob/main/src/sum/sum.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Sum[k, {k, 1, 10}]
Out[1]= 55
```

```mathematica
In[1]:= Sum[k^2, {k, 1, n}]
Out[1]= 1/6 n (1 + n) (1 + 2 n)
```

```mathematica
In[1]:= Sum[k^3, {k, 1, n}]
Out[1]= 1/4 n^2 (1 + n)^2
```

```mathematica
In[1]:= Sum[2^k, {k, 0, 5}]
Out[1]= 63
```

### Notes

`Sum` evaluates numeric ranges directly and closes symbolic finite ranges in closed form through the polynomial, geometric, and Gosper (`Method`) routines — so `Sum[k^2, {k, 1, n}]` returns Faulhaber's polynomial. The Gosper backend handles hypergeometric summands over a symbolic upper bound. Infinite sums are **not** evaluated: `Sum[1/k^2, {k, 1, Infinity}]` and `Sum[1/2^k, {k, 0, Infinity}]` both return unevaluated, so convergent series like ζ(2) and geometric series stay symbolic. `Sum` is `HoldAll`, so the iterator variable is not evaluated before the range is set up.
