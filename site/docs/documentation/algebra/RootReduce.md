# RootReduce

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
RootReduce[expr] canonicalises an algebraic expression: a constant algebraic number becomes a rational, a quadratic radical, or a Root object; a rational function over a radical tower has its denominator rationalised. Threads over lists, equations, inequalities and logic. Option: Method -> "Automatic" | "Recursive" | "NumberField".
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= RootReduce[Sqrt[2] + Sqrt[3]]
Out[1]= Root[1 - 10 #1^2 + #1^4 &, 4]

In[2]:= RootReduce[(Sqrt[18] + Sqrt[27]) / Sqrt[5 + 2 Sqrt[6]]]
Out[2]= 3

In[3]:= RootReduce[1/(1 + Sqrt[2])]
Out[3]= -1 + Sqrt[2]

In[4]:= RootReduce[1/(1 + 2^(1/3) + 2^(2/3))]
Out[4]= Root[-1 + 3 #1 + 3 #1^2 + #1^3 &, 1]

In[5]:= RootReduce[1/(1 + k^(1/3))]        (* parametric tower *)
Out[5]= (1 - k^(1/3) + k^(2/3))/(1 + k)
```

## Implementation notes

- `Protected`, `Listable`. Threads over lists, and over equations, inequalities

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/rootreduce.c`](https://github.com/stblake/mathilda/blob/main/src/rootreduce.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
