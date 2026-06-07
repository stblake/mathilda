# Root

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Root[Function[t, p[t]], k]
    Represents the k-th root of the univariate polynomial p in the
    variable t. k is canonical: real roots first ascending, then
    complex roots ordered by Re ascending, |Im| ascending, with the
    negative-Im member of each conjugate pair first. N[Root[..]]
    and N[Root[..], prec] return a numerical approximation via a
    companion-matrix + Sturm + Newton pipeline.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= N[Root[Function[#^3 - 2 # - 5], 1], 30]
Out[1]= 2.094551481542326591482386540579

In[2]:= N[Root[Function[#^3 + # + 1], 1], 20]    (* real root first *)
Out[2]= -0.682327803828019327372

In[3]:= N[Root[Function[#^3 + # + 1], 2], 20]    (* conj pair: -Im first *)
Out[3]= 0.341163901914009663686 + -1.16154139999725193609*I

In[4]:= N[Root[Function[#^3 + # + 1], 3], 20]
Out[4]= 0.341163901914009663686 + 1.16154139999725193609*I
```

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/root.c`](https://github.com/stblake/mathilda/blob/main/src/root.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
