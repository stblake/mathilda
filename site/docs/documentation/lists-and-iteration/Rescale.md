# Rescale

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Rescale[x, {min, max}]
    gives x rescaled to run from 0 to 1 over the range min to max, equivalent to (x - min)/(max - min).
Rescale[x, {min, max}, {ymin, ymax}]
    gives x rescaled to run from ymin to ymax over the range min to max.
Rescale[list]
    rescales each element of list to run from 0 to 1 over the range Min[list] to Max[list].
Rescale threads over a list first argument and works with exact, real, complex, and symbolic quantities.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Rescale[2.5, {-10, 10}]
Out[1]= 0.625

In[2]:= Rescale[-3/2, {-2, 2}]
Out[2]= 1/8

In[3]:= Rescale[3, {-9, 7}, {11, 28}]
Out[3]= 95/4

In[4]:= Rescale[{-2, 0, 2}]
Out[4]= {0, 1/2, 1}

In[5]:= Rescale[1 + 2 I, {0, 1 + I}]
Out[5]= 3/2 + 1/2*I
```

## Implementation notes

**Attributes:** `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
