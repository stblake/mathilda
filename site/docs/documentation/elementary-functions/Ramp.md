# Ramp

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Ramp[x]
    gives x for x >= 0 and 0 for x < 0 -- the positive part of x, and
    the standard spelling of a rectified linear unit.

The zero returned for a negative argument carries the argument's own
exactness: Ramp[-1.] is 0. and Ramp[-3] is the exact 0, so a Real
list maps to a Real list and an integer list to an integer one.
Non-real arguments, and symbolic ones whose sign cannot be decided,
are left unevaluated. Ramp is Listable and a NumericFunction.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Ramp[{-1., 0., 2.5}]
Out[1]= {0.0, 0.0, 2.5}

In[2]:= Ramp[{-3, 0, 4}]
Out[2]= {0, 0, 4}

In[3]:= Ramp[{-1/2, 3/4}]
Out[3]= {0, 3/4}

In[4]:= Ramp[1 - Sqrt[2]]
Out[4]= 0

In[5]:= Ramp[1. + 2. I]
Out[5]= Ramp[1.0 + 2.0*I]
```

## Implementation notes

- `Listable`, `NumericFunction`, `Protected`.
- The zero returned for a negative argument carries the **argument's own

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
