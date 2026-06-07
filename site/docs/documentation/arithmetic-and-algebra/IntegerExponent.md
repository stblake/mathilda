# IntegerExponent

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
IntegerExponent[n, b] gives the highest power of b that divides n.
IntegerExponent[n] is equivalent to IntegerExponent[n, 10] and gives the number of trailing zeros in the decimal digits of n.
IntegerExponent ignores the sign of n; IntegerExponent[0, b] is Infinity.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= IntegerExponent[1230000]
Out[1]= 4

In[2]:= IntegerExponent[2^10 + 2^7, 2]
Out[2]= 7

In[3]:= IntegerExponent[144, 2]
Out[3]= 4

In[4]:= IntegerExponent[100!, 2]
Out[4]= 97

In[5]:= IntegerExponent[0]
Out[5]= Infinity
```

```mathematica
In[1]:= IntegerExponent[]
Out[1]= IntegerExponent[]

In[2]:= IntegerExponent[1, 2, 3, 4]
Out[2]= IntegerExponent[1, 2, 3, 4]
```

```mathematica
In[1]:= IntegerExponent[1.123]
Out[1]= IntegerExponent[1.123]
```

## Implementation notes

- `Protected`, `Listable`. Threads element-wise over a list of integers in

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
