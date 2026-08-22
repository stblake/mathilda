# Hypergeometric0F1

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Hypergeometric0F1[b, z]`**

is the confluent hypergeometric 0F1, equal to HypergeometricPFQ\[{}, {b}, z\].

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= Hypergeometric0F1[1/2, z]
Out[1]= Cosh[2 Sqrt[z]]
```

### Applications (3)

```mathematica
In[2]:= Hypergeometric0F1[1/2, z]
Out[2]= Cosh[2 Sqrt[z]]

In[3]:= Hypergeometric0F1[3/2, z]
Out[3]= (1/2 Sinh[2 Sqrt[z]])/Sqrt[z]

In[4]:= N[Hypergeometric0F1[3/2, -(Pi^2/16)]*Pi/2, 40]
Out[4]= 0.99999999999999999999999999999999999999991 + 2.3165577250480442772379354523138034341839e-41*I
```

## Implementation notes

**Attributes:** `NumericFunction`, `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_hypergeopfq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergeopfq.c)

## Notes & additional examples

### Notes

`Hypergeometric0F1[b, z]` is the confluent limit `HypergeometricPFQ[{}, {b}, z]`. It converges for all `z` and underlies the Bessel functions: `0F1[1/2, z] = Cosh[2 Sqrt[z]]` and `0F1[3/2, z] = Sinh[2 Sqrt[z]]/(2 Sqrt[z])`. The tiny imaginary residue in the last result is numerical noise from the radical of a negative argument.
