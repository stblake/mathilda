# IntegerQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
IntegerQ[expr]
    gives True if expr is an Integer or BigInt, False otherwise.
Returns False on rationals with denominator > 1, reals, and symbolic
expressions (even those that are integer-valued, e.g. 2 Pi / Pi).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= MachineNumberQ[Sin[1000.]]
Out[1]= True

In[2]:= MachineNumberQ[Exp[1000.]]      (* overflows to +inf *)
Out[2]= False

In[3]:= MachineNumberQ[-29037945.290347]
Out[3]= True

In[4]:= MachineNumberQ[N[Pi, 30]]       (* MPFR, not machine *)
Out[4]= False

In[5]:= MachineNumberQ[1.0 + 2.0 I]
Out[5]= True

In[6]:= MachineNumberQ[1 + 2 I]         (* exact Gaussian integer *)
Out[6]= False
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
