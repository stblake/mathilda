# AtomQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
AtomQ[expr]
    gives True if expr is an atomic object (Integer, Real, BigInt,
    Rational, Complex, Symbol, or String), and False if expr is a
    compound expression of the form head[...].
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

`builtin_atomq` (`src/core.c`) returns `True` for any non-`EXPR_FUNCTION` node (integers, reals, bigints, symbols, strings) and for the two function heads Mathilda treats as atomic, `Rational` and `Complex`; every other `EXPR_FUNCTION` yields `False`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
