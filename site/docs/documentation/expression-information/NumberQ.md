# NumberQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NumberQ[expr]
    gives True if expr is an explicit number (Integer, BigInt, Rational,
    Real, MPFR, or Complex), and False otherwise.  Symbolic constants
    such as Pi give False; use NumericQ for those.
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

`builtin_numberq` (`src/core.c`) returns `True` for an explicit number — `EXPR_INTEGER`, `EXPR_REAL`, `EXPR_BIGINT`, `EXPR_MPFR` (under `USE_MPFR`), or a `Rational`/`Complex` head — and `False` otherwise. (Contrast `NumericQ`, whose `is_numeric_quantity` helper also accepts symbolic constants like `Pi` and numeric-function calls.)

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
