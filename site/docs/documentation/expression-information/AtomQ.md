# AtomQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AtomQ[expr]`**

gives True if expr is an atomic object (Integer, Real, BigInt, Rational, Complex, Symbol, or String), and False if expr is a compound expression of the form head\[...\].

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringQ["AbC"]
Out[1]= True

In[2]:= StringQ[""]
Out[2]= True

In[3]:= StringQ[123]
Out[3]= False

In[4]:= StringQ[] StringQ::argx: StringQ called with 0 arguments; 1 argument is expected.
```

### Scope (6)

```mathematica
In[5]:= MachineNumberQ[Sin[1000.]]
Out[5]= True
```

Overflows to +inf

```mathematica
In[6]:= MachineNumberQ[Exp[1000.]]
Out[6]= False
```

```mathematica
In[7]:= MachineNumberQ[-29037945.290347]
Out[7]= True
```

MPFR, not machine

```mathematica
In[8]:= MachineNumberQ[N[Pi, 30]]
Out[8]= False
```

```mathematica
In[9]:= MachineNumberQ[1.0 + 2.0 I]
Out[9]= True
```

Exact Gaussian integer

```mathematica
In[10]:= MachineNumberQ[1 + 2 I]
Out[10]= False
```

## Implementation notes

`builtin_atomq` (`src/core.c`) returns `True` for any non-`EXPR_FUNCTION` node (integers, reals, bigints, symbols, strings) and for the two function heads Mathilda treats as atomic, `Rational` and `Complex`; every other `EXPR_FUNCTION` yields `False`.

**Attributes:** `Protected`.

## References

**See also:** [NumberQ](../../expression-information/NumberQ/), [IntegerQ](../../expression-information/IntegerQ/), [StringQ](../../expression-information/StringQ/), [MachineNumberQ](../../expression-information/MachineNumberQ/), [Complex](../../arithmetic/Complex/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
