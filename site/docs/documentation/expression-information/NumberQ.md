# NumberQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NumberQ[expr]`**

gives True if expr is an explicit number (Integer, BigInt, Rational, Real, MPFR, or Complex), and False otherwise.  Symbolic constants such as Pi give False; use NumericQ for those.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= StringQ["AbC"]
Out[1]= True

In[2]:= StringQ[""]
Out[2]= True

In[3]:= StringQ[123]
Out[3]= False
```

### Scope (4)

```mathematica
In[4]:= MachineNumberQ[Sin[1000.]]
Out[4]= True
```

Overflows to +inf

```mathematica
In[5]:= MachineNumberQ[Exp[1000.]]
Out[5]= False
```

```mathematica
In[6]:= MachineNumberQ[-29037945.290347]
Out[6]= True
```

MPFR, not machine

```mathematica
In[7]:= MachineNumberQ[N[Pi, 30]]
Out[7]= False
```

## Implementation notes

`builtin_numberq` (`src/core.c`) returns `True` for an explicit number — `EXPR_INTEGER`, `EXPR_REAL`, `EXPR_BIGINT`, `EXPR_MPFR` (under `USE_MPFR`), or a `Rational`/`Complex` head — and `False` otherwise. (Contrast `NumericQ`, whose `is_numeric_quantity` helper also accepts symbolic constants like `Pi` and numeric-function calls.)

**Attributes:** `Protected`.

## References

**See also:** [AtomQ](../../expression-information/AtomQ/), [IntegerQ](../../expression-information/IntegerQ/), [StringQ](../../expression-information/StringQ/), [MachineNumberQ](../../expression-information/MachineNumberQ/), [Complex](../../arithmetic/Complex/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_bigint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bigint.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
