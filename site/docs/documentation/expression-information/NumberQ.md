# NumberQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NumberQ[expr]`**

gives True if expr is an explicit number (Integer, BigInt, Rational, Real, MPFR, or Complex), and False otherwise.  Symbolic constants such as Pi give False; use NumericQ for those.

## Examples (15)

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

### Applications (5)

```mathematica
In[11]:= NumberQ[3]
Out[11]= True

In[12]:= NumberQ[2/3]
Out[12]= True

In[13]:= NumberQ[1 + 2 I]
Out[13]= True
```

A symbolic constant is not an explicit number

```mathematica
In[14]:= NumberQ[Pi]
Out[14]= False
```

```mathematica
In[15]:= NumberQ[x]
Out[15]= False
```

## Implementation notes

`builtin_numberq` (`src/core.c`) returns `True` for an explicit number — `EXPR_INTEGER`, `EXPR_REAL`, `EXPR_BIGINT`, `EXPR_MPFR` (under `USE_MPFR`), or a `Rational`/`Complex` head — and `False` otherwise. (Contrast `NumericQ`, whose `is_numeric_quantity` helper also accepts symbolic constants like `Pi` and numeric-function calls.)

**Attributes:** `Protected`.

## References

**See also:** [AtomQ](../../expression-information/AtomQ/), [IntegerQ](../../expression-information/IntegerQ/), [StringQ](../../expression-information/StringQ/), [MachineNumberQ](../../expression-information/MachineNumberQ/), [Complex](../../arithmetic/Complex/), [ExactNumberQ](../../other-advanced/ExactNumberQ/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_bigint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bigint.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)

## Notes & additional examples

### Notes

`NumberQ[expr]` is `True` for an **explicit** number — an integer, bigint, machine
or arbitrary-precision real, rational, or complex. It draws the line exactly where
`NumericQ` does not: `NumberQ[Pi]` is `False` because `Pi` is a symbol that merely
*has* a numeric value, whereas `NumericQ[Pi]` is `True`. Use `NumberQ` when you
need an already-evaluated literal number, and `NumericQ` when a symbolic constant
or a numeric-function call should also qualify.
