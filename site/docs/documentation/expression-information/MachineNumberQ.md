# MachineNumberQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MachineNumberQ[expr] gives True if expr is a machine-precision real or complex number, and False otherwise.`**

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

## Algorithm

Mathilda — numeric evaluation implementation.

See numeric.h for the module-level overview and extensibility notes.

This file implements `N[expr]` / `N[expr, prec]`. Phase 1 targets machine-precision IEEE doubles; Phase 2 (gated behind USE_MPFR) adds MPFR arbitrary precision. Phase-2 extension points are marked with an inline "Phase 2" marker so the eventual additions are obvious.

## Implementation notes

`builtin_machinenumberq` (`src/numeric.c`) returns `True` when the argument is a finite `EXPR_REAL` (via `is_machine_real_leaf`, which checks `EXPR_REAL` and `isfinite`), or a `Complex` whose real and imaginary parts are both finite machine reals; otherwise `False`. Exact integers/rationals and arbitrary-precision `EXPR_MPFR` values are not machine numbers.

**Attributes:** `Protected`.

## References

**See also:** [AtomQ](../../expression-information/AtomQ/), [NumberQ](../../expression-information/NumberQ/), [IntegerQ](../../expression-information/IntegerQ/), [StringQ](../../expression-information/StringQ/), [Complex](../../arithmetic/Complex/)

- Source: [`src/numeric.c`](https://github.com/stblake/mathilda/blob/main/src/numeric.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_hankelmatrix.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hankelmatrix.c)
- Tests: [`tests/test_hilbertmatrix.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hilbertmatrix.c)
- Tests: [`tests/test_machine_number_q.c`](https://github.com/stblake/mathilda/blob/main/tests/test_machine_number_q.c)
- Tests: [`tests/test_toeplitzmatrix.c`](https://github.com/stblake/mathilda/blob/main/tests/test_toeplitzmatrix.c)
