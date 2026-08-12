# IntegerQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IntegerQ[expr]`**

gives True if expr is an Integer or BigInt, False otherwise.

<details>
<summary>Notes</summary>

Returns False on rationals with denominator \> 1, reals, and symbolic expressions (even those that are integer-valued, e.g. 2 Pi / Pi).

</details>

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

In[5]:= MachineNumberQ[Exp[1000.]]      (* overflows to +inf *)
Out[5]= False

In[6]:= MachineNumberQ[-29037945.290347]
Out[6]= True

In[7]:= MachineNumberQ[N[Pi, 30]]       (* MPFR, not machine *)
Out[7]= False
```

## Implementation notes

`builtin_integerq` (`src/core.c`) returns `True` exactly when `expr_is_integer_like(arg)` holds (an `EXPR_INTEGER` or `EXPR_BIGINT`), and `False` otherwise.

**Attributes:** `Protected`.

## See also

[AtomQ](../../expression-information/AtomQ/), [NumberQ](../../expression-information/NumberQ/), [StringQ](../../expression-information/StringQ/), [MachineNumberQ](../../expression-information/MachineNumberQ/), [Complex](../../arithmetic/Complex/)

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_bigint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bigint.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_divisible.c`](https://github.com/stblake/mathilda/blob/main/tests/test_divisible.c)
- Tests: [`tests/test_random.c`](https://github.com/stblake/mathilda/blob/main/tests/test_random.c)
