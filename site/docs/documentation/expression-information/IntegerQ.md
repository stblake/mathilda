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

## Examples (14)

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

### Applications (4)

An exact integer

```mathematica
In[11]:= IntegerQ[5]
Out[11]= True
```

A Real is not an integer, even at an integral value

```mathematica
In[12]:= IntegerQ[5.0]
Out[12]= False
```

A Rational is not an integer

```mathematica
In[13]:= IntegerQ[1/2]
Out[13]= False
```

A symbol is not known to be one, so False

```mathematica
In[14]:= IntegerQ[x]
Out[14]= False
```

## Implementation notes

`builtin_integerq` (`src/core.c`) returns `True` exactly when `expr_is_integer_like(arg)` holds (an `EXPR_INTEGER` or `EXPR_BIGINT`), and `False` otherwise.

**Attributes:** `Protected`.

## References

**See also:** [AtomQ](../../expression-information/AtomQ/), [NumberQ](../../expression-information/NumberQ/), [StringQ](../../expression-information/StringQ/), [MachineNumberQ](../../expression-information/MachineNumberQ/), [Complex](../../arithmetic/Complex/), [ExactNumberQ](../../other-advanced/ExactNumberQ/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_bigint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bigint.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_divisible.c`](https://github.com/stblake/mathilda/blob/main/tests/test_divisible.c)
- Tests: [`tests/test_meminfo.c`](https://github.com/stblake/mathilda/blob/main/tests/test_meminfo.c)

## Notes & additional examples

### Notes

`IntegerQ[e]` is `True` exactly when `e` is an exact integer — a machine `Integer` or an
arbitrary-precision bigint — and `False` for everything else, including reals at integral
values, rationals, and symbols. Like the other `*Q` predicates it always returns a
boolean; it never stays unevaluated.

Use it to test exactness: `IntegerQ[5.0]` is `False` because `5.0` is a machine real, not
an integer.
