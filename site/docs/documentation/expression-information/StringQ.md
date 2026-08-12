# StringQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringQ[expr]`**

gives True if expr is a string, and False otherwise. The empty

<details>
<summary>Notes</summary>

string "" gives True. Called with any number of arguments other than one it leaves the expression unevaluated (StringQ::argx).

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

**Attributes:** `Protected`.

## See also

[AtomQ](../../expression-information/AtomQ/), [NumberQ](../../expression-information/NumberQ/), [IntegerQ](../../expression-information/IntegerQ/), [MachineNumberQ](../../expression-information/MachineNumberQ/), [Complex](../../arithmetic/Complex/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_pred_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_pred_compile.c)
