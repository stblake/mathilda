# CompileDiagnostics

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
CompileDiagnostics[argspec, expr] reports whether expr compiles for the given Compile[] argument specification, and if not, the innermost subexpression that could not be lowered. For a compiled body it also gives the result type and the instruction count with and without the optimiser.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= CompileDiagnostics[{{x, _Real}}, Sin[x] + x^2]
Out[1]= {"Compiled" -> True, "ResultType" -> "Real", "Instructions" -> 4, "CommonSubexpressions" -> 0, "InstructionsUnoptimized" -> 4}

In[2]:= CompileDiagnostics[{{x, _Real}}, Sin[x] + BarnesG[x]]
Out[2]= {"Compiled" -> False, "Reason" -> "no machine lowering for this head at these argument types", "Subexpression" -> "BarnesG[x]"}

In[3]:= CompileDiagnostics[{{x, _Real}}, Sin[x] + y]
Out[3]= {"Compiled" -> False, "Reason" -> "symbol is not a declared argument and holds no machine value", "Subexpression" -> "y"}

In[4]:= CompileDiagnostics[{{z, _Complex}}, Zeta[z]]   (* real kernel, no complex one *)
Out[4]= {"Compiled" -> False, "Reason" -> "no machine lowering for this head at these argument types", "Subexpression" -> "Zeta[z]"}
```

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/compile/compiled_function.c`](https://github.com/stblake/mathilda/blob/main/src/compile/compiled_function.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
