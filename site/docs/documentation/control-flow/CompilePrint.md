# CompilePrint

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CompilePrint[cf] prints the bytecode of the CompiledFunction cf: its argument and result registers with their types, the scalar/array/tile register banks, and one line per instruction giving both the raw operands and a readable rendering. For an object whose body did not compile it reports the bail reason instead. Returns Null.`**

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= CompilePrint[Compile[{x}, Integrate[x, x]]] Signature   CompiledFunction[{x : Real}, Integrate[x, x]] Program     not compiled — every call runs the interpreter Reason      no machine lowering for this head at these argument types Bailed on   Integrate[x, x]
Out[1]= 1/2 Bailed Null Program Reason Signature argument at call compiled every for head interpreter lowering machine no not on runs the these this types u2014 x^2 CompiledFunction[{Pattern[x, Real]}, 1/2 x^2]
```

## Implementation notes

**Attributes:** `Protected`.

## See also

[CompileDiagnostics](../../control-flow/CompileDiagnostics/), [HoldAll](../../expression-information/HoldAll/), [Compile](../../control-flow/Compile/), [Attributes](../../expression-information/Attributes/)

## References

- Source: [`src/compile/compiled_function.c`](https://github.com/stblake/mathilda/blob/main/src/compile/compiled_function.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
