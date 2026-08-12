# CompileDiagnostics

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CompileDiagnostics[argspec, expr] reports whether expr compiles for the given Compile[] argument specification, and if not, the innermost subexpression that could not be lowered. Accepts the same WorkingPrecision -> n / "BigIntegers" -> True options as Compile[], so it also reports whether the arbitrary-precision subset lowers (ResultType MPFRReal/MPFRComplex/BigInteger). For a compiled body it also gives the result type and the instruction count with and without the optimiser.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## See also

[Compile](../../control-flow/Compile/), [HoldAll](../../expression-information/HoldAll/), [Plot](../../graphics/Plot/), [NIntegrate](../../numerical-calculus/NIntegrate/), [NSum](../../numerical-calculus/NSum/), [ContourPlot](../../graphics/ContourPlot/)

## References

- Source: [`src/compile/compiled_function.c`](https://github.com/stblake/mathilda/blob/main/src/compile/compiled_function.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_compile_arbprec.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_arbprec.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
- Tests: [`tests/test_compile_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_linalg.c)
- Tests: [`tests/test_compile_transforms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_transforms.c)
