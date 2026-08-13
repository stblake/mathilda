# SetOptions

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SetOptions[s, name -> value, ...] sets default options for the symbol`**

s and returns the new Options\[s\].  It can change Protected (but not Locked) symbols, and only changes existing options -- an unknown name raises SetOptions::optnf.  Use AppendTo\[Options\[s\], ...\] to add one.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Options[LinearSolve]
Out[1]= {Method -> Automatic, Modulus -> 0, ZeroTest -> Automatic}

In[2]:= SetOptions[f, c -> 3] SetOptions::optnf: c is not a known option for f. AppendTo[Options[f], c -> 3]
Out[2]= Optional[{SetOptions::optnf (a -> 1), SetOptions::optnf (b -> 2), SetOptions::optnf (c -> 3)}, a c Dot[f, {a -> 1, b -> 2, c -> 3}] for is known not option]
```

## Algorithm

options_builtin.c — Options, SetOptions, OptionValue and the registry of default option settings for option-accepting builtins.

Mathilda stores a symbol's default options as a List[Rule[name, val], ...] on SymbolDef.default_options (the DefaultValues-equivalent), reached through symtab_set_options / symtab_get_options. This file implements:

```text
  Options[...]      query a symbol's defaults or an expression's explicit
                    options, optionally selected by name.
  SetOptions[...]   redefine individual default options of a symbol.
  OptionValue[...]  resolve a single option value from explicit options plus
                    defaults; the bare/2-arg forms are resolved inside a rule
                    by optionvalue_inject_context (see apply_down_values).
  options_register_defaults  the comprehensive table wired from core_init.
```

Memory: every result is freshly built. Sub-expressions taken from `res` or from the stored options are duplicated with expr_copy (a refcount bump), so nothing aliased is mutated in place and the evaluator remains free to release

```text
`res` after a builtin returns.
```

## Implementation notes

- `Options`, `SetOptions`, and `OptionValue` all have attribute `{Protected}`.
- Default options survive `Clear[f]` (only rules are cleared) and are removed
  with the symbol by `Remove[f]`.
- The registered defaults mirror the options each builtin's evaluator actually
  honors, with the value it falls back to when the option is absent. The sweep
  covers, among others: `Integrate` (`Method`), `Limit` (`Direction`,
  `Assumptions`), `Series`/`PowerExpand` (`Assumptions`), `D`/`Dt`
  (`NonConstants`), `Sum`/`Product` (`Method`), `Simplify`/`FullSimplify`
  (`Assumptions`, `ComplexityFunction`, `TransformationFunctions`),
  `GroebnerBasis` (`MonomialOrder -> Lexicographic`, `CoefficientDomain ->
  Rationals`, `Method`, `Sort`, `Modulus`), `Factor`/`Together`/`Cancel`/`Apart`
  and the `Polynomial*` family (`Extension`, `Modulus`), `IrreduciblePolynomialQ`
  /`SquareFreeQ` (`GaussianIntegers`, ...), the eigen/linear-algebra heads
  (`Eigenvalues`/`Eigenvectors` → `Cubics`, `Quartics`; `LinearSolve`,
  `LeastSquares`, `NullSpace`, `MatrixRank`, `PseudoInverse`,
  `SingularValueDecomposition`), `PrimeQ`/`CoprimeQ`/`FactorInteger`
  (`GaussianIntegers`), the `Random*` heads (`WorkingPrecision`), the numerical
  calculus/root-finding heads (`NIntegrate`, `NSum`, `NProduct`, `NLimit`, `ND`,
  `NSeries`, `NResidue`, `FindRoot`, `FindMinimum`, `FindMaximum`, `NRoots`,
  `NSolve`, `Solve`, `Fit`), `Plot`, and the structural family that reads
  `Heads` (`Cases`/`Count`/`DeleteCases`/`MemberQ`/`Map`/`Apply`/`MapAll`/`Level`
  /`Depth`/`LeafCount` default `Heads -> False`; `Position`/`FreeQ` default
  `Heads -> True`).

**Attributes:** `Protected`.

## References

**See also:** [Options](../../assignment-and-rules/Options/), [OptionValue](../../assignment-and-rules/OptionValue/), [Set](../../assignment-and-rules/Set/), [Hold](../../expression-information/Hold/), [Integrate](../../calculus/Integrate/), [Limit](../../calculus/Limit/), [Series](../../power-series/Series/), [PowerExpand](../../algebra/PowerExpand/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_numberform.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numberform.c)
- Tests: [`tests/test_options.c`](https://github.com/stblake/mathilda/blob/main/tests/test_options.c)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
