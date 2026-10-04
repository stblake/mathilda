# Options

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Options[f] gives the list of default option rules for the symbol f.`**

Options\[expr\] gives the options explicitly set in an expression such as a graphics object.  Options\[obj, name\] gives the setting for the named option; Options\[obj, {names}\] gives a list of settings.  Assign to Options\[f\] to redefine all default options at once.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Options[LinearSolve]
Out[1]= {Method -> Automatic, Modulus -> 0, ZeroTest -> Automatic}

In[2]:= Options[f] = {a -> 1, b -> 2}; f[OptionsPattern[]] := {OptionValue[a], OptionValue[b]} {f[], f[a -> 17], f[b -> 18], f[a -> 17, b -> 18]}

In[3]:= SetOptions[f, c -> 3] SetOptions::optnf: c is not a known option for f. AppendTo[Options[f], c -> 3]
Out[3]= Optional[{SetOptions::optnf (a -> 1), SetOptions::optnf (b -> 2), SetOptions::optnf (c -> 3)}, a c Dot[f, {a -> 1, b -> 2, c -> 3}] for is known not option]
```

### Applications (2)

A single named option as a rule

```mathematica
In[4]:= Options[NumberForm, DigitBlock]
Out[4]= {DigitBlock -> Infinity}
```

A symbol with no registered options has none

```mathematica
In[5]:= Options[gsym]
Out[5]= {}
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

**Algorithm.** `builtin_options` (`src/options_builtin.c`) returns an object's
option list, always as a freshly built `List` (never `NULL` for a recognised
shape). `options_of_object` branches on the first argument: a **symbol** yields a
copy of its registered defaults via `symtab_get_options` (or `{}` if it has
none); a **compound expression** yields just the option rules explicitly present
among its own arguments (each `Rule`/`RuleDelayed` whose left side is a
symbol/string, detected by `is_option_rule`); anything else yields `{}`.

The two-argument forms `Options[obj, name]` and `Options[obj, {names}]` build the
full list and then select the matching whole rules with `lookup_rule`, whose name
comparison is context-insensitive (`strip_context` drops any `` ` ``-qualified
prefix). Selected rules are returned in the requested order; a name with no
setting simply contributes nothing.

**Data structures & limits.** A symbol's defaults are stored on
`SymbolDef.default_options` as `List[Rule[name, val], ...]` — the
`DefaultValues`-equivalent — reached through `symtab_get_options`/`_set_options`.
Every element handed back is `expr_copy`'d, so the stored list is never aliased
or mutated. `Options` is `Protected`.

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

**See also:** [SetOptions](../../assignment-and-rules/SetOptions/), [OptionValue](../../assignment-and-rules/OptionValue/), [Set](../../assignment-and-rules/Set/), [Hold](../../expression-information/Hold/), [Integrate](../../calculus/Integrate/), [Limit](../../calculus/Limit/), [Series](../../power-series/Series/), [PowerExpand](../../algebra/PowerExpand/)

- Source: [`src/options_builtin.c`](https://github.com/stblake/mathilda/blob/main/src/options_builtin.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_algebraicnumbernorm.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumbernorm.c)
- Tests: [`tests/test_algebraicnumbertrace.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumbertrace.c)
- Tests: [`tests/test_array_pad.c`](https://github.com/stblake/mathilda/blob/main/tests/test_array_pad.c)

## Notes & additional examples

### Notes

`Options[s]` returns a symbol's default option settings as a list of rules
`{name -> value, ...}`; a symbol with none gives `{}`. `Options[s, name]` or
`Options[s, {names}]` selects just the requested rules, in the order asked for.

Applied to a compound expression rather than a symbol, `Options[expr]` returns
only the option rules that appear explicitly among `expr`'s arguments. Name
matching ignores any context prefix, and the returned list is a fresh copy, so
nothing stored is aliased.
