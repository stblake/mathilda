# OptionValue

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`OptionValue[name] gives the value of an option named name in the`**

options matched by OptionsPattern\[\] in the enclosing rule. OptionValue\[f, name\] uses options associated with the head f; OptionValue\[f, opts, name\] resolves from the explicit rules opts then the defaults from f; a trailing Hold wraps the result in Hold.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Options[LinearSolve]
Out[1]= {Method -> Automatic, Modulus -> 0, ZeroTest -> Automatic}

In[2]:= Options[f] = {a -> 1, b -> 2}; f[OptionsPattern[]] := {OptionValue[a], OptionValue[b]} {f[], f[a -> 17], f[b -> 18], f[a -> 17, b -> 18]}

In[3]:= SetOptions[f, c -> 3] SetOptions::optnf: c is not a known option for f. AppendTo[Options[f], c -> 3]
Out[3]= Optional[{SetOptions::optnf (a -> 1), SetOptions::optnf (b -> 2), SetOptions::optnf (c -> 3)}, a c Dot[f, {a -> 1, b -> 2, c -> 3}] for is known not option]
```

### Applications (3)

The default setting for one option

```mathematica
In[4]:= OptionValue[NumberForm, DigitBlock]
Out[4]= Infinity
```

Explicit options, no defaults to fall back on

```mathematica
In[5]:= OptionValue[Automatic, {a -> 2}, a]
Out[5]= 2
```

Explicit setting wins over the default

```mathematica
In[6]:= OptionValue[Plot, {PlotRange -> All}, PlotRange]
Out[6]= All
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

**Algorithm.** `builtin_optionvalue` (`src/options_builtin.c`) resolves a single
option setting. It accepts `OptionValue[f, name]`, `OptionValue[f, opts, name]`,
and a four-argument `... , Hold` form; the bare `OptionValue[name]` has no
context at top level and is left unevaluated (it is only meaningful inside a
fired rule, where `optionvalue_inject_context` rewrites it to the three-argument
form carrying the enclosing head and its explicit options). The name must be a
symbol or string.

Resolution (`ov_resolve`) checks **explicit options first**, then defaults
derived from `f`: `lookup_value` scans the explicit `opts` (a single rule or a
list), and `defaults_lookup` interprets `f` as `Automatic` (no defaults), a
symbol (→ its `Options`), a single rule, or a list of such specs, first match
winning. All name comparisons strip any context prefix (`strip_context`) so a
symbol and its qualified string name compare equal. An unresolved option returns
`NULL` (the call stays unevaluated); the `Hold` form wraps the resolved value in
`Hold[...]`.

**Data structures & limits.** Values are read from the same
`List[Rule[name, val]]` option representation used by `Options`, and every
returned value is `expr_copy`'d. `OptionValue` is `Protected`.

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

**See also:** [Options](../../assignment-and-rules/Options/), [SetOptions](../../assignment-and-rules/SetOptions/), [Set](../../assignment-and-rules/Set/), [Hold](../../expression-information/Hold/), [Integrate](../../calculus/Integrate/), [Limit](../../calculus/Limit/), [Series](../../power-series/Series/), [PowerExpand](../../algebra/PowerExpand/)

- Source: [`src/options_builtin.c`](https://github.com/stblake/mathilda/blob/main/src/options_builtin.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)
- Tests: [`tests/test_options.c`](https://github.com/stblake/mathilda/blob/main/tests/test_options.c)
- Tests: [`tests/test_parallelmixedtower.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parallelmixedtower.c)

## Notes & additional examples

### Notes

`OptionValue[f, name]` returns the setting of option `name` for `f`, taking `f`'s
default when nothing overrides it. `OptionValue[f, opts, name]` consults the
explicit `opts` first and only then the defaults from `f` (which may be
`Automatic` for none, a symbol whose `Options` are used, a rule, or a list of
such specs). Name matching ignores context prefixes.

The bare `OptionValue[name]` has no context at top level and is left unevaluated;
it is meaningful inside the right-hand side of a rule (for example an
`OptionsPattern[]` definition), where the enclosing head and its options are
injected automatically. An unknown option is left unevaluated rather than
defaulted.
