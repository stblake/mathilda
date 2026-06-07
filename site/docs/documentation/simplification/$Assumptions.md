# $Assumptions

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
$Assumptions
    is the default setting for the Assumptions option used in Simplify and other functions that take assumptions.
$Assumptions defaults to True (no assumptions). Functions like Assuming temporarily extend $Assumptions for the duration of their body.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= $Assumptions
Out[1]= True
```

## Implementation notes

**Algorithm.** `$Assumptions` is a global assumption store, not a builtin. In
`simp_init` it is given an OwnValue defaulting to `True`
(`symtab_add_own_value("$Assumptions", ...)`). Simplify and Element read it via
`read_dollar_assumptions`, which fetches the OwnValue's replacement *directly*
(`symtab_get_own_values`) and deep-copies it rather than evaluating it — evaluating
would re-fire the rule and, for a bound `Element[x, Reals]`, cause Element to
recurse. `Assuming` extends it by desugaring to `Block[{$Assumptions =
$Assumptions && a}, body]`. The accumulated expression is parsed into an
`AssumeCtx` fact set by `assume_ctx_from_expr` (`simp_assume.c`), which flattens
`And`/`List` conjunctions and splits `Element[{...}, dom]` into per-variable facts.

**Data structures.** A single OwnValue `Rule` on the `$Assumptions` symbol in the
symbol table; consumed as an `AssumeCtx` (flat `Expr*` fact array) during
simplification.

- A system symbol with default `OwnValue` `True` (no assumptions). `Assuming`

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/simp/simp.c`](https://github.com/stblake/mathilda/blob/main/src/simp/simp.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)
