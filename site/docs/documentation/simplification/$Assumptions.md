# $Assumptions

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`$Assumptions`**

is the default setting for the Assumptions option used in Simplify and other functions that take assumptions.

<details>
<summary>Notes</summary>

$Assumptions defaults to True (no assumptions). Functions like Assuming temporarily extend $Assumptions for the duration of their body.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= $Assumptions
Out[1]= True
```

### Applications (6)

The default: no assumptions in force

```mathematica
In[2]:= $Assumptions
Out[2]= True
```

Assuming extends $Assumptions for the body

```mathematica
In[3]:= Assuming[x > 0, Refine[Sqrt[x^2]]]
Out[3]= x
```

A positional assumption is conjoined with $Assumptions

```mathematica
In[4]:= Refine[Sqrt[x^2], x > 0]
Out[4]= x
```

Simplify reads the extended $Assumptions

```mathematica
In[5]:= Assuming[Element[n, Integers], Simplify[Sin[n Pi]]]
Out[5]= 0
```

A direct rebinding, which is what Assuming desugars to

```mathematica
In[6]:= Block[{$Assumptions = x > 0}, Refine[Abs[x]]]
Out[6]= x
```

Conjunction of facts, consumed per symbol

```mathematica
In[7]:= Assuming[a > 0 && b > 0, Refine[Sqrt[a^2 b^2]]]
Out[7]= a b
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
  temporarily extends `$Assumptions` for the duration of its body.

**Attributes:** none registered.

## References

**See also:** [Simplify](../../simplification/Simplify/), [Assuming](../../simplification/Assuming/)

- Source: [`src/simp/simp.c`](https://github.com/stblake/mathilda/blob/main/src/simp/simp.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)

## Notes & additional examples

### Notes

`$Assumptions` is a global symbol — the default setting for the `Assumptions` option used by
`Simplify`, `Refine`, `Element` and the other functions that take assumptions. It is not a
builtin; it carries an `OwnValue` that defaults to `True` (no assumptions). The consumers
read its value directly rather than evaluating it, so a bound fact such as
`Element[x, Reals]` does not re-fire while being read.

`Assuming[fact, body]` temporarily extends it — effectively
`Block[{$Assumptions = $Assumptions && fact}, body]` — so nested `Assuming` calls compose and
the rebinding is restored on exit. The accumulated value is flattened into a per-variable fact
set when a consumer needs it. The precedence a consumer applies is uniform: a positional
assumption is conjoined with `$Assumptions`, whereas an explicit `Assumptions -> X` option
replaces it.
