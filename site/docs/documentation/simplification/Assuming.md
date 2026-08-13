# Assuming

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Assuming[assum, expr]`**

evaluates expr with assum appended to $Assumptions, so that assum is included in the default assumptions used by functions such as Simplify.

**`Assuming[assum, expr] is effectively equivalent to Block[{$Assumptions = $Assumptions && assum}, expr], so nested invocations compose and the rebinding of $Assumptions is restored on exit.`**

<details>
<summary>Notes</summary>

Assuming converts lists of assumptions to conjunctions.

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= Assuming[x > 0, Simplify[Sqrt[x^2 y^2], y < 0]]
Out[1]= -x y
```

### Applications (5)

```mathematica
In[2]:= Simplify[Sqrt[x^2]]
Out[2]= Sqrt[x^2]

In[3]:= Assuming[x > 0, Simplify[Sqrt[x^2]]]
Out[3]= x

In[4]:= Assuming[x > 0, Simplify[Sqrt[x^2] + Abs[x]]]
Out[4]= 2 x

In[5]:= Assuming[Element[k, Integers], Simplify[Sin[k Pi]]]
Out[5]= 0

In[6]:= Assuming[a > 0 && b > 0, Simplify[Log[a b] - Log[a] - Log[b]]]
Out[6]= 0
```

## Implementation notes

**Algorithm.** `builtin_assuming` (`Assuming[assum, body]`, `ATTR_HOLDREST`)
desugars to `Block[{$Assumptions = $Assumptions && assum}, body]` and evaluates
that block. A `List` of assumptions is first normalised to an `And` conjunction
(the standard convention). Building it as `Set[$Assumptions, And[$Assumptions,
assum]]` inside the `Block` variable list reuses Block's existing scope save /
restore machinery, so `$Assumptions` is temporarily extended for the dynamic
extent of `body` and restored afterward. Nested `Assuming` calls compose
naturally because each Block reads the current `$Assumptions` value before
extending it.

**Data structures.** No state of its own — it constructs a `Block[...]` `Expr*`
and hands it to the evaluator; the assumption set lives in the `$Assumptions`
OwnValue.

- `HoldRest`, `Protected` (the assumption argument evaluates; the body is held
  until the assumption is in scope).
- Effectively `Block[{$Assumptions = $Assumptions && assum}, expr]`, so nested
  `Assuming` calls compose and the rebinding of `$Assumptions` is restored on
  exit. Lists of assumptions are converted to conjunctions.

**Attributes:** `HoldRest`, `Protected`.

## References

**See also:** [$Assumptions](../../simplification/$Assumptions/), [Simplify](../../simplification/Simplify/), [HoldRest](../../other-advanced/HoldRest/)

- Source: [`src/simp/simp_builtins.c`](https://github.com/stblake/mathilda/blob/main/src/simp/simp_builtins.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)
- Tests: [`tests/test_assuming.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assuming.c)
- Tests: [`tests/test_element.c`](https://github.com/stblake/mathilda/blob/main/tests/test_element.c)
- Tests: [`tests/test_invtrig_simplify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_invtrig_simplify.c)
- Tests: [`tests/test_limit_assumptions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_limit_assumptions.c)

## Notes & additional examples

### Notes

`Assuming[assum, expr]` evaluates `expr` with `assum` appended to `$Assumptions`, so the assumption is visible to functions such as `Simplify` and `Refine`. Lists of assumptions are combined into a conjunction. It behaves like `Block[{$Assumptions = $Assumptions && assum}, expr]`: nested invocations compose and the rebinding of `$Assumptions` is restored on exit.
