# Assuming

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Assuming[assum, expr]
    evaluates expr with assum appended to $Assumptions, so that assum is included in the default assumptions used by functions such as Simplify.
Assuming converts lists of assumptions to conjunctions.
Assuming[assum, expr] is effectively equivalent to Block[{$Assumptions = $Assumptions && assum}, expr], so nested invocations compose and the rebinding of $Assumptions is restored on exit.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Assuming[x > 0, Simplify[Sqrt[x^2 y^2], y < 0]]
Out[1]= -x y
```

## Implementation notes

**Algorithm.** `builtin_assuming` (`Assuming[assum, body]`, `ATTR_HOLDREST`)
desugars to `Block[{$Assumptions = $Assumptions && assum}, body]` and evaluates
that block. A `List` of assumptions is first normalised to an `And` conjunction
(matching Mathematica). Building it as `Set[$Assumptions, And[$Assumptions,
assum]]` inside the `Block` variable list reuses Block's existing scope save /
restore machinery, so `$Assumptions` is temporarily extended for the dynamic
extent of `body` and restored afterward. Nested `Assuming` calls compose
naturally because each Block reads the current `$Assumptions` value before
extending it.

**Data structures.** No state of its own — it constructs a `Block[...]` `Expr*`
and hands it to the evaluator; the assumption set lives in the `$Assumptions`
OwnValue.

- `HoldRest`, `Protected` (the assumption argument evaluates; the body is held

**Attributes:** `HoldRest`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/simp/simp_builtins.c`](https://github.com/stblake/mathilda/blob/main/src/simp/simp_builtins.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)

## Notes & additional examples

### Worked examples

Without assumptions `Sqrt[x^2]` cannot be reduced; supplying `x > 0` resolves it:

```mathematica
In[1]:= Simplify[Sqrt[x^2]]
Out[1]= Sqrt[x^2]

In[2]:= Assuming[x > 0, Simplify[Sqrt[x^2]]]
Out[2]= x

In[3]:= Assuming[x > 0, Simplify[Sqrt[x^2] + Abs[x]]]
Out[3]= 2 x
```

Domain assumptions feed Simplify's decision procedures; integer `k` kills the sine, positive `a, b` collapse the logarithm:

```mathematica
In[1]:= Assuming[Element[k, Integers], Simplify[Sin[k Pi]]]
Out[1]= 0

In[2]:= Assuming[a > 0 && b > 0, Simplify[Log[a b] - Log[a] - Log[b]]]
Out[2]= 0
```

### Notes

`Assuming[assum, expr]` evaluates `expr` with `assum` appended to `$Assumptions`, so the assumption is visible to functions such as `Simplify` and `Refine`. Lists of assumptions are combined into a conjunction. It behaves like `Block[{$Assumptions = $Assumptions && assum}, expr]`: nested invocations compose and the rebinding of `$Assumptions` is restored on exit.
