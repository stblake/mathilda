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

- `HoldRest`, `Protected` (the assumption argument evaluates; the body is held

**Attributes:** `HoldRest`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)
