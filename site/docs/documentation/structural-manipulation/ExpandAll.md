# ExpandAll

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ExpandAll[expr]`**

expands out all products and integer powers in any part of expr.

**`ExpandAll[expr, patt]`**

avoids expanding parts of expr that do not contain terms matching patt.

<details>
<summary>Notes</summary>

ExpandAll effectively maps Expand and ExpandDenominator onto every part of expr, including function heads, arguments, exponents, and denominators. ExpandAll automatically threads over lists, as well as equations, inequalities, and logic functions.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= ExpandAll[1/(1+x)^3 + Sin[(1+x)^3]]
Out[1]= 1/(1 + 3 x + 3 x^2 + x^3) + Sin[1 + 3 x + 3 x^2 + x^3]

In[2]:= ExpandAll[(x+z)^2/(x+y)^2]
Out[2]= x^2/(x^2 + 2 x y + y^2) + 2 (x z)/(x^2 + 2 x y + y^2) + z^2/(x^2 + 2 x y + y^2)

In[3]:= ExpandAll[E^(I a (t-b))]
Out[3]= E^(-I a b + I a t)

In[4]:= ExpandAll[((1+a) (1+b))[x]]
Out[4]= (1 + a + b + a b)[x]

In[5]:= ExpandAll[(f[(x+y)^2] + g[(y+z)^2])^2, x]
Out[5]= f[x^2 + 2 x y + y^2]^2 + g[(y + z)^2]^2 + 2 f[x^2 + 2 x y + y^2] g[(y + z)^2]
```

## Implementation notes

- `Protected`.
- A thin recursive driver over the accelerated `Expand`: it descends into every
  part of `expr` — function heads, arguments, exponents, and the bases of
  denominators — and applies `Expand` (and `ExpandDenominator`) at each node,
  where a top-level `Expand` reaches none of them.
- `ExpandAll[expr, patt]` avoids expanding parts of `expr` that do not contain
  terms matching the pattern `patt`.
- Threads over `List`, equations, inequalities, and logic functions.
- Called with a number of arguments other than 1 or 2, emits `ExpandAll::argt`
  and stays unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [Expand](../../algebra/Expand/), [ExpandDenominator](../../algebra/ExpandDenominator/), [List](../../other-advanced/List/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_expand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expand.c)
