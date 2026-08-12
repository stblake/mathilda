# ComplexExpand

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ComplexExpand[expr]`**

expands expr assuming that all variables are real.

**`ComplexExpand[expr, {x1, x2, ...}]`**

expands expr assuming that variables matching any of the xi are complex; the xi may be patterns. ComplexExpand rewrites expr into explicit real and imaginary parts, propagating through Plus, Times, Power, Exp, Log, the circular and hyperbolic functions and their inverses, and the Re/Im/Abs/Arg/Conjugate/Sign/ReIm heads. The option TargetFunctions -\> {Re, Im} (default), {Abs, Arg}, or Conjugate chooses the output basis. ComplexExpand automatically threads over lists, equations, inequalities, and logic functions.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= ComplexExpand[Sin[x + I y]]
Out[1]= Sin[x] Cosh[y] + I Cos[x] Sinh[y]

In[2]:= ComplexExpand[Re[z^2], {z}]
Out[2]= -Im[z]^2 + Re[z]^2

In[3]:= ComplexExpand[Tan[x + I y]]
Out[3]= Sin[2 x]/(Cos[2 x] + Cosh[2 y]) + I Sinh[2 y]/(Cos[2 x] + Cosh[2 y])
```

### Options (1)

```mathematica
In[4]:= ComplexExpand[Re[z^2], {z}, TargetFunctions -> Conjugate]
Out[4]= 1/2 (z^2 + Conjugate[z]^2)
```

## Algorithm

```text
complex_expand.c  --  ComplexExpand
```

See complex_expand.h for the user-facing contract.

Architecture ------------ The engine is one recursive routine, cx_decompose(e, ctx, &re, &im), that

```text
writes real-valued expressions re, im with  e == re + I*im  under the
```

decomposition context ctx (which symbols are complex, and the output

```text
TargetFunctions basis).  Every other operation is a thin wrapper: the
```

builtin front-end evaluates its argument, threads over lists / relations, runs cx_decompose, and assembles Expand[re + I*im] (or pulls out one component for a Re/Im/Abs/Arg/... wrapper, which cx_decompose already handles as ordinary nodes).

TargetFunctions -> {Re, Im} and {Abs, Arg} both flow through the (re, im) engine; they differ only in how a *complex atom* is decomposed (the single

```text
substitution point cx_atom_reim()).  TargetFunctions -> Conjugate is a
```

separate, simpler path: it conjugates the whole expression (I -> -I, z -> Conjugate[z]) and averages, which reproduces the z^2/2 + Conjugate[z]^2/2 family directly.

Memory ------ All cx_* helpers BORROW their Expr* arguments and return freshly-owned,

```text
evaluated Expr*.  The builtin never frees `res` (the evaluator owns it).
```

## Implementation notes

**Attributes:** `Protected`.

## See also

[Plus](../../arithmetic/Plus/), [Times](../../arithmetic/Times/), [Power](../../arithmetic/Power/), [Abs](../../arithmetic/Abs/), [Arg](../../arithmetic/Arg/), [Exp](../../elementary-functions/Exp/), [Log](../../elementary-functions/Log/), [Re](../../arithmetic/Re/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_cherry_ei.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_ei.c)
- Tests: [`tests/test_complexexpand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_complexexpand.c)
