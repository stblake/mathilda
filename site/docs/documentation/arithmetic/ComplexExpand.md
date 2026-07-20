# ComplexExpand

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ComplexExpand[expr]
    expands expr assuming that all variables are real.
ComplexExpand[expr, {x1, x2, ...}]
    expands expr assuming that variables matching any of the xi are complex; the xi may be patterns.
    ComplexExpand rewrites expr into explicit real and imaginary parts, propagating through Plus, Times, Power, Exp, Log, the circular and hyperbolic functions and their inverses, and the Re/Im/Abs/Arg/Conjugate/Sign/ReIm heads.
    The option TargetFunctions -> {Re, Im} (default), {Abs, Arg}, or Conjugate chooses the output basis.
    ComplexExpand automatically threads over lists, equations, inequalities, and logic functions.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ComplexExpand[Sin[x + I y]]
Out[1]= Sin[x] Cosh[y] + I Cos[x] Sinh[y]

In[2]:= ComplexExpand[Re[z^2], {z}]
Out[2]= -Im[z]^2 + Re[z]^2

In[3]:= ComplexExpand[Re[z^2], {z}, TargetFunctions -> Conjugate]
Out[3]= 1/2 (z^2 + Conjugate[z]^2)

In[4]:= ComplexExpand[Tan[x + I y]]
Out[4]= Sin[2 x]/(Cos[2 x] + Cosh[2 y]) + I Sinh[2 y]/(Cos[2 x] + Cosh[2 y])
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
