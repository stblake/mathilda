# ComplexExpand

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ComplexExpand[expr]`**

expands expr assuming that all variables are real.

**`ComplexExpand[expr, {x1, x2, ...}]`**

expands expr assuming that variables matching any of the xi are complex; the xi may be patterns. ComplexExpand rewrites expr into explicit real and imaginary parts, propagating through Plus, Times, Power, Exp, Log, the circular and hyperbolic functions and their inverses, and the Re/Im/Abs/Arg/Conjugate/Sign/ReIm heads. The option TargetFunctions -\> {Re, Im} (default), {Abs, Arg}, or Conjugate chooses the output basis. ComplexExpand automatically threads over lists, equations, inequalities, and logic functions.

## Examples (12)

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

### Applications (8)

Every free symbol is assumed real

```mathematica
In[5]:= ComplexExpand[Sin[x + I y]]
Out[5]= Sin[x] Cosh[y] + I Cos[x] Sinh[y]
```

Binomial over real a, b

```mathematica
In[6]:= ComplexExpand[(a + I b)^2]
Out[6]= a^2 + (2*I) a b - b^2
```

The modulus is real-valued

```mathematica
In[7]:= ComplexExpand[Abs[1 + I x]]
Out[7]= Sqrt[1 + x^2]
```

Euler's formula

```mathematica
In[8]:= ComplexExpand[Exp[I t]]
Out[8]= Cos[t] + I Sin[t]
```

With z declared complex, split via Re[z], Im[z]

```mathematica
In[9]:= ComplexExpand[z^2, {z}]
Out[9]= -Im[z]^2 + Re[z]^2 + (2*I) Im[z] Re[z]
```

A variable times its conjugate is the squared modulus

```mathematica
In[10]:= ComplexExpand[z Conjugate[z], {z}]
Out[10]= Im[z]^2 + Re[z]^2
```

Polar output basis

```mathematica
In[11]:= ComplexExpand[z^2, {z}, TargetFunctions -> {Abs, Arg}]
Out[11]= Abs[z]^2 Cos[Arg[z]]^2 - Abs[z]^2 Sin[Arg[z]]^2 + (2*I) Abs[z]^2 Cos[Arg[z]] Sin[Arg[z]]
```

Threads over the equation

```mathematica
In[12]:= ComplexExpand[Cos[x + I y] == 0]
Out[12]= Cos[x] Cosh[y] - I Sin[x] Sinh[y] == 0
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

**Algorithm.** `builtin_complex_expand` first peels any trailing
`TargetFunctions -> ...` option rule off the positionals (choosing the output
basis `{Re, Im}` by default, `{Abs, Arg}`, or `Conjugate`), then builds the
complex-variable list from the optional second positional — symbols or
*patterns*; anything not matching one of them is treated as real. The work is a
single recursive routine `cx_decompose(e, ctx, &re, &im)` that writes two
real-valued expressions with `e == re + I*im`. `Plus` sums the pairs
componentwise; `Times` folds them by complex multiplication (`cx_cmul`); `Power`
goes through `cx_power` — `E^z` and a literal `Exp[z]` use
`Exp[u](Cos[v] + I Sin[v])`, a real base with an integer exponent stays real, a
complex base with an integer exponent is raised by square-and-multiply on the
`(re, im)` pair, a positive-real base with a real exponent stays real, and
everything else uses the polar master formula `r^u Exp[-v*theta]` with
`theta = Arg[base]`. `Log[w]` becomes `(1/2) Log[Expand[u^2+v^2]] + I Arg[w]`;
the twelve circular/hyperbolic heads expand through the standard addition-formula
closed forms in `cx_apply_elem` (Tan/Cot/Tanh/Coth via double-angle rational
forms, Sec/Csc/Sech/Csch as reciprocals), and the twelve inverse heads are
rewritten to logarithmic form in `cx_inverse_to_log` and recursed.
`Re/Im/Abs/Arg/Conjugate/Sign/ReIm` are ordinary nodes; an unknown head whose
arguments are all real is assumed real, otherwise it is returned as
`Re[f[...]] + I Im[f[...]]`. A complex atom decomposes to `Re[z]`/`Im[z]` (or
`Abs[z] Cos[Arg z]` / `Abs[z] Sin[Arg z]` under `{Abs, Arg}`); every other atom
is real. `TargetFunctions -> Conjugate` is a separate path: `cx_conjugate_of`
maps `I -> -I` and each complex variable `z -> Conjugate[z]`, then averages the
expression with its conjugate. The front end threads over `List`, the relational
heads and the logical heads (collapsing `Equal`/`Unequal` to `True`/`False` via
`Simplify` when it can), and each leaf returns `Expand[re + I*im]`.

**Data structures.** Everything is `Expr` trees. Every `cx_*` helper borrows its
arguments and returns a freshly-owned, already-evaluated `Expr` built through the
`mk_*`/`ev` (= `eval_and_free`) shorthand; the builtin never frees `res` (the
evaluator owns it). Complex-variable membership is decided by running the pattern
matcher (`match`) against each supplied pattern, so `x_` and richer patterns work
as the variable spec.

**Complexity / limits.** Dominated by the repeated `eval_and_free`/`Expand` at
every node, so cost grows with the expanded size of the real and imaginary parts
(the `Times` fold is quadratic in the number of complex factors). The head is
purely symbolic and structural: no NDArray or `Compile[]` path, and `Protected`
is its only attribute. `Arg[...]`/`Abs[...]` on a symbolic complex atom are kept
inert, so an expression carrying those is as far as the decomposition reduces.

**Attributes:** `Protected`.

## References

**See also:** [Plus](../../arithmetic/Plus/), [Times](../../arithmetic/Times/), [Power](../../arithmetic/Power/), [Abs](../../arithmetic/Abs/), [Arg](../../arithmetic/Arg/), [Exp](../../elementary-functions/Exp/), [Log](../../elementary-functions/Log/), [Re](../../arithmetic/Re/)

- Source: [`src/complex_expand.c`](https://github.com/stblake/mathilda/blob/main/src/complex_expand.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_cherry_ei.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_ei.c)
- Tests: [`tests/test_complexexpand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_complexexpand.c)

## Notes & additional examples

### Notes

`ComplexExpand` rewrites an expression into explicit real and imaginary parts.
Its governing assumption is that **every free symbol is real** unless the second
argument says otherwise: `ComplexExpand[expr, {x1, x2, ...}]` treats variables
matching any `xi` as complex, and the `xi` may be patterns. A real symbol
contributes `symbol + 0 I`; a complex symbol `z` contributes `Re[z] + I Im[z]`
(the inert `Re`/`Im` heads are how an unknown complex atom is carried).

The result is always `Expand[re + I*im]`, so the imaginary unit `I` appears
explicitly and the real and imaginary parts are each fully distributed. The
decomposition propagates through `Plus`, `Times`, `Power`, `Exp`, `Log`, the
circular and hyperbolic functions and their inverses, and the
`Re`/`Im`/`Abs`/`Arg`/`Conjugate`/`Sign`/`ReIm` heads.

`TargetFunctions` chooses the output basis: `{Re, Im}` (default), `{Abs, Arg}`
(a polar form, with the complex atom written as `Abs[z] (Cos[Arg z] + I Sin[Arg
z])`), or `Conjugate` (which reports the parts in terms of `z` and
`Conjugate[z]`). `ComplexExpand` automatically threads over lists, equations,
inequalities and logical combinations, and will collapse an `Equal`/`Unequal`
leaf to `True`/`False` when `Simplify` can decide it. It is the real-variable
companion to `ExpToTrig`/`TrigToExp` and to the bare `Re`/`Im`/`Abs`/`Arg`
heads, which stay inert on a symbolic argument rather than splitting it.
