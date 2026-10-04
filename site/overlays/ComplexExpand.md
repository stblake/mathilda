### Worked examples

```mathematica
In[1]:= ComplexExpand[Sin[x + I y]]  (* every free symbol is assumed real *)
```

```mathematica
In[1]:= ComplexExpand[(a + I b)^2]  (* binomial over real a, b *)
```

```mathematica
In[1]:= ComplexExpand[Abs[1 + I x]]  (* the modulus is real-valued *)
```

```mathematica
In[1]:= ComplexExpand[Exp[I t]]  (* Euler's formula *)
```

```mathematica
In[1]:= ComplexExpand[z^2, {z}]  (* with z declared complex, split via Re[z], Im[z] *)
```

```mathematica
In[1]:= ComplexExpand[z Conjugate[z], {z}]  (* a variable times its conjugate is the squared modulus *)
```

```mathematica
In[1]:= ComplexExpand[z^2, {z}, TargetFunctions -> {Abs, Arg}]  (* polar output basis *)
```

```mathematica
In[1]:= ComplexExpand[Cos[x + I y] == 0]  (* threads over the equation *)
```

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
