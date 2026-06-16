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
