### Worked examples

By default `Simplify` applies its built-in transformations, collapsing the
Pythagorean identity:

```mathematica
In[1]:= Simplify[Cos[x]^2 + Sin[x]^2, TransformationFunctions -> {Automatic}]
Out[1]= 1
```

Passing an empty list disables every transformation, so the same expression is
left untouched — a direct way to see which step the built-in collection was
responsible for:

```mathematica
In[1]:= Simplify[Cos[x]^2 + Sin[x]^2, TransformationFunctions -> {}]
Out[1]= Cos[x]^2 + Sin[x]^2
```

A user-supplied transformation can stand in for the built-ins entirely: here a
single rewrite rule recovers the simplification without `Automatic`:

```mathematica
In[1]:= f = Function[e, e /. Sin[a_]^2 + Cos[a_]^2 -> 1];
In[2]:= Simplify[Cos[x]^2 + Sin[x]^2, TransformationFunctions -> {f}]
Out[2]= 1
```

Built-in transformers may also be named explicitly and combined with `Automatic`,
letting `TrigToExp` participate in the search:

```mathematica
In[1]:= Simplify[1 + Tan[x]^2, TransformationFunctions -> {Automatic, TrigToExp}]
Out[1]= Sec[x]^2
```
