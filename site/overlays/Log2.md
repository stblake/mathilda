### Worked examples

```mathematica
In[1]:= Log2[1024]  (* exact powers of two stay exact *)
```

```mathematica
In[1]:= Log2[{1, 2, 4, 8, 16}]  (* Listable over a vector of powers *)
```

```mathematica
In[1]:= N[Log2[10], 30]  (* the base-2 logarithm of ten to thirty digits *)
```

```mathematica
In[1]:= Log2[0.1]  (* a positive machine real goes straight to libm log2 *)
```

```mathematica
In[1]:= D[Log2[x], x]  (* derivative through the Log[2, z] definition *)
```

### Notes

`Log2[z]` is `Log[2, z] = Log[z]/Log[2]`, the base-2 counterpart of `Log10`.
Information-theoretic and computer-science uses — bits of entropy, tree depth,
complexity exponents — want a base-2 logarithm, and spelling it `Log2` both reads
clearly and avoids the one-ulp drift of `Log[z]/Log[2]` at exact powers of two.

Like `Log10`, it carries an `NDArray` kernel and a `Compile[]` lowering, so a
packed buffer or a compiled body evaluates `Log2` on the buffer directly.
