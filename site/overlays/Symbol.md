### Worked examples

```mathematica
In[1]:= Symbol["x"]
```

```mathematica
In[1]:= Head[Symbol["abc"]]  (* the result is a genuine symbol *)
```

```mathematica
In[1]:= Symbol["a`x"]  (* an embedded backtick gives an absolutely-qualified name *)
```

```mathematica
In[1]:= {f[x], f["x"], f[2]} /. f[s_Symbol] :> g[s]  (* x_Symbol matches only the symbol *)
```

### Notes

`Symbol["name"]` returns the symbol with the given name, creating it if it does
not yet exist. The name must satisfy the standard symbol-name syntax: each
backtick-delimited context segment starts with a letter or `$` and continues with
letters, digits, or `$`. A leading backtick makes the name relative to the current
`$Context`, an embedded backtick gives an absolutely-qualified name, and a bare
name is resolved through `$Context` / `$ContextPath`.

Because every symbol's `Head` is `Symbol`, an `x_Symbol` pattern matches any
symbol and nothing else. An invalid name emits `Symbol::symname` and leaves the
call unevaluated; a non-string argument does likewise.
