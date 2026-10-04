### Worked examples

```mathematica
In[1]:= LoadModule["simp/FullSimplify.m"]  (* load an internal module; True on success *)
```

```mathematica
In[1]:= LoadModule["simp/FullSimplify.m"]; LoadModule["simp/FullSimplify.m"]  (* already loaded: still True, not re-read *)
```

```mathematica
In[1]:= Quiet[LoadModule["no/such/module.m"]]  (* not found -> False *)
```

### Notes

`LoadModule["relpath"]` loads a Mathilda `.m` source module, with `relpath`
relative to the source tree's `src/internal` directory (e.g.
`"simp/FullSimplify.m"`). It returns `True` if the module was found and loaded (or
had already been loaded) and `False` otherwise. This is the mechanism the lazy
per-family loading in `FullSimplify` uses.

Resolution is deliberately **independent of the working directory**: it tries
`$MATHILDA_HOME`, then paths relative to the running executable (so a relocated or
installed binary still finds its bundled modules), then a compile-time prefix, then
a CWD ladder. Each module is loaded **at most once**, so repeated calls never
re-register rules — which is why the second call above is still `True` without
re-reading the file. It shares its file-reading core with `Get`.
