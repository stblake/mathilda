### Worked examples

```mathematica
In[1]:= $AutoCompilation  (* on by default *)
```

```mathematica
In[1]:= $AutoCompilation = False  (* force every numeric body through the interpreter *)
```

```mathematica
In[1]:= $AutoCompilation  (* the change is visible when read back *)
```

```mathematica
In[1]:= $AutoCompilation = True  (* restore the default *)
```

### Notes

`$AutoCompilation` turns Mathilda's behind-the-scenes bytecode compilation on and
off. It covers both automatic mechanisms — the adapter that compiles a held body
once for many sample points (`Plot`, `Table`, `NIntegrate`, ...) and the
numeric-loop compiler for `Do`/`For`/`While`/`Map`/`Nest`/`Fold`/`FixedPoint` — but
not a `CompiledFunction` you built yourself with `Compile[]`. Because a compiled
body is contracted to return the interpreter's answer, the switch changes speed and
nothing else; it exists to compare the two paths. Only `True` or `False` is
accepted, and it reads back `False` when the session was started with
`MATHILDA_NO_AUTOCOMPILE` set.
