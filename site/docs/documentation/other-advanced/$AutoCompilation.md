# $AutoCompilation

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`$AutoCompilation`**

controls whether Mathilda compiles numeric bodies to bytecode behind the scenes. True by default; set it to False to force every such body through the interpreter.

<details>
<summary>Notes</summary>

Covers both automatic mechanisms: the adapter that compiles a held body once for many sample points (Plot, Plot3D, Table, NIntegrate, NSum, FindRoot, the plot samplers) and the numeric-loop compiler for Do, For, While, Map, Nest, Fold and FixedPoint bodies. Compile\[\] and any CompiledFunction the user built explicitly are NOT affected -- those were asked for. A compiled body is contracted to give the interpreter's answer, so this changes speed and nothing else; it exists so the two paths can be compared. Reads back False in a session started with the environment variable MATHILDA\_NO\_AUTOCOMPILE set. Only True or False is accepted.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

On by default

```mathematica
In[1]:= $AutoCompilation
Out[1]= True
```

Force every numeric body through the interpreter

```mathematica
In[2]:= $AutoCompilation = False
Out[2]= False
```

The change is visible when read back

```mathematica
In[3]:= $AutoCompilation
Out[3]= False
```

Restore the default

```mathematica
In[4]:= $AutoCompilation = True
Out[4]= True
```

## Implementation notes

**Definition.** `$AutoCompilation` is a boolean system variable that controls
whether Mathilda compiles numeric bodies to bytecode behind the scenes. It is
`True` by default; set it to `False` to force every such body through the
interpreter. It covers **both** automatic mechanisms at once: the auto-compile
adapter that compiles a held body once for many sample points (`Plot`, `Plot3D`,
`Table`, `NIntegrate`, `NSum`, `FindRoot`, the plot samplers) and the numeric-loop
compiler for `Do`, `For`, `While`, `Map`, `Nest`, `Fold` and `FixedPoint` bodies. A
`CompiledFunction` the user built with `Compile[]` is *not* affected — that was
asked for.

**Representation.** It is one row of the `EVAL_SYSFLAGS` table in `src/eval.c`,
which pairs the name with a setter and getter. The setter
(`eval_set_autocompilation`) flips the two underlying C flags —
`autocompile_set_enabled` and `numloop_set_enabled` — so one user-facing switch
drives both paths. The value is surfaced as an OwnValue (registered by
`eval_init_sysflags` with the live default) so it can be read back as well as
assigned; the assignment hook in `apply_assignment` routes `$`-prefixed writes
through `eval_sync_sysflag`. Only `True` and `False` are accepted; any other value
raises `$AutoCompilation::flagset` and leaves the flag unchanged.

**Usage & limits.** A compiled body is contracted to give the interpreter's answer,
so this changes **speed and nothing else** — it exists so the two paths can be
compared. It reads back `False` in a session started with the environment variable
`MATHILDA_NO_AUTOCOMPILE` set (the OwnValues are registered after the environment
overrides are read, so the reported value is never a lie).

**Attributes:** none registered.

## References

- Source: [`src/eval.c`](https://github.com/stblake/mathilda/blob/main/src/eval.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

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
