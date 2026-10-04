# $RaylibVerbose

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`$RaylibVerbose`**

controls whether the Raylib graphics backend prints its diagnostic log (window and OpenGL initialisation lines) to the terminal. False by default, so opening a Show/Plot window or exporting an image stays quiet; set it to True to see Raylib's full trace log.

<details>
<summary>Notes</summary>

Affects only the terminal chatter from the renderer, never the graphics produced. Has no effect in a build without graphics (USE\_GRAPHICS=0). Only True or False is accepted.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

Quiet by default

```mathematica
In[1]:= $RaylibVerbose
Out[1]= False
```

Ask the Raylib backend to print its init log

```mathematica
In[2]:= $RaylibVerbose = True
Out[2]= True
```

Back to quiet

```mathematica
In[3]:= $RaylibVerbose = False
Out[3]= False
```

## Implementation notes

**Definition.** `$RaylibVerbose` is a boolean system variable that controls whether
the Raylib graphics backend prints its diagnostic log — window and OpenGL
initialisation lines — to the terminal. It is `False` by default, so opening a
`Show`/`Plot` window or exporting an image stays quiet; set it to `True` to see
Raylib's full trace log.

**Representation.** It is one row of the `EVAL_SYSFLAGS` table in `src/eval.c`,
alongside `$AutoCompilation` and `$AutoArrayPacking`. Its setter/getter
(`raylib_verbose_set` / `raylib_verbose_enabled`) live in
`src/graphics/plot_common.c` and toggle a single raylib-free backing flag, so the
switch exists and is readable even in a build compiled without graphics. The value
is surfaced as an OwnValue by `eval_init_sysflags` so it can be read back; writes go
through the `$`-assignment hook (`eval_sync_sysflag`), which accepts only `True` or
`False` and otherwise raises `$RaylibVerbose::flagset` and keeps the current value.

**Usage & limits.** Affects only the terminal chatter from the renderer — never the
graphics produced. It has no visible effect in a build without graphics
(`USE_GRAPHICS=0`), where there is no Raylib log to gate, though the variable still
reads and assigns.

**Attributes:** none registered.

## References

- Source: [`src/eval.c`](https://github.com/stblake/mathilda/blob/main/src/eval.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`$RaylibVerbose` toggles the Raylib graphics backend's diagnostic log — the window
and OpenGL initialisation lines it would otherwise print when a `Show`/`Plot`
window opens or an image is exported. It is `False` by default. It affects only the
terminal chatter, never the graphics themselves, and has no visible effect in a
build compiled without graphics (`USE_GRAPHICS=0`). Only `True` or `False` is
accepted.
