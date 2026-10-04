### Worked examples

```mathematica
In[1]:= $RaylibVerbose  (* quiet by default *)
```

```mathematica
In[1]:= $RaylibVerbose = True  (* ask the Raylib backend to print its init log *)
```

```mathematica
In[1]:= $RaylibVerbose = False  (* back to quiet *)
```

### Notes

`$RaylibVerbose` toggles the Raylib graphics backend's diagnostic log — the window
and OpenGL initialisation lines it would otherwise print when a `Show`/`Plot`
window opens or an image is exported. It is `False` by default. It affects only the
terminal chatter, never the graphics themselves, and has no visible effect in a
build compiled without graphics (`USE_GRAPHICS=0`). Only `True` or `False` is
accepted.
