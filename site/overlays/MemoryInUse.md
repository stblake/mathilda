### Worked examples

```mathematica
In[1]:= Head[MemoryInUse[]]  (* an Integer count of bytes *)
```

```mathematica
In[1]:= MemoryInUse[] > 0  (* always some memory resident *)
```

### Notes

`MemoryInUse[]` gives the number of bytes of memory currently resident for the
Mathilda process. The actual figure **changes from run to run**, so examples check
only its shape (`Head[MemoryInUse[]]` is `Integer`, `MemoryInUse[] > 0`) rather
than a literal byte count.

This is the process **resident set size**, which is not the quantity Mathematica's
`MemoryInUse[]` reports: Wolfram counts only the current session's data, whereas
this also includes the binary, the shared libraries (GMP, MPFR, LAPACK, Readline,
and Raylib on a graphics build), the stacks, and whatever the allocator holds
without returning it. On a freshly started kernel the two are not interchangeable —
but RSS is the number Activity Monitor and `top` show, which is what a notebook
status bar wants to agree with. On a platform that offers no way to ask, the call
returns unevaluated rather than reporting a misleading `0`.
