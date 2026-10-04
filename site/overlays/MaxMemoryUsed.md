### Worked examples

```mathematica
In[1]:= Head[MaxMemoryUsed[]]  (* an Integer count of bytes *)
```

```mathematica
In[1]:= MaxMemoryUsed[] > 0  (* a positive high-water mark *)
```

### Notes

`MaxMemoryUsed[]` gives the peak number of bytes resident for the Mathilda process
over its lifetime. The actual figure **changes from run to run**, so examples check
only its shape (`Head[MaxMemoryUsed[]]` is `Integer`, `MaxMemoryUsed[] > 0`) rather
than a literal byte count.

The value is a genuine high-water mark from the operating system (`getrusage`'s
`ru_maxrss`), not the largest figure some earlier `MemoryInUse[]` call happened to
observe. That distinction is the point: a polled maximum would miss any spike
falling between two polls, and a status bar polling once a second would miss nearly
every spike worth knowing about. Because the peak and the current figure are read
from two different OS counters (`getrusage` versus the per-task resident size),
`MaxMemoryUsed[]` and `MemoryInUse[]` are not guaranteed to agree to the byte at a
given instant — near startup the peak can even read a page below the current RSS.
