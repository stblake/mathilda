### Worked examples

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_streams.txt"]; MemberQ[Streams[], str]  (* the open stream is listed *)
In[2]:= Close[str]; Streams["/tmp/mathilda_streams.txt"]  (* after Close, none for that file *)
```

### Notes

`Streams[]` returns the list of currently open `InputStream`/`OutputStream`
objects; `Streams["file"]` lists only those open for the named file. Each entry is
the same inert object `OpenRead`/`OpenWrite` handed
back, so it can be passed straight to `Read`, `Write` or
`Close`.

The list reflects the live registry, so it shrinks as streams are closed and is
empty once everything is closed (and at program start). The integer in each object
is an internal handle whose exact value depends on how many streams have been
opened in the session, so membership (`MemberQ[Streams[], str]`) is the robust way
to test for a stream rather than matching its printed form.
