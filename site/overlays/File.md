### Worked examples

```mathematica
In[1]:= Head[File["data.txt"]]  (* an inert wrapper head, kept unevaluated *)
```

```mathematica
In[1]:= FullForm[File["data.txt"]]  (* a one-argument File[...] holding the path string *)
```

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_word.txt"]; WriteString[str, "alpha beta gamma\n"]; Close[str];
In[2]:= ReadList[File["/tmp/mathilda_word.txt"], Word]  (* a File wrapper works like a bare path *)
```

### Notes

`File["name"]` is a symbolic wrapper for a file name, accepted wherever `Read`, `ReadList`,
`OpenRead`, `OpenWrite`, `OpenAppend` and `Close` take a file. It is an inert, `Protected`
head: `File["x"]` evaluates to itself, and `Head[File["x"]]` is `File`.

The I/O layer unwraps it on the way in — `stream_filename_arg` accepts either a bare string
or a one-argument `File[...]` and takes the path from inside — so
`ReadList[File["data.txt"], Word]` is exactly `ReadList["data.txt", Word]`. It is a name
wrapper only: it holds no handle and does not open or test the file (an open stream is an
`InputStream`/`OutputStream` object).
