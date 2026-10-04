### Worked examples

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_ws.txt"];
In[2]:= WriteString[str, "a", "b", "c"]; Close[str];  (* raw: no quotes, no separators, no newline *)
In[3]:= ReadList["/tmp/mathilda_ws.txt", String]  (* the whole line is "abc" *)
```

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_ws2.txt"]; WriteString[str, "1 2 3\n4 5 6\n"]; Close[str];  (* embed newlines yourself *)
In[2]:= ReadList["/tmp/mathilda_ws2.txt", Number]
```

### Notes

`WriteString[stream, s1, s2, ...]` writes the strings **verbatim** — no surrounding
quotes, no separators between arguments, and no trailing newline. It is the tool
for laying out an exact byte stream (a CSV line, a header, a line you terminate
with an explicit `"\n"`), in contrast to `Write`, which prints input
form and appends a newline.

A non-string argument is written in input form. Output is flushed after each call.
`stream` may be an `OutputStream`, a `"file"`, or `File["file"]`; a named file that
is not already open is auto-opened (truncating) and left open. `WriteString`
returns `Null`.
