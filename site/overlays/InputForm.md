### Worked examples

```mathematica
In[1]:= InputForm[1/2]
Out[1]= 1/2

In[2]:= InputForm[a + b]
Out[2]= a + b

In[3]:= InputForm[{1, 2, 3}]
Out[3]= {1, 2, 3}
```

```mathematica
In[1]:= InputForm[1/2 + 3/4 I]
Out[1]= 1/2 + 3/4*I
```

```mathematica
In[1]:= InputForm[{1, 1/2, "a", x}]
Out[1]= {1, 1/2, "a", x}
```

### Notes

`InputForm` prints an expression in a form the parser can read back in, unlike `FullForm` which exposes the internal tree. It is the form to use when you need to copy a result back into the REPL or store it as text.
