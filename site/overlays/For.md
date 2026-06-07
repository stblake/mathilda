### Worked examples

```mathematica
In[1]:= (s = 0; For[i = 1, i <= 5, i++, s = s + i]; s)
Out[1]= 15

In[2]:= (f = 1; For[i = 1, i <= 4, i++, f = f*i]; f)
Out[2]= 24
```

### Notes

`For[start, test, incr, body]` runs `start` once, then repeatedly evaluates `body` and `incr` while `test` is `True`. Like `Do`, it returns `Null` and is used for its side effects.
