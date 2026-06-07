### Worked examples

```mathematica
In[1]:= ReplaceRepeated[x, {x -> y, y -> z}]
Out[1]= z

In[2]:= (1 + 0) + (a*1) //. {x_ + 0 -> x, x_*1 -> x}
Out[2]= 1 + a

In[3]:= Sin[x]^2 + Cos[x]^2 //. Sin[a_]^2 + Cos[a_]^2 -> 1
Out[3]= 1
```

### Notes

`ReplaceRepeated` (infix `//.`) reapplies the rules until the expression stops changing, whereas `ReplaceAll` (`/.`) makes a single pass. It is the natural choice for simplification rule sets where one rewrite exposes the next.
