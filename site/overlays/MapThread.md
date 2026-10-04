### Worked examples

```mathematica
In[1]:= MapThread[f, {{a, b, c}, {x, y, z}}]  (* f is applied to corresponding elements *)
```

```mathematica
In[1]:= MapThread[Plus, {{1, 2, 3}, {10, 20, 30}}]  (* add two lists componentwise *)
```

```mathematica
In[1]:= MapThread[#1^#2 &, {{2, 3, 4}, {2, 2, 3}}]  (* a pure function of two arguments *)
```

```mathematica
In[1]:= MapThread[f, {{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}, 2]  (* thread two levels deep *)
```

### Notes

`MapThread[f, {l1, l2, ...}]` threads `f` over corresponding elements of the lists
`li`, so `f` receives one argument per list. It is the several-variable
generalisation of `Map`, and unlike `Thread` it takes the function and the argument
lists separately. The lists must have the same length; with the optional level `n`
they must agree in shape down through level `n`.

The pure-function form `#1^#2 &` shows the first and second threaded arguments as
`#1` and `#2`. Threading over packed arrays keeps the result packed.
