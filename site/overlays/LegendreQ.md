### Worked examples

```mathematica
In[1]:= LegendreQ[0, x]  (* the second-kind solution carries a logarithm *)
```

```mathematica
In[1]:= LegendreQ[1, x]
```

```mathematica
In[1]:= LegendreQ[2, x]
```

```mathematica
In[1]:= LegendreQ[1/2, 0]  (* an exact value at the origin for non-integer order *)
```

```mathematica
In[1]:= N[LegendreQ[0, 1/2], 20]
```

```mathematica
In[1]:= N[LegendreQ[2.5, 0.4]]  (* non-integer order on the cut *)
```

### Notes

`LegendreQ[n, x]` is the Legendre function of the second kind `Q_n(x)`. Unlike
`P_n`, it carries a logarithm even at integer order: for integer `n` the result
is `P_n(x) L(x) + v_n(x)` with `L(x) = (1/2)(Log[1+x] - Log[1-x]) = ArcTanh[x]`
and a polynomial correction `v_n` (orders up to `n = 2000`).

A non-integer order is handled through two Frobenius `2F1` series on the cut
`|x| < 1`; the exact special value `Q_v(0)` is returned for an exact zero
argument. The associated functions `LegendreQ[n, m, x]` and
`LegendreQ[n, m, a, x]` (type `a` in 1, 2, 3) are obtained by differentiating
`Q_n`.
