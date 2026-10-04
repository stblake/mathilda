### Worked examples

```mathematica
In[1]:= AtomQ[5]  (* integers are atoms *)
```

```mathematica
In[1]:= AtomQ[x]  (* so is a bare symbol *)
```

```mathematica
In[1]:= AtomQ[1/2]  (* a Rational is treated as atomic *)
```

```mathematica
In[1]:= AtomQ[3 + 4 I]  (* and so is a Complex *)
```

```mathematica
In[1]:= AtomQ[1 + x]  (* a Plus has parts, so it is not an atom *)
```

```mathematica
In[1]:= AtomQ[{1, 2, 3}]  (* a list is a List expression, not an atom *)
```

### Notes

`AtomQ[e]` is `True` exactly for the expressions with no subparts: integers, bigints,
reals, symbols and strings, plus the two compound heads Mathilda treats as atomic,
`Rational` and `Complex`. Every other `f[...]` is `False`.

`AtomQ` is the complement of "has a head you can take `Part`s of". It never returns
unevaluated — the answer is a structural property of the already-evaluated argument.
