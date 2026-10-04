### Worked examples

```mathematica
In[1]:= n = 1; While[n < 4, n++]; n  (* repeat the body while the test stays True *)
```

```mathematica
In[1]:= m = 0; While[m < 100, m += 7]; m  (* read the running value out afterward *)
```

```mathematica
In[1]:= gcd2[a0_, b0_] := Module[{a = a0, b = b0, t}, While[b != 0, t = b; b = Mod[a, b]; a = t]; a]; gcd2[48, 36]  (* the Euclidean algorithm as a While loop *)
```

### Notes

`While[test, body]` evaluates `test`, then `body`, repeatedly, until `test` first
fails to give `True`; `While[test]` runs an empty body, useful when `test` itself
has the side effect. Both arguments are held (`HoldAll`) and re-evaluated each
pass.

`While` returns `Null`, so — like `Do` — it is used for its side effects and the
result is read out of a variable afterward. `Break[]` exits the loop, `Continue[]`
skips to the next test, and `Return[v]` makes the loop yield `v`; `Throw`, `Abort`
and `Quit` propagate unchanged. If the first test is not `True`, the body never
runs.
