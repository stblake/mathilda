### Worked examples

```mathematica
In[1]:= Reap[Sow[1]; Sow[2]; 99]  (* {final value, {sown values in order}} *)
```

```mathematica
In[1]:= Reap[Do[Sow[i^2], {i, 4}]]  (* Do iterates faithfully, so Sow fires each step *)
```

```mathematica
In[1]:= Reap[Sow[1, x]; Sow[2, y]; Sow[3, x], _]  (* grouped by tag, in first-seen order *)
```

```mathematica
In[1]:= Reap[42]  (* nothing sown -> an empty collection *)
```

### Notes

`Reap[expr]` evaluates `expr` and returns `{value, {sown...}}`, gathering everything `Sow`
deposits during the evaluation; `Reap[expr, patt]` keeps only matching tags, and the
third-argument form applies `f[tag, {e...}]` to each tag group. Values come back in exact
`Sow` order, grouped by tag.

One subtlety: Mathilda's `Sum` evaluates its summand **symbolically once** and closed-forms
the answer, so a `Sow` inside `Sum` fires a single time, not once per index. Use `Do`,
`Table` or an explicit loop when you want one sown value per iteration — hence the `Do`
example rather than a `Sum`.
