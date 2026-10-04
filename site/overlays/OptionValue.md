### Worked examples

```mathematica
In[1]:= OptionValue[NumberForm, DigitBlock]  (* the default setting for one option *)
```

```mathematica
In[1]:= OptionValue[Automatic, {a -> 2}, a]  (* explicit options, no defaults to fall back on *)
```

```mathematica
In[1]:= OptionValue[Plot, {PlotRange -> All}, PlotRange]  (* explicit setting wins over the default *)
```

### Notes

`OptionValue[f, name]` returns the setting of option `name` for `f`, taking `f`'s
default when nothing overrides it. `OptionValue[f, opts, name]` consults the
explicit `opts` first and only then the defaults from `f` (which may be
`Automatic` for none, a symbol whose `Options` are used, a rule, or a list of
such specs). Name matching ignores context prefixes.

The bare `OptionValue[name]` has no context at top level and is left unevaluated;
it is meaningful inside the right-hand side of a rule (for example an
`OptionsPattern[]` definition), where the enclosing head and its options are
injected automatically. An unknown option is left unevaluated rather than
defaulted.
