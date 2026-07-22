# Module

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Module[{x, y, ...}, expr] specifies that x, y, ... are local variables.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= x = 1; Module[{x = 2}, x + 1]
Out[1]= 3

In[2]:= x
Out[2]= 1
```

## Implementation notes

**Algorithm.** `builtin_module` (in `src/modular.c`) implements lexical scoping by **alpha-renaming** locals to unique temporaries. `Module` carries `HoldAll | Protected` (set in `src/attr.c`), so the variable list and body arrive unevaluated. The handler reads (and post-increments) the global `$ModuleNumber` counter and forms a per-invocation suffix; each local `x` (or `x = init`) becomes a fresh symbol `x$<n>` (e.g. `x$7`). Each temporary is tagged `ATTR_TEMPORARY`, and if it had an initializer (evaluated in the outer scope) that value is installed as an `OwnValue` on the renamed symbol.

The rename is applied to the body by `substitute_scoping`, a recursive tree walk over a `ScopingEnv` linked list mapping old name → replacement symbol. Crucially this walk is shadow-aware: when it descends into a *nested* scoping construct (`Module`/`Block`/`With`/`Function`/`Table`) that rebinds one of the same names, that name is dropped from the environment passed downward, so inner bindings are not corrupted; and binding RHSs are substituted with the outer environment (so `With[{q=...}, With[{k=q},...]]` resolves correctly) while binding LHS names are left intact. The renamed body is then `evaluate`d. `Return[v]` (or `Return[v, Module]`) targeting this boundary is trapped via `eval_classify_return`. Finally the `ScopingEnv`, the `VarInfo` temporaries, and the evaluated initializers are freed (the renamed symbols' OwnValues persist in the symbol table).

**Limit.** Body is taken as a single expression (arg_count must be 2).

- `HoldAll`, `Protected`.
- Variables are renamed to `name$nnn` using `$ModuleNumber`.
- Created symbols have the `Temporary` attribute.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Harold Abelson and Gerald Jay Sussman, *Structure and Interpretation of Computer Programs*, 2nd ed., §3.1 (local state and lexical scoping).
- Source: [`src/modular.c`](https://github.com/stblake/mathilda/blob/main/src/modular.c)
- Specification: [`docs/spec/builtins/scoping-constructs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/scoping-constructs.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Module[{x = 5}, x^2 + 1]
Out[1]= 26
```

```mathematica
In[1]:= Module[{a = 2, b = 3}, a*b + a + b]
Out[1]= 11
```

```mathematica
In[1]:= f[n_] := Module[{s = 0}, s = n^2 + n; s]; f[4]
Out[1]= 20
```

```mathematica
In[1]:= g[n_] := Module[{f}, f[0] = 1; f[k_] := k*f[k - 1]; f[n]]; g[6]
Out[1]= 720
```

```mathematica
In[1]:= Module[{x = 1}, Do[x = x + 1/x, {5}]; x]
Out[1]= 969581/272890
```

### Notes

`Module[{vars}, body]` introduces lexically scoped local variables, optionally
with initial values (`{x = 5}`). The locals are renamed to unique symbols so they
never collide with global definitions of the same name. Inside the body, locals
can be reassigned (`s = n^2 + n`) and a sequence of statements separated by `;`
is evaluated left to right, with the last expression returned. This makes
`Module` the standard tool for writing multi-step function definitions. The local
symbols can even carry their own `DownValues`: in the factorial example a *local*
`f` is given a base case and a recursive rule, so the recursion runs entirely
inside the module without polluting any global `f`. The continued-fraction
example iterates `x -> x + 1/x` five times with `Do`, returning the exact rational
result thanks to arbitrary-precision arithmetic.
