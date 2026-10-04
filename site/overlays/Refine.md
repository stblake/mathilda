### Worked examples

```mathematica
In[1]:= Refine[Sqrt[x^2], x > 0]  (* the square root of a square collapses once the sign is known *)
```

```mathematica
In[1]:= Refine[Abs[x], x < 0]  (* Abs picks up the sign fact *)
```

```mathematica
In[1]:= Refine[Log[x^2], x > 0]  (* Log[x^p] -> p Log[x] for positive x *)
```

```mathematica
In[1]:= Refine[Sign[x], x > 0]  (* Sign resolved to +1 under the assumption *)
```

```mathematica
In[1]:= Refine[Conjugate[x], Element[x, Reals]]  (* a real symbol is its own conjugate *)
```

```mathematica
In[1]:= Refine[x > 0, x > 1]  (* an inequality decided by a Reduce/CAD entailment *)
```

```mathematica
In[1]:= Refine[Element[x, Reals], x > 0]  (* algebraic quantities in an inequality are assumed real *)
```

```mathematica
In[1]:= Refine[a - b, a == b]  (* equality substitution drives the difference to zero *)
```

```mathematica
In[1]:= Assuming[a > 0 && b > 0, Refine[Sqrt[a^2 b^2]]]  (* Refine with no positional assumption reads $Assumptions *)
```

### Notes

`Refine[expr, assum]` gives the form `expr` would take if its symbols were replaced by
explicit values satisfying `assum`. It shares its whole assumption engine with `Simplify`
(`AssumeCtx`, `apply_assumption_rules`): every rewrite `Refine` applies is one `Simplify`
also applies. The difference is that `Refine` is the rewrite-and-decide pass *without*
`Simplify`'s complexity-minimising search — so it will turn `Sqrt[x^2]` into `x` under
`x > 0`, but it does not, for instance, factor.

A predicate target (`Element[...]`, an equation, an inequality, or a logical combination) is
decided to `True`/`False` when it provably follows from or contradicts the assumptions:
domain membership through the assumption-aware prover, equations through the assumption-aware
zero test plus equality substitution, and inequalities through a `Reduce`/CAD entailment over
the reals (`P` is `True` exactly when `assum && !P` is unsatisfiable). Quantities appearing
algebraically in inequalities are assumed real.

`Refine[expr]` with no positional assumption uses the default assumptions — the current
`$Assumptions`, as extended by any enclosing `Assuming`. An explicit `Assumptions -> X`
option replaces `$Assumptions`; a positional assumption is conjoined with it. The
`TimeConstraint` option (default 30 s) caps the time spent on any single condition check
before that transformation is abandoned.
