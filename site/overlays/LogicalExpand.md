### Worked examples

```mathematica
In[1]:= LogicalExpand[(a || b) && c]  (* distribute And over Or into disjunctive normal form *)
```

```mathematica
In[1]:= LogicalExpand[(a || b) && (c || d)]  (* full distribution gives four conjunctive clauses *)
```

```mathematica
In[1]:= LogicalExpand[Implies[a, b]]  (* material implication becomes !a || b *)
```

```mathematica
In[1]:= LogicalExpand[Xor[a, b]]  (* exclusive-or expands to its two odd-parity clauses *)
```

```mathematica
In[1]:= LogicalExpand[Equivalent[a, b]]  (* equivalence is both-true-or-both-false *)
```

```mathematica
In[1]:= LogicalExpand[a && !a]  (* complementation collapses to False *)
```

### Notes

`LogicalExpand[expr]` rewrites a logical combination into **disjunctive normal
form** — an `Or` of `And`s. It distributes `And` over `Or`, applies De Morgan to
`Not`, and expands `Implies`, `Xor` and `Equivalent` into their `And`/`Or`
definitions, simplifying by idempotence, complementation and absorption as it
builds each clause.

A tautology collapses to `True` and a contradiction to `False`. Every non-logical
subexpression is treated as an opaque Boolean atom — there is **no domain
reasoning**, so `LogicalExpand` never decides whether `x > 0` is true — but a pair
like `x == a` and `x != a` (or `Element` and `NotElement` of the same arguments)
is recognised as complementary literals and cancels. For the full solution set of
equations and inequalities, with domain reasoning, use `Reduce`; `LogicalExpand`
is the purely structural Boolean normaliser.
