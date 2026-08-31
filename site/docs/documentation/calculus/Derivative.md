# Derivative

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Derivative[n1, n2, ...][f] is the general form, representing a function`**

<details>
<summary>Notes</summary>

f' represents the derivative of a function f of one argument. obtained from f by differentiating n1 times with respect to the first argument, n2 times with respect to the second argument, and so on. f' is equivalent to Derivative\[1\]\[f\]; f'' evaluates to Derivative\[2\]\[f\]. Derivative is a functional operator acting on functions to give derivative functions. Derivative is generated when D is applied to functions whose derivatives the system does not know. Mathilda attempts to convert Derivative\[n1,...,nm\]\[f\] to a pure function. When f is a symbol carrying DownValues, the evaluator rewrites the head as Function\[{t1,...,tm}, f\[t1,...,tm\]\] with the rule expanded into the body, then differentiates that pure function. If no DownValue matches, the original Derivative form is returned. Attributes: Protected.

</details>

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (8)

```mathematica
In[1]:= Derivative[2][f][x]
Out[1]= Derivative[2][f][x]

In[2]:= D[%, x]
Out[2]= 0

In[3]:= D[f[x, y], y]
Out[3]= Derivative[0, 1][f][x, y]

In[4]:= f[x_] := x^5 + 6 x^3

In[5]:= f'[x]
Out[5]= 18 x^2 + 5 x^4

In[6]:= f'[5]
Out[6]= 3575

In[7]:= g[x_, y_] := x^2 y^3

In[8]:= Derivative[1, 1][g][a, b]
Out[8]= 6 a b^2
```

### Applications (7)

```mathematica
In[9]:= f'[x]
Out[9]= Derivative[1][f][x]

In[10]:= D[f[x], x]
Out[10]= Derivative[1][f][x]

In[11]:= Derivative[2][Cos]
Out[11]= Derivative[2][Cos]

In[12]:= D[f[g[x]], x]
Out[12]= Derivative[1][g][x] Derivative[1][f][g[x]]

In[13]:= D[f[g[x], h[x]], x]
Out[13]= Derivative[1][g][x] Derivative[1, 0][f][g[x], h[x]] + Derivative[1][h][x] Derivative[0, 1][f][g[x], h[x]]

In[14]:= D[(f[x])^2, x]
Out[14]= 2 f[x] Derivative[1][f][x]

In[15]:= Derivative[1][#^2 + 1 &]
Out[15]= 2 #1 &
```

## Algorithm

deriv.c -- Native C implementation of Mathematica-style differentiation.

This module replaces the fragile rule-based bootstrap in src/internal/deriv.m with a direct, dispatch-driven implementation.

Overview -------- The two key entry points are the builtins D (partial derivative) and Dt (total derivative). Both ultimately funnel through a single recursive core, ``compute_deriv``, parameterised by an optional differentiation variable. When the variable is non-NULL we compute a partial derivative treating everything else as constant (using a fast FreeQ-style walk to short-circuit constant sub-trees). When the variable is NULL we compute a total derivative -- unknown symbols then participate as ``Dt[sym]`` terms.

Why this is faster than the rule-based implementation ----------------------------------------------------- The old deriv.m relied on ~60 DownValues. Each call to D[f, x] would:

```text
  * scan the DownValues list for D linearly,
  * attempt pattern matching against every rule head (Plus, Times,
    Power, every elementary function, ...),
  * run ``/;`` side-conditions such as FreeQ,
  * perform attempt-evaluate/backtrack cycles in the matcher,
  * recursively re-evaluate the result through the full rule engine.
```

In contrast, this module performs a single head-symbol strcmp dispatch per call, constructs the derivative expression directly, and lets the outer evaluator simplify arithmetic. Crucially, the constant-detection step uses a tailored structural traversal (expr_free_of) that avoids calling out to the generic FreeQ builtin.

Returned expressions -------------------- Every builder below produces plain un-reduced expression trees (e.g. Plus[0, x] or Times[1, x]). The outer Mathilda evaluator runs a full fixed-point reduction on the value we return, so Plus[0, ...], Times[1, ...], and all subsequent chain-rule simplifications fold automatically. This keeps the code readable and avoids duplicating the arithmetic simplifier.

Memory ownership ---------------- Every helper that returns an ``Expr*`` returns a freshly allocated tree owned by the caller. Input expressions are never mutated; sub-expressions that need to be reused are always deep-copied.

## Implementation notes

**Algorithm.** `Derivative[n1, ..., nm]` is primarily a *tag* head: the actual
differentiation work is done inside the `D` dispatch (src/calculus/deriv.c).
`builtin_derivative` itself returns NULL, leaving `Derivative[n]` in canonical
unevaluated form; the builtin exists chiefly so attributes can be registered on
the symbol. Two pieces of real logic apply the tag:

- `derivative_of_pure_function(deriv_head, pure_fn)` differentiates
  `Derivative[n1,...,nm][Function[{t1,...,tm}, body]]` by partial-differentiating
  the body `ni` times in each slot `ti` via the shared `compute_deriv` core.
- `derivative_of_symbol(deriv_head, fsym)` reduces `Derivative[...][f]` when the
  symbol `f` carries DownValues: it mints fresh temporary slot symbols, builds
  and evaluates `f[t1,...,tm]` (triggering the DownValue rewrite), wraps the
  substituted body in a synthetic `Function`, and delegates to
  `derivative_of_pure_function`. If the call did not rewrite (no matching
  DownValue) it aborts to NULL. All `ni` must be nonnegative integers.

For unknown functions, `compute_deriv`'s chain rule emits `Derivative[...][f]`
factors, so the tag composes naturally through the rest of differentiation.

**Data structures.** `Expr*` trees; uses a static counter to generate
collision-free temporary slot-variable names (`Derivative$<id>$<k>`) without
registering them in the symbol table.

**Complexity / limits.** Only nonnegative integer derivative orders are
reduced; symbolic or negative orders stay unevaluated.

- `Protected`.
- Acts primarily as a tag carried through the differentiation
  pipeline: `D` and `Dt` produce `Derivative[...]` heads for
  unknown functions and advance their indices when differentiating
  an expression that already contains one.
- **Operator-form reduction**: `Derivative[n1, ..., nm][f]` is reduced
  at evaluation time when `f` has DownValues that match. The evaluator
  synthesises `Function[{t1, ..., tm}, f[t1, ..., tm]]` (with the
  DownValue rule expanded into the body) and differentiates that pure
  function via the existing pure-function pipeline. This makes
  `f'[x]` (i.e. `Derivative[1][f][x]`) compute the derivative of a
  user-defined `f`. When `f` has no matching DownValue the form is
  left unevaluated, matching Mathematica.

**Attributes:** `Protected`.

## References

**See also:** [D](../../calculus/D/), [Dt](../../calculus/Dt/)

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (Kluwer, 1992), ch. 2.
- Source: [`src/calculus/deriv.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/deriv.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_beta.c`](https://github.com/stblake/mathilda/blob/main/tests/test_beta.c)
- Tests: [`tests/test_collect_corpus.c`](https://github.com/stblake/mathilda/blob/main/tests/test_collect_corpus.c)
- Tests: [`tests/test_deriv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_deriv.c)
- Tests: [`tests/test_findroot.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findroot.c)

## Notes & additional examples

### Notes

`Derivative[n][f]` is the functional operator representing `f` differentiated `n` times; the surface forms `f'` and `f''` parse to `Derivative[1][f]` and `Derivative[2][f]`. It is the object `D` generates whenever it differentiates an unknown function head, which is why `D[f[x], x]` returns `Derivative[1][f][x]` and the chain rule on `f[g[x]]` yields a product of `Derivative[1]` operators. Note that `Derivative` does not auto-resolve against the known elementary table here: `Derivative[1][Sin]` and `Derivative[2][Cos]` stay in operator form rather than collapsing to `Cos` or `-Cos`. Apply the operator to an explicit argument (e.g. `Derivative[1][f][a]`) to obtain the evaluated-at form.
