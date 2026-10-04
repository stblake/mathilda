# Trace

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Trace[expr]`**

Generates a nested list of the expressions produced while evaluating expr. Each argument sub-evaluation that takes a step appears as a sublist, mirroring the structure of the evaluation. Returns {} when expr needs no rewriting.

**`Trace[expr, form]`**

Includes only the steps whose expression matches the pattern form (e.g. Trace\[expr, \_Integer\]).

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= Trace[1 + 1]
Out[1]= {1 + 1, 2}

In[2]:= Trace[5]
Out[2]= {}

In[3]:= Trace[2^3 + 4^2 + 1]
Out[3]= {{2^3, 8}, {4^2, 16}, 8 + 16 + 1, 25}

In[4]:= Trace[x^Range[10]]
Out[4]= {{Range[10], {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}}, x^{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {x, x^2, x^3, x^4, x^5, x^6, x^7, x^8, x^9, x^10}}

In[5]:= g[y_] := y^2; Trace[g[1 + 1]]
Out[5]= {{1 + 1, 2}, g[2], 2^2, 4}

In[6]:= Trace[1 + 2 + 3, _Integer]
Out[6]= {6}

In[7]:= Trace[Nest[f, x, 3], _f]
Out[7]= {f[f[f[x]]]}
```

### Applications (5)

```mathematica
In[8]:= Trace[1 + 1]
Out[8]= {1 + 1, 2}
```

An inert atom needs no rewriting

```mathematica
In[9]:= Trace[5]
Out[9]= {}
```

Nested, mirroring the evaluation structure

```mathematica
In[10]:= Trace[2^3 + 4^2 + 1]
Out[10]= {{2^3, 8}, {4^2, 16}, 8 + 16 + 1, 25}
```

Filtered to steps matching a pattern

```mathematica
In[11]:= Trace[1 + 2 + 3, _Integer]
Out[11]= {6}
```

```mathematica
In[12]:= Trace[Nest[f, x, 3], _f]
Out[12]= {f[f[f[x]]]}
```

## Algorithm

Trace[expr] / Trace[expr, form] — user-facing builtin (beads-planning-b3k Phase 3; the two-arg form is beads-planning-h4u).

A thin wrapper over the evaluator-side collector eval_collect_trace() (src/eval.c). Trace[expr] returns a NESTED List that mirrors the structure of expr's evaluation: each argument sub-evaluation that takes a step appears as a sublist, and the reassembled intermediate form appears as a step (matching Mathematica). E.g. Trace[2^3 + 4^2 + 1] -> {{2^3,8},{4^2,16},8+16+1,25}.

```text
  - >=1 rewrite -> the nested list of forms
  - no rewrite (inert atom / normal form) -> {}
```

Trace[expr, form] filters that nested trace to the step leaves whose expression matches the pattern `form` (same structural matcher as MatchQ, src/match.c), flattening the nesting away into a plain List of matches. Because Trace is HoldAll, `form` reaches the builtin unevaluated, so pattern literals like `_Integer` or `f[_]` work directly. A form that matches nothing yields {}.

All the nesting/clock/ownership subtlety lives in eval_collect_trace; this file only handles arity dispatch, form filtering, HoldForm wrapping, registration, attributes, and the docstring. Arities other than 1 or 2 return NULL so the call stays unevaluated rather than silently misbehaving.

## Implementation notes

**Algorithm.** `builtin_trace` (`src/trace.c`) is a thin wrapper over the
evaluator-side collector `eval_collect_trace` (`src/eval.c`), which re-runs the
argument's evaluation while recording every form it passes through. Because
`Trace` carries `HoldAll`, the argument reaches the builtin unevaluated, so its
rewrite sequence is observed from the start (and the evaluation clock is bumped
once so an already-evaluated argument is still traced in full). The collector
returns a **nested** `List` mirroring the evaluator's own recursion: each argument
or head sub-evaluation that takes a step becomes a sublist, one that takes no step
contributes nothing, and the reassembled intermediate form appears as a step.
`Trace[expr, form]` then filters that tree with `trace_collect_matches`, keeping
only step leaves that structurally `match` the held pattern `form` and flattening
the nesting into a plain `List`.

**Data structures.** The raw collector output is an `Expr` tree whose `List` nodes
are the nesting markers. `trace_holdform_tree` walks it and wraps every non-`List`
leaf in `HoldForm`, so the returned structure is inert under the evaluator's
fixed-point re-pass (a bare `1 + 1` would otherwise reduce again to give `{2, 2}`)
while still printing transparently. The two-argument filter accumulates matches in
a growable `Expr**` buffer, copying each matched node.

**Complexity / limits.** Proportional to the number of evaluation steps `expr`
takes. A builtin's internal computation and `Listable` threading show as a single
atomic rewrite (matching Mathematica — `Range[10]` is one step). `Trace` is
reentrant (an inner `Trace` appears as one reduced value to the outer). Arities
other than 1 or 2 return `NULL` (stay unevaluated); `TraceDepth` is not
implemented. Attributes `HoldAll`, `Protected`.

- `HoldAll`, `Protected`. The argument is held so its rewrite sequence can be
  observed from the start.
- **Nested.** Evaluation semantics are unchanged — the collector merely observes
  the evaluator's own recursion. Argument and head sub-evaluations that take a
  step become sublists; a sub-evaluation that takes no step contributes nothing.
- **Automatic-rewrite steps are atomic.** A builtin's internal computation (e.g.
  `Range` folding) and `Listable` threading of the threaded elements are shown as
  a single rewrite, not decomposed — so `Range[10]` and `x^{1..10}` each appear
  as one step, matching Mathematica.
- Each step leaf is returned wrapped in `HoldForm` (printed transparently), so
  the result stays inert and does not re-evaluate.
- Not memoized: `Trace` bumps the evaluation clock once per call so an
  already-evaluated argument is still traced in full.
- Reentrant: `Trace` may appear inside a traced expression; the inner list is
  produced independently and appears to the outer trace as one reduced value.
- `TraceDepth` and arities other than 1 or 2 are not implemented; e.g.
  `Trace[expr, form, extra]` stays unevaluated.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [MatchQ](../../pattern-matching/MatchQ/), [HoldAll](../../expression-information/HoldAll/), [Range](../../lists-and-iteration/Range/), [HoldForm](../../expression-information/HoldForm/)

- Source: [`src/trace.c`](https://github.com/stblake/mathilda/blob/main/src/trace.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_print.c`](https://github.com/stblake/mathilda/blob/main/tests/test_print.c)
- Tests: [`tests/test_trace.c`](https://github.com/stblake/mathilda/blob/main/tests/test_trace.c)

## Notes & additional examples

### Notes

`Trace[expr]` returns a **nested** list of the forms `expr` passes through while
evaluating, mirroring the evaluator's own recursion: each argument sub-evaluation
that takes a step becomes a sublist, and the reassembled intermediate form appears
as a step. An expression that needs no rewriting (an inert atom, a normal form)
traces to `{}`. A builtin's internal work and `Listable` threading are shown as a
single atomic step, matching Mathematica — `Range[10]` is one step, not ten.

`Trace[expr, form]` keeps only the step leaves whose expression matches the
pattern `form`, flattening the nesting into a plain list. `form` is held, so
pattern literals such as `_Integer` or `f[_]` may be written directly. Each step
is returned wrapped in `HoldForm`, so the result prints transparently yet stays
inert and does not re-evaluate.
