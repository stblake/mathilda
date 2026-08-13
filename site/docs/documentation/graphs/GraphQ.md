# GraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphQ[g] gives True if g is a valid graph, and False otherwise.`**

## Examples

_No verified examples yet for this function._

## Algorithm

graphq.c - GraphQ[g]: is g a valid graph?

A thin wrapper over graph_is_valid (graph_util.c): returns the symbol True when the (already-evaluated) argument is a canonical, valid graph, and False otherwise. Non-unary calls are left unevaluated (NULL).

Memory (SPEC section 4): returns a freshly-allocated symbol; the evaluator frees `res`.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
