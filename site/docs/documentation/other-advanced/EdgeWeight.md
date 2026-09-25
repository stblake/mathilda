# EdgeWeight

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeWeight[g] gives the weights of g's edges, in EdgeList order. Defaults to all 1s if g was built without an EdgeWeight option.`**

## Examples

_No verified examples yet for this function._

## Algorithm

edgeweight.c - EdgeWeight[g]: the graph's per-edge weights, in EdgeList order. Defaults to List[1, 1, ..., 1] (one per edge) when g carries no EdgeWeight -- matching Wolfram Language's own behavior for an unweighted graph, and giving WeightedAdjacencyMatrix[g] a well-defined answer for every valid graph, not just ones explicitly built with weights.

Memory (SPEC section 4): returns a fresh list; the evaluator frees res.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
