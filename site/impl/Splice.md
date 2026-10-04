---
source: src/assoc_ops.c
---
**Algorithm.** `Splice` is an evaluator hook, not an ordinary builtin. During the
evaluator's Sequence-flattening step (which runs for heads without `SequenceHold` /
`HoldAllComplete`, and only once a `Splice` argument has actually been seen, so no
other call pays for it), `eval_splice_args` replaces each qualifying `Splice[{e…}]`
with its elements. `splice_applies` gates this: the one-argument `Splice[list]`
splices only inside a `List` or an `Association`; the two-argument `Splice[list, h]`
splices into any head matching the pattern `h` (checked by `MatchQ`), but never into
an association. Anywhere else `Splice` stays inert.

**Data structures.** A `hit` bitmap marks the arguments to expand; the enclosing
node's argument array is rebuilt in place. A packed inner list is materialised to a
nested `List` first.

**Complexity / limits.** O(total arguments after expansion). Because the hook keys on
a `Splice` head actually appearing among the arguments, calls with no `Splice` are
untouched. A bare `Splice[{…}]` at top level (not inside a List/Association, no head
spec) is left unevaluated.
