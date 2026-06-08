# Scoping Constructs

10 built-in function(s) in this category.

- [`$Context`]($Context.md) — $Context is a string giving the current context. New symbols are created  _(Stable)_
- [`$ContextPath`]($ContextPath.md) — $ContextPath is a list of contexts used (in order) to resolve bare  _(Stable)_
- [`Begin`](Begin.md) — Begin["ctx`"] sets the current context ($Context) to "ctx`", saving  _(Stable)_
- [`BeginPackage`](BeginPackage.md) — BeginPackage["ctx`"] sets the current context to "ctx`" and restricts  _(Stable)_
- [`Block`](Block.md) — Block[{x, y, ...}, expr] evaluates expr with local values for x, y, ....  _(Stable)_
- [`Context`](Context.md) — Context[] gives the current context ($Context).  _(Stable)_
- [`End`](End.md) — End[] restores the context that was active before the matching Begin[]  _(Stable)_
- [`EndPackage`](EndPackage.md) — EndPackage[] restores the state saved by BeginPackage and prepends the  _(Stable)_
- [`Module`](Module.md) — Module[{x, y, ...}, expr] specifies that x, y, ... are local variables.  _(Stable)_
- [`With`](With.md) — With[{x = x0, ...}, expr] specifies that x should be replaced by x0 throughout expr.  _(Stable)_
