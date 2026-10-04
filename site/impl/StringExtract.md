---
source: src/strings/stringextract.c
---
**Algorithm.** `builtin_stringextract` splits a string into blocks and selects by position, treating each argument after the subject as one level (fewer than two arguments emits `StringExtract::argm`). The split at every level is delegated to `StringSplit` (synthesised as `StringSplit[str, sep]` and evaluated), so the whole string-pattern engine, whitespace-run collapsing, and empty-end trimming are reused verbatim — making `StringExtract[s, sep -> All]` exactly `StringSplit[s, sep]`. `apply_position` resolves `n`/`-n` (`pos_single`), `All` (every block), `{n…}` (a collection, deferring any `Span`/`All` element to `Part` semantics), and `m ;; n` (a `Span`, via `expr_part`). An explicit `sep -> pos` rule carries its own separator; a bare position gets a depth-default separator — whitespace at the lowest level, growing runs of `"\n"` above. Multi-level specs recurse into each selected block.

**Data structures.** A `LevelSpec` table of borrowed-or-owned separator/position pairs; results nest as `List`s.

**Complexity / limits.** Dominated by the delegated `StringSplit` evaluations. An out-of-range single index yields `Missing["PartAbsent", n]` rather than leaving the call unevaluated; a non-string subject (or list element) does leave it unevaluated. A list subject threads.
