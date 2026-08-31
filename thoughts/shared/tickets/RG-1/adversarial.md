# Adversarial review — RG-1 (`RandomGraph[{n,m},k]`)

**Date:** 2026-08-30 22:54:35 (written); HIGH re-verified and resolved 2026-08-30 23:03:49
**Scope:** uncommitted working-tree diff vs `7701df5b` — `src/graph/generators.c`,
`src/graph/graph.c`, `src/graph/graph.h`, `tests/test_graph.c`,
`docs/spec/builtins/graphs.md`, `docs/spec/changelog/2026-08-24.md`
**Plan:** `thoughts/shared/tickets/RG-1/plan.md`

Ticket scope (AC-1 … AC-17, AC-19) is implemented and matches the docs, docstring, and
changelog. Ownership on the new paths is correct: `vertex_list` frees the memcpy'd array,
`cand` is freed after `expr_new_function`, `sampled` is freed on the reject path, the
`k`-loop frees accumulated graphs and `gs` on partial failure, `res` is never freed, and
`NULL` is returned for every decline.

## ~~HIGH~~ → RESOLVED — the `maxe` overflow guard still wraps; out-of-bounds heap write

> **RESOLVED 2026-08-30 23:03:49** (`src/graph/generators.c` mtime; working tree on top of
> HEAD `7701df5b` — the fix is uncommitted, so it has no SHA of its own yet).
> Re-verified by reading the current file: `src/graph/generators.c:174` now reads
> `if (n > 2147483647L) return NULL;` and sits **before** the multiply at `:175-176`, with
> the comment block at `:166-173` naming this exact `n = 2^32 + 1` witness as its reason.
> With `n < 2^31` the product `n(n-1)` cannot exceed `2^62`, so the `unsigned long long`
> arithmetic is exact and the downstream `SIZE_MAX` bound is meaningful. The wrap described
> below is unreachable, and the 17 GB `calloc` is never attempted.
> Live check (now safe, since both witnesses decline before allocating):
> `RandomGraph[{4294967297, 1}]` and `RandomGraph[{2147483648, 1}]` both return
> unevaluated; `RandomGraph[{5,4}]` still returns `Graph[<5 vertices, 4 edges>]`.
> The finding below is retained as written for the record — it describes the pre-fix file.

`src/graph/generators.c:169-171`, exploited at `:136-142`.

`maxe = (unsigned long long)n * (n-1) / 2` is computed in the same 64-bit width it is
meant to protect, so it wraps silently. Witness `n = 2^32 + 1 = 4294967297`, `m = 1`:

- `n*(n-1) = 2^64 + 2^32` → wraps to `2^32`, so `maxe = 2^31`.
- `2^31 < SIZE_MAX / sizeof(Expr*)` (≈2^61), so the bound at `:171` passes.
- `m = 1 <= maxe`, so `:172` passes.
- `:137` allocates `calloc(2^31, 8)` = 17.2 GB — under Linux overcommit this routinely
  succeeds.
- `:140-142` then iterates over the **true** `n(n-1)/2 ≈ 9.2e18` candidate pairs, writing
  `cand[k++]` past the `2^31`-th slot. Heap corruption, not a `NULL`-checkable failure.

The loop bound is `n`, not `ncand`; nothing reconciles the two. Where the `calloc` happens
to fail the result is a benign `NULL`, so the outcome is allocator- and platform-dependent,
which is worse than a deterministic decline.

AC-18's chosen witness (`n = 2^62`) wraps to a value *above* the `SIZE_MAX` bound and so
is rejected by coincidence — the test passes while the bug class is untouched.

Fix: bound `n` itself before multiplying (e.g. reject `n` above ~`2^31`, well past any
allocatable candidate set), or detect the wrap in the product, rather than checking the
product after the fact.

**Verification status:** CONFIRMED by code reading (arithmetic and the loop bound both
checked in the file). Deliberately **not** exercised on the live binary — running the
witness risks a 17 GB allocation followed by an out-of-bounds write.

## MEDIUM — `int_vertices` returns `NULL` on allocation failure and is dereferenced

`src/graph/generators.c:46-47`, newly reached via `vertex_list` at `:108`.

The `calloc` result is unchecked and `v[i] = ...` dereferences it. Pre-existing code, but
the diff adds a caller and the changelog advertises that "both `calloc`s are
`NULL`-checked" — there are three allocation sites on this path and this is the unchecked
one.

Still present as of 23:03:49 (`:45-49` in the current file — `Expr** v = (n > 0) ? calloc(...)`
followed by an unguarded `v[i] = ...`). The reachability argument has changed with the HIGH's
fix: the multi-billion-vertex route via a wrapped `maxe` is now closed, so this is reachable
only through a genuine allocation failure at a large-but-legal `n` (up to `2^31 - 1`
vertices), not through the overflow. Severity stays MEDIUM.

## LOW — the changelog's leak claim has no artifact behind it

`docs/spec/changelog/2026-08-24.md`: "`Do[RandomGraph[{1, 0}], {50}]` is now completely
leak-free". No valgrind output or scripted leak check is in the diff, so if the claim
regresses nothing in the tree notices. `tests/scripts/geometry_leakcheck.sh` — cited in the
same changelog file — is the existing pattern for exactly this.

## Could not assess

- Whether the C suite passes, and whether AC-18's witness returns unevaluated on this
  build — the 466-binary suite was explicitly out of scope for this run.
- AC-17 (RNG stream position unchanged) — needs both binaries built and a seeded draw
  diffed. Inspection shows no RNG consumption added or removed on the `n >= 2` path, which
  is consistent with the claim but is not the measurement AC-17 asks for.
- No domain checklist configured (`.claude/GUIDANCE_ROLES.md` absent), so this ran on the
  generic list alone.
