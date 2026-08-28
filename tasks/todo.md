# Todo: Book additions + Reduce list-form fix

Plan: `/Users/user/.claude/plans/immutable-knitting-puddle.md`

## Task D — Fix `Reduce[{eqns...}, vars]` (list = conjunction)  [do first]
- [ ] `src/solve/reduce.c`: List→And conversion in `builtin_reduce`, free `owned_list` on all return paths
- [ ] Docstring line for the list form
- [ ] `docs/spec/builtins/solutions-of-equations.md` note
- [ ] `tests/reduce_corpus.m` case
- [ ] Rebuild; verify list form + 1-elem list; valgrind

## Task A — CRT example (§4.6.5)
- [ ] Rewrite `book/examples/number-theory/crt.m` (lead with ChineseRemainder)
- [ ] Rewrite prose in `book/chapters/math/number-theory.tex`; add `\usagebox{ChineseRemainder}`

## Task B — New §3.6.2 "Systems of equations"
- [ ] New `book/examples/03-introduction/solve-systems.m`
- [ ] New subsection in `book/chapters/03-introduction.tex`

## Task C — EliminationOrder explicit
- [ ] §3.5 prose (`03-introduction.tex`, pair algebra-groebner/2)
- [ ] §4.2.7 clause (`chapters/math/algebra.tex`)

## Finalize
- [ ] `docs/spec/changelog/2026-08-24.md` entry
- [ ] `cd book && make examples && make usage && make pdf`
- [ ] Spot-check rendered PDF

## Review — all tasks complete & verified

**Task D — `Reduce` list form (bug fix).** `src/solve/reduce.c`: a `List` in the
statement slot is rewritten to `And` (mirrors `reduce_int.c`), owned temp freed
on all 6 return paths. Valgrind: list-form and `&&`-form leak totals identical
(13,440/420) → zero new leaks. Docstring + `solutions-of-equations.md` note.
`reduce_corpus.m` +2 cases; verifier `reduce_check_prelude.m` normalizes List→And
for back-sampling. Tests: reduce_tests + reduce_corpus (160/160) + 5 others all
pass; `make check-c99` clean.

**Task A — CRT example (§4.6.5).** Rewrote `crt.m` + prose: leads with the
builtin, keeps Sun Zi by-hand recipe as "under the hood", adds offset form (128)
and inconsistent decline; `\usagebox{ChineseRemainder}` card renders (PDF p.98).

**Task B — new §3.6.2 "Systems of equations".** `solve-systems.m` + subsection:
circle∩line attacked by LinearSolve/RowReduce/Solve/GroebnerBasis/Resultant/
Reduce(list form)/ChineseRemainder, threaded by elimination. Renders PDF p.24-25;
cross-refs (§4.4, §4.6, §4.2.9-10) resolve; Numerical/Continuum renumbered 3.6.3/4.

**Task C — EliminationOrder explicit.** §3.5 caption + §4.2.7 clause name
`MonomialOrder -> EliminationOrder`. (Finding: the 3-arg form already forces it;
a bare `EliminationOrder` without naming the elim var does NOT eliminate — so the
example was correct, only the prose needed it.)

**Finalize.** Changelog `2026-08-24.md` (2 new sections, newest-first). Book PDF
rebuilt (128 pp), index regenerated with new entries (Chinese Remainder Theorem,
congruence!system of, monomial order!elimination, system of equations!linear/
polynomial, Sun Zi). No LaTeX errors.

**Flagged, out of scope:** `ChineseRemainder` absent from
`site/docs/assets/builtins.json`/site docs (so `\B{ChineseRemainder}` has no
hyperlink) — separate site-generation subsystem; left for a follow-up.
