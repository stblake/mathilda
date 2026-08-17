# Task: Normalise parser precedences onto a clean 1–10000 ladder

## Plan
- [x] Map the full precedence surface (parse.c operator table + inline literals;
      print.c `get_expr_prec` + inline `parent_prec` args; confirm print_latex.c
      is an independent scale to leave alone).
- [x] Choose scheme (user: clean ladder in [1,10000], Part capped at 10000).
- [x] Capture a behavioural baseline (FullForm/InputForm/TeXForm over an
      operator-dense corpus) from the pre-change binary.
- [x] Apply the order-preserving remap to src/parse.c and src/print.c.
- [x] Fix stale in-code comments citing old numbers.
- [x] Build + `make check-c99`.
- [x] Differential: rerun corpus, diff vs baseline (must be empty).
- [x] Run parse_tests / parser_precision_tests / print_tests + eval smoke test.
- [x] Update docs: docs/spec/operators.md, SPEC.md §3.2, changelog.

## The ladder (OLD → NEW)
`;`10→100 · Put 30→300 · Set-family 40→500 · Postfix 70→800 · `&`90→1000 ·
ReplaceAll 110→1200 · Rule-family 120→1500 · Condition 130→1700 · Optional 140→1900 ·
StringExpr 155→2100 · Alternatives 160→2300 · Repeated 170→2500 · And/Or 215→2800 ·
Not 230→3000 · comparisons/Span 290→3200 · Plus 310→3500 · Times 400→4500 ·
Divide/Rational 470→5000 · unary-minus 480→5100 · Dot 490→5300 · Power 590→6500 ·
(Power+1) 591→6501 · StringJoin 600→6700 · Prefix/Apply/Map 620→7000 · Composition 625→7200 ·
Incr/Decr 660→7500 · Derivative 670→7700 · PatternTest 680→7900 · Factorial 710→8200 ·
MessageName 780→8600 · Call/atom-default 1000→9500 · Part 1100→10000.
Plus computed forms: top-level min_prec 11→101, Span prec+1 291→3201.

## Review
**What changed:** A pure, order-preserving renumbering of every parser precedence
onto a 1–10000 ladder, applied in lockstep to the parser table (`get_operator`,
src/parse.c) and the printer's parenthesiser (`get_expr_prec` + ~30 inline
`parent_prec` args, src/print.c). `src/print_latex.c` uses a separate
self-contained `PREC_*` scale and was intentionally left untouched.

**Why it's safe:** the mapping is order-isomorphic (old_a<old_b ⟺ new_a<new_b,
equals→equals). The Pratt loop and the parenthesiser make only `<`/`==`
comparisons on these numbers, so all parse trees and printed forms are unchanged.

**Verification (all green):**
- Differential over a 73-case operator-dense corpus (FullForm + InputForm +
  TeXForm): baseline vs after **diff is empty**.
- `make check-c99` clean; full `make` links.
- parse_tests, parser_precision_tests, print_tests all pass.
- End-to-end eval smoke test (arithmetic assoc, Power, Part, Solve, Map, D,
  Expand, Rational) correct.

**Result:** adjacent precedence levels now sit ≥200 apart (vs as little as 5
before), leaving room to insert new operators without renumbering.

**Optional follow-up (not done):** replace the magic numbers with shared
`#define PREC_*` constants so future operator insertions can't silently drift the
parser and printer apart.
