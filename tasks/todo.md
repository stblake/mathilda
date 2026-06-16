# Task: Trivial + highly-non-trivial examples for all non-trivial builtins

Goal: For every non-trivial builtin, the doc page should show BOTH a simple/trivial
example AND one or more highly non-trivial ("research level") examples that show off
Mathilda's capabilities. Examples live in `site/overlays/<Name>.md` (the established
"Worked examples" channel) and MUST be verified against the built `./Mathilda`.

## Decisions (confirmed with user)
- Scope: all non-trivial functions — 320 candidates (448 total minus 128 trivial
  accessors / `$`-sysvars / simple File I/O / RNG / basic comparisons & assignment).
- Placement: `site/overlays/<Name>.md` (NOT spec files). Overlays render under
  "Notes & additional examples". Not auto-verified by generate.py → each example
  verified by hand against the binary during authoring.

## Pipeline
- `site/generate.py` merges `overlays/<Name>.md` onto each page. Front matter
  (`status:`, `references:`) overrides; body appended. Examples are plain
  ```mathematica In[n]:= / Out[n]= ``` blocks.
- After all overlays written: `python3 site/generate.py` (needs `./Mathilda`),
  then `make docs-build` (strict MkDocs) to confirm the site still builds.

## Plan
- [x] Understand generator + overlay format
- [x] Confirm scope + placement with user
- [x] Compute 320 candidates -> 30 batches (/tmp/mathilda_doc_candidates.json)
- [ ] PILOT: one agent on batch 06; review quality + round-trip
- [ ] Fan out remaining 29 batches (subagents, ~11 fns each)
- [ ] Regenerate site (generate.py) + strict build
- [ ] Spot-check a sample of rendered pages for fabricated outputs
- [ ] Commit + push

## Constraints for authoring agents
- Run `./Mathilda` FOREGROUND only (printf ... | ./Mathilda). No backgrounding,
  no until-loop pollers (OOM risk), no `timeout` (macOS lacks it).
- NEVER fabricate Out[]; paste the exact verified output.
- Preserve existing overlay front matter (status/references) and Notes prose.
- Only edit `site/overlays/<Name>.md`. Do not touch generate.py / spec / src.

## Review
- Done: 320 candidate builtins processed across 30 batches (1 pilot + 29 fan-out
  subagents). Overlays merged grew 215 -> 396. `generate.py` re-verified examples;
  strict MkDocs build passes (1134 spec-verified + many overlay examples).
- Every overlay example was run through `./Mathilda` by the authoring agent; I
  additionally reproduced the wrong-result bug claims before logging them.
- Bugs/discrepancies consolidated in `MATHILDA_BUGS.md` (8 confirmed wrong-result/
  broken items incl. Coefficient-Laurent, EvenQ-of-symbol corrupting Sum, matrix
  Norm, Flatten[_,Infinity], Root[Function]-N, Variables[eqn]; plus a discrepancy
  list and unimplemented-function list).
- Known cosmetic nit (not fixed): a few multi-step overlay blocks reuse `In[1]`
  across lines instead of 1,2,3; outputs are all correct/verified.
