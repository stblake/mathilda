# Session 3 Rehearsal: `docs/TEACHING.md`'s Own Weakness List, Checked

**Kit version checked against: `9.0.7`, commit `288bdd55f7076bfc9a993225e80affcb88a8e738`**
(same clone as `teaching-session1-rehearsal.md`, in this same repo). Read-only against the
kit repo throughout, as instructed.

**First correction to the request itself, stated up front rather than buried**: Session 3
currently lists **five** weak items, not four. A fifth — idea-stage fabrication
(`feature-discovery` → `idea-to-prd` → `prd-to-epic`) — is present in the live text and is
not in the four summarized in the request. Checked below alongside the other four.

---

## Item 1 — handoff citations not verified

**Session 3's current text**: states the weak claim as "nothing mechanically checks a
citation against the file it names before writing it," tells the reader to "confirm this
yourself by reading `commands/handoff.md` at the current commit," and separately names the
provenance-tagging fix as the thing to demo instead.

**Checked**: `grep` for any citation-verification logic in `commands/handoff.md` — zero
hits. The weak claim is still true, verified directly, not assumed from the doc's own
wording.

**On the rewording specifically, since this is the one the request flagged as most
important**: **it landed, and it landed correctly.** The live Session 3 text does **not**
assert the "nine of nine wrong" anecdote to the room at all — it never appears in
`docs/TEACHING.md`'s own prose. It points instead at `OPEN-QUESTIONS.md` item G-10 as
supporting material and, separately, at the directly-checkable claim (no verification step
exists). A room following the text as written will not be handed an anecdote to repeat.

**But the anecdote itself is worth knowing about, one level down, because another
session's own audit already found the exact problem the request was worried about — and I
verified that audit's own claim rather than just relaying it.** `OPEN-QUESTIONS.md` (around
line 1167) contains a self-audit, dated the same night as this one: it states plainly that
the two incidents the kit's own artifacts cite as motivating cases — this "nine of nine"
handoff and a research doc's "confirmed directly with the maintainer" line — are **not
independently locatable** in the repo's current state, and that `ADR-0006` "names no file
path for the case it describes." I re-verified this claim myself rather than take it as
settled: `grep -rn "9 of 9"` across every tracked `.md` file in the repo returns exactly the
one line inside that same self-audit paragraph — the incident is not cited anywhere else as
a standalone artifact.

**One thing worth telling Michael that this self-audit gets slightly wrong, or at least
states more broadly than the evidence supports.** The audit treats both cited incidents as
equally unlocatable. They are not the same case. The "confirmed directly with the
maintainer" incident is **this session's own GR-12 finding**, from earlier tonight — a real
event, verified directly against a real `AskUserQuestion` exchange in this session's own
transcript, that led to a real, shipped fix (the four-state provenance vocabulary, confirmed
present in all three named files below). It is not written to disk as a standalone
`thoughts/shared/research/` artifact, which is exactly what the self-audit checked for and
correctly found absent — but "not committed as a standalone file" and "not independently
verifiable" are different claims, and this one is verifiable, just not by grepping the repo.
The "nine of nine" handoff incident may be the same situation (a real test case run live,
never saved as a standalone artifact) or it may genuinely be unlocatable even in that
sense — I have no way to tell which from here, and neither does the self-audit, which is
exactly why it hedges rather than asserts. Worth saying to whoever owns this: "not in the
repo" is not automatically "did not happen," and conflating the two in the audit's own
prose is a smaller version of the identical problem the audit exists to catch.

**Provenance vocabulary landing, checked directly rather than trusted from that same
self-audit**: `grep -l` for the four terms (`stated-by-human`, `chosen-from-options`,
`accepted-default`, `model-inferred`) across `commands/handoff.md`, `commands/create-plan.md`,
and `skills/grill-me/SKILL.md` — all three match. Confirmed independently.

**Verdict: still true, correctly reworded, demo substitute genuine.**

---

## Item 2 — `verify-implementation`'s diff range can be empty at the moment it's first run

**Checked directly against `commands/verify-implementation.md`**: the diff-range logic
(`git merge-base HEAD origin/main`, `HEAD~1` fallback, `EMPTY DIFF` reported honestly rather
than a false pass) is unchanged from what `OPEN-QUESTIONS.md`'s G-13/G-14 entry describes.
That entry is itself still open (no "fixed" marker, ends with "Fix needs either..."),
including a sharper live-reproduced sub-case beyond the original claim: on a branch already
caught up with `origin/main`, `merge-base` resolves to `HEAD` itself, so even the documented
`HEAD~1` fallback never fires, and the tool still doesn't name "already caught up" as its own
empty-range case.

**Demo substitute** ("run it after the first real commit exists... show the ladder's
per-phase output there"): genuine and matches this session's own experience exactly — every
`/verify-implementation`-style check run tonight (both graph tickets) was run after a real
commit existed, and the per-phase ladder output was the actual load-bearing signal both
times, diff-range timing never came up as an issue in practice.

**Verdict: still true, substitute genuine and independently corroborated by this session's
own unrelated work tonight.**

---

## Item 3 — telemetry measures kit adoption, not all agentic work

**Checked**: `docs/AGENTIC_LADDER.md:62` states this limit in almost these exact words — "The
scope limit: this measures adoption within the toolkit, and nothing beyond it" — and line 42
repeats it for one specific excluded-events count. `commands/telemetry-report.md` exists as
the named demo substitute.

**Verdict: still true, stated plainly in the kit's own material, demo substitute real.**

---

## Item 4 — `docs-site`'s dangling-reference detector, known bug, pinned not fixed

**Checked**: `OPEN-QUESTIONS.md`'s G-28 entry confirms this is still open, still
deliberately unfixed (the reasoning given — rewording the source to dodge a detector bug
would hide the bug rather than leave it visible — matches Session 3's own framing exactly).

**One thing Session 3's phrasing undersells, worth flagging rather than letting pass**: the
live G-28 entry is headed *"Widened 2026-08-22, still not fixed... Now 5 pinned false
positives, not 2."* The bug's surface has grown since it was first named, not stayed flat.
Session 3's text ("a good example of the kit choosing a visible known issue over a quiet
workaround") is accurate about the *policy* being sound, but doesn't convey that the
*count* is trending up, not just sitting pinned. A room that hears "known, pinned, stable"
and later learns it went from 2 to 5 in one day will read that as something being
downplayed, even if the underlying decision (don't dodge the detector) was the right one
both times.

**Verdict: still true; the growth trend is a real, small understatement worth correcting in
the material.**

---

## Item 5 — idea-stage fabrication (not in the four the request named)

**Checked**: `OPEN-QUESTIONS.md` around line 1077 has an entry headed *"OPEN (2026-08-23,
cold end-to-end run of the idea stage) — a fourth failure axis: fabricated"* — dated the same
night as everything else tonight, describing a fresh, independently-run cold test that
produced exactly what Session 3's text now describes: invented Success Criteria numbers (a
90% figure, a 100% bar) with no traceable source, and a full tech-stack assertion against a
directory confirmed empty. Marked `OPEN`, i.e. not yet fixed as of this commit.

**Verdict: real, current, freshly added — and the single most severe item on the list by
this session's own judgment, see below.**

---

## The judgment asked for directly, not relayed from anyone else

**Is it the right list?** Four of the five items are accurate, verifiable, and each names a
real demo substitute that actually works — that's a genuinely well-maintained list, checked
line by line rather than trusted. The fifth (idea-stage fabrication) is the right item to
lead with if this session had to rank them: it is the only one where the kit doesn't just
have an unverified mechanism (item 1) or a timing edge case (item 2) or a scope boundary
(item 3) or a cosmetically-noisy-but-harmless detector bug (item 4) — it is the only one
where the kit **actively invents specific, false numbers and presents them with full
confidence**, in the artifact type (a PRD/epic) most likely to reach a client's hands
unreviewed. Session 3's own text already treats it this way ("the one most worth naming
before anyone runs it live") — that ranking is correct, independently confirmed.

**What I would add, from my own night's findings, that this list does not mention at all:**

1. **The kit's own commands are not invocable by an autonomous agent, at all, mid-session,
   after install — not just the `/reload-plugins` footnote item 1 alludes to, the whole
   thing.** This session installed the plugin successfully, confirmed it enabled at the
   correct version, and then could not dispatch a single `ais`-provided command or skill for
   the rest of the session — not `/research-codebase`, not `/create-plan`, not
   `/implement-plan`, not `/verify-implementation`, nothing. Every one of those was instead
   read from the plugin's own cache directory as a markdown file and followed by hand. This
   is now documented in the kit's own `README.md` (found during a separate re-check
   tonight) with the exact right instruction for an agent that hits it — "there is no retry
   that fixes it from inside the session... escalate to a human at an actual terminal" — but
   it is not in Session 3's own list of named weaknesses, and it is arguably a bigger one
   than any single item currently there: it means the entire RPI+plan-review loop this
   session ran twice, tonight, on a real codebase, only works with a human physically
   re-invoking each step at an actual keyboard. A room evaluating this kit for "automatic
   RPI" (the goal named in this project's own orchestrator handoff document) deserves to
   hear that plainly, in the same room where the other five limits get named.

2. **A confidently-wrong proposal, found this afternoon, in the exact demo this list is
   meant to set up for.** `docs/exemplars/teaching-session1-rehearsal.md` (this repo,
   committed earlier today) found that Session 1's own opening demo — DETECT reading a bare
   Makefile's real targets back to the room — does not fire on this specific repo, and
   instead proposes CMake/ctest/clang-tidy commands that do not build or test it at all,
   sourced from a real but structurally-misleading manifest match. That is the same failure
   *shape* as item 5 (confident, well-cited, wrong) happening one session earlier in the
   pipeline than where Session 3 currently looks for it. If item 5 is worth a dedicated
   weakness entry because it fabricates with confidence, this deserves the same treatment,
   not just a rehearsal note filed separately.

3. **"Assertions nothing verifies" is broader than item 1's citation-checking scope — it
   recurs in the plan-review pipeline too, with no citation involved at all.** Item 1 is
   scoped to `handoff.md`'s `file:line` citations specifically. This session found the
   identical failure *shape* somewhere item 1 doesn't look: a plain factual claim ("8
   builtins share this structure") stated in one ticket's research doc, repeated in its
   plan, **surviving that plan's own adversarial `plan-reviewer` pass** (which checked the
   shape of the defect it was sent to find, not every incidental number riding along with
   it), then carried forward unverified into a second ticket's research and plan, and only
   caught when a second, independent review pass happened to re-derive the number from
   source by chance. Nothing here involves a citation at all — it's a stated fact that
   simply never got re-checked, twice, across two already-reviewed documents. Worth its own
   line: the pattern is "an asserted fact survives review because review checks what it was
   pointed at, not everything riding along with it," and citation-verification (item 1) is
   only the narrowest instance of it.

**Is anything on the list presented as smaller than it actually is?** One clear case, named
above: item 4's phrasing reads as "known, pinned, stable" when the live entry it's sourced
from says the count grew from 2 to 5 in one day. Everything else checked out at the size the
list already gives it — this was the one place where the live source told a slightly
different story than the summary in front of the room would.
