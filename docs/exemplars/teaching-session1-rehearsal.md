# Session 1 Rehearsal: `docs/TEACHING.md` Walked Against This Repo

**Kit version actually run against: `9.0.7`, commit `288bdd55f7076bfc9a993225e80affcb88a8e738`**
(`ms-bain/ai-sdlc-starterkit`, fresh `git clone`, not the plugin cache — the cache install
in this session remains pinned at `8.0.0`, per `KIT-FEEDBACK-GRAPH.md` GR-01/GR-15). A
second live re-check against whatever commit is current now was attempted before writing
this up and was blocked by a transient tool-classifier outage (`claude-sonnet-5[1m]
temporarily unavailable`) on every `Bash` call for several minutes running; not retried
further once the findings below were already complete and reproducible against a real,
verified commit. If the kit has moved again, that is expected — the version pinned above is
exact and everything below was run against it, not against a summary of it.

**What this is not**: a fix, a PR, or a judgment about the kit's overall trajectory. Per the
request that produced this: record what happened, do not repair it. The kit repo itself was
treated as read-only throughout.

---

## Session 1, as written, walked step by step

Quoting `docs/TEACHING.md`'s own Session 1 section throughout, then recording what actually
happened on this repo.

### Step 1 — the opener: bare-Makefile target parsing

**What the facilitator is told to do**: "Start with `/ais:setup-kit` against a real
repository the team already owns — ideally one with no CMakeLists, no obvious build
tooling... watch the tool propose exactly those commands back."

**Caveat on repo fit, stated up front**: this repo does not perfectly match the "ideally...
no CMakeLists" case the script describes — it has a real, primary top-level `makefile` *and*
a real, secondary `tests/CMakeLists.txt` (a CMake-based unit-test build, not the project's
main build). That combination — common in real C/C++ repos, where a hand-written Makefile
drives the main build and CMake is used only for a test target — turns out to be exactly
where this demo step breaks, as below.

**What actually happened**, running `detect_ladder.py --repo . --json` fresh (no
`--existing`, matching a first-time client demo) against this repo:

```json
{
  "toolchains_detected": ["c-cmake", "node-ts", "python"],
  "manifests": {
    "python": ["benchmarks/requirements.txt"],
    "node-ts": ["frontend/package.json"],
    "c-cmake": ["tests/CMakeLists.txt"]
  },
  "proposals": {
    "static": [
      {"toolchain": "c-cmake", "command": "clang-tidy -p build $(git ls-files '*.c' '*.cc' '*.cpp' '*.h' '*.hpp')", "source": "hint:c-cmake (manifest: tests/CMakeLists.txt)"},
      {"toolchain": "c-cmake", "command": "cppcheck --error-exitcode=1 .", "source": "hint:c-cmake (manifest: tests/CMakeLists.txt)"},
      ...
    ],
    "typecheck": [
      {"toolchain": "c-cmake", "command": "cmake --build build", "source": "hint:c-cmake (manifest: tests/CMakeLists.txt)"},
      ...
    ],
    "unit": [
      {"toolchain": "c-cmake", "command": "ctest --test-dir build --output-on-failure", "source": "hint:c-cmake (manifest: tests/CMakeLists.txt)"},
      ...
    ]
  },
  "detection_status": "ok"
}
```
(node-ts/python entries omitted above for length — `frontend/package.json` and
`benchmarks/requirements.txt` are real, tracked files in this repo; a genuine web frontend
and a Python benchmark harness exist alongside the C99 core, so those two detections are
correct, not noise.)

**This is real progress over 8.0.0** — at that version this exact repo produced
`toolchains_detected: []` across the board (this session's own GR-03 finding). At 9.0.7,
`detect_manifests()` correctly finds `tests/CMakeLists.txt` via a now-recursive shallow
search (the fix for GR-03/GR-15's root-only-manifest sub-bug, confirmed live). **But the
specific, celebrated moment — "watch the tool propose exactly those [Makefile] commands
back" — does not happen.** Not one of the `c-cmake` proposals above (`clang-tidy`,
`cppcheck`, `cmake --build build`, `ctest`) is a command this repo's own `makefile` defines
or that a maintainer of this repo would actually run. This repo's real, canonical build/test
commands (`make check-c99`, `make check-packed-aware`, `cd tests/build && cmake .. && make
&& for t in *_tests; do ./$t; done`) never appear anywhere in the output.

**Root cause, traced directly in `detect_ladder.py`, not inferred from behavior**:

1. `propose_from_makefile()` — the function that actually parses real Makefile target names
   via `_MAKE_TARGET_RE` and proposes `make <target>` commands, which is the mechanism
   Session 1's script is describing — is only ever called when
   **`if "c-cmake" not in manifests:`** (`detect_ladder.py:869`). Its own docstring gives the
   reasoning: *"a CMake project's own generated or checked-in Makefile is a build artifact
   CMake wrote, not a hand-authored signal of intent."* That reasoning is correct for a repo
   where CMake generated the Makefile. It is **wrong** for this repo's actual shape: the
   `tests/CMakeLists.txt` GR-03/GR-15 fix now (correctly) detects is a small, separate,
   test-only CMake project nested one directory down — it did not generate the real,
   hand-written, 600-line top-level `makefile`, which predates it and is the thing this
   project's own `SPEC.md` names as the canonical build. Finding *any* CMakeLists.txt
   anywhere in the tree is being read as "this whole repo is CMake's," which is false here.
2. Verified directly, bypassing the gate: `propose_from_makefile()` called on this repo in
   isolation returns `{}` anyway. `MAKE_PHASE_HINTS` only recognizes generic target names
   (`"lint"`, `"check-style"`, `"test"`, `"tests"`, `"check"`, `"typecheck"`, ...) — this
   repo's real targets (`check-c99`, `check-packed-aware`, `check-array-exactness`,
   `check-nd-surfaces`, ...) are all domain-specific compound names that match none of the
   hints. So even with the gate removed, this specific repo's Makefile still would not have
   produced the demo moment — a second, independent reason, not a duplicate of the first.

**A smaller, genuinely separate bug found while tracing this**: `parse_makefile_targets()`
produces a spurious literal `"\\"` in its target set for this repo. Root cause: this repo's
`.PHONY:` declaration spans three physical lines via trailing backslash continuation
(`makefile:600-602`); the parser only reads the first physical line of a `.PHONY:` block and
splits it on whitespace, so the line's own trailing `\` becomes a "target name." Low
severity — the targets listed only on the continuation lines (`check-array-exactness`,
`check-nd-surfaces`, etc.) are still captured correctly via their own separate `target:`
rule lines elsewhere in the file, so nothing is actually lost from the final proposal set —
but it is a real, reproducible data-quality defect in the parser, not a hypothetical.

**Falsifiable-in-the-room claim: does it hold?** **No, for this repo, for two independently
confirmed reasons** — not "it depends," not "sometimes." A room given this repo to try live
would watch the tool propose `cmake --build build` and `ctest --test-dir build`, then watch
someone in the room try `cmake --build build` from the repo root and get an immediate,
visible error (no root `CMakeLists.txt` exists to build). That is a worse outcome for the
"survives scrutiny" bar Session 1 explicitly sets for itself than an honest `FAILED` would
have been.

**Where a room would get confused rather than blocked**: this is the sharper finding than
"it doesn't work." A clean `FAILED` (see Step 2 below) is not confusing — it is Session 1's
own second-half demo, and the room is primed to expect and respect it. What actually happens
here is worse: `detection_status: "ok"`, five plausible-looking, correctly-sourced-and-cited
proposals, and a *wrong* toolchain confidently in front of the room. Nothing in the output
signals "these commands don't actually build this repo" — that only surfaces if someone in
the room actually runs one. A facilitator who does not already know this repo's real build
shape would have no reason to suspect anything before that moment, and the failure, when it
lands, lands as "the tool was wrong," not "the tool was honest about not knowing" — precisely
the distinction Session 1's own text says matters.

### Step 2 — the second half: a genuinely unrecognizable directory

**What the facilitator is told to do**: run detection again against "a directory with
genuinely nothing recognizable in it," and expect a loud `FAILED`, never a silent
`not-configured`-looking empty result.

**What actually happened**, against a freshly created empty directory:
```json
{
  "toolchains_detected": [], "manifests": {}, "proposals": {"static": [], "typecheck": [], "unit": [], "integration": []},
  "detection_status": "FAILED — no manifest, CI citation, recognized build file, or spec/docs hint matched anything in this repo. Checked for: *.tf, .flake8, CMakeLists.txt, Cargo.toml, Gemfile, build.gradle, build.gradle.kts, go.mod, mypy.ini, package.json, pom.xml, pyproject.toml, pytest.ini, requirements.txt, ruff.toml, setup.cfg, setup.py, tox.ini, Makefile, makefile, GNUmakefile. Searched the repo root and its immediate subdirectories... This is NOT the same as a correctly empty repo — it means this script has no name for whatever toolchain is actually here... Do not present this to a human as an ordinary not-configured answer."
}
```

**This half works exactly as advertised.** The message is loud, names precisely what it
checked, states its own search depth, and explicitly instructs against presenting it as a
benign result. Confirmed by execution, not by reading the source's stated intent. No
confusion risk here — this is a genuinely good, room-safe moment.

### Step 3 — `/ais:guide-me` cold

**What the facilitator is told to do**: run it with no arguments and expect one line naming
version, install path, and freshness before anything else.

**What actually happened**: `ais v9.0.7 — installed at /private/tmp/ais_demo (source:
unknown) — up to date with origin/main (as of last fetch)`. Works exactly as advertised — a
single, immediate, informative line. `(source: unknown)` is an artifact of testing against a
bare `git clone` rather than a real plugin install; a genuine client install would resolve
this to `plugin` or `git`, not a defect in the tool itself.

---

## Wall-clock, against the fifteen minutes Session 1 claims

Running the two commands Session 1 actually asks a facilitator to run (`detect_ladder.py`
twice, `version_info.py` once) takes seconds — well inside the fifteen-minute budget on
raw execution time alone, exactly as the doc implies. The fifteen minutes this session
actually spent went entirely into tracing *why* Step 1 didn't do what the script says it
would, which a live facilitator would not do in the room — they would either not notice
(the confusion-not-blocking risk above) or would notice and have no ready explanation,
since nothing in the tool's own output names the CMake-detection gate as the reason its
Makefile-target-reading branch never ran.

---

## The question that matters most: is the zero-signal case from 8.0.0 still this repo's case?

**No, partially — and the partial answer is the finding, not either extreme.**

`toolchains_detected` is no longer empty; `detect_manifests()` genuinely and correctly finds
real signal in this repo now (`tests/CMakeLists.txt`, `frontend/package.json`,
`benchmarks/requirements.txt`), which is real, verified progress since 8.0.0 and directly
attributable to the GR-03/GR-15 fix. The specific "five minutes when DETECT has nothing to
draft" framing this session wrote yesterday needs one more branch, not a retraction: DETECT
can now produce **a draft that exists but is wrong for this specific repo's actual build
system**, sourced from a real file, correctly cited, and confidently marked `"ok"` — which
is arguably a harder case for a live demo to recover from gracefully than the honest, loud
`FAILED` this repo used to produce at 8.0.0. The training material's durability argument
still holds either way (a human confirming or correcting a wrong draft is still cheaper long
-term than nothing at all, and still produces a durable, re-runnable, honest record) — but
"detection produces nothing" and "detection produces a confident, wrong, well-cited draft"
are two different demo failure modes, and this repo, today, at 9.0.7, is now firmly the
second one, not the first.
