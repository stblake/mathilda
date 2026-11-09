# G&R pipeline outputs (snapshot of the 2026-11-09 run)

These are the intermediate results of the Gradshteyn & Ryzhik pipeline
(`tools/gr/`, contract in `../TRANSCRIBE.md`), kept so that no stage has to be
re-run from the page images. The large JSON files are gzipped (`gzip -9 -n`);
run `gunzip -k` on one before you pass it to a script.

The data is split into three batches by chapter span:

| Directory | G&R chapters | Formulas transcribed | Records sent to the oracle |
|---|---|---|---|
| `ch2/` | 2 (indefinite, elementary) | 2336 | 1803 |
| `ch34/` | 3-4 (definite, elementary) | 4733 | 4319 |
| `ch5to17/` | 5, 6, 7 and 17 (special functions, transforms) | 2775 | 2491 |

Each directory holds these files:

- `entries.json.gz`: the output of `merge.py` over the batch's `../raw/p*.json`
  files. It contains every transcribed formula, including the skipped ones and
  the reason each was skipped.
- `mverified.json.gz`: the verdict of the Mathilda numeric oracle
  (`verify.py verify`, which runs `verify.m`). This is the latest re-verified
  run; `verify2.log` has its status counts.
- `py.json.gz` (and `py2.json.gz` in `ch2/`): the output of `pycheck.py`, a
  Python mpmath+scipy cross-check run on top of the oracle verdict.
  `build_corpus.py` reads this file. In `ch34/`, `py.json` was computed from
  the *earlier* `mverified`. Its re-run on the current `mverified` was cut off
  part-way, and `py2.json.cache.jsonl.gz` is that partial cache: un-gzip it to
  `py2.json.cache.jsonl` and `pycheck.py` resumes from it.
- `ch2/run.jsonl.gz` and `ch2/triage.json.gz`: the Mathilda `Integrate` run
  over the ch2 corpus (`tests/test_mathar_corpus.c` with `MATHAR_CORPUS_JSON`),
  and its `triage.py` buckets.
- `ch2/out/`: the `build_corpus.py` output for ch2. This is the corpus `.m`,
  `GRSkipped.md`, `GRNIntegrateGaps.md` and `GRNumericDefects.md`.
