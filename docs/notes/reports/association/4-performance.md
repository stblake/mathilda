## 4. Performance vs Mathematica 15 and Python

**Summary.** On bulk operations Mathilda keeps pace with Mathematica 15 and is often faster: KeyDrop, Keys, Values, KeySort, Merge, Normal and bulk Lookup take 0.06–0.6× Mathematica's time. Every bulk operation scales linearly. Single-key lookup costs O(1) and matches Mathematica. The one serious failure is **incremental update**. An `AssociateTo` loop or an `a[k] = v` loop is **O(n²)** in Mathilda, against O(n) in both reference systems. At n = 10⁴ that makes Mathilda 300–1400× slower than Mathematica, and a 10⁵-entry build is out of reach: we extrapolate about 10 minutes. The other two clear losses are `PositionIndex` (18× slower) and `GroupBy` (1.7× slower). Python's `dict` is 5–40× faster than both CAS systems on almost everything.

### 4.1 Setup

| | |
|---|---|
| Machine | Apple M5 Max (18 cores), 128 GB RAM, macOS 26.6.1 |
| Mathilda | v0.187 (`src/version.h`), the `./mathilda` binary built 2026-09-24, run as `./mathilda file.m` |
| Mathematica | 15.0.0 for Mac OS X ARM (64-bit), May 19 2026, run through `wolframscript -file` |
| Python | CPython 3.14.7: `dict`, `collections.Counter`/`defaultdict`, `itertools.groupby` |

**Method.**
- **Same script in both CAS systems.** Mathilda and Mathematica run the *same* `.m` source (`/tmp/assocreport/perf-bench/body.m` and `loopbody.m`). Its header sets only `n`.
- **Warm medians.** Each case evaluates once to warm up, which also captures the result. It is then timed with `AbsoluteTiming[e;]`: the median of 5 runs, or 3 runs at n = 10⁶ and for the insert loops.
- **Timing scope.** Times are in-process, so neither system's startup is included.
- **Python.** `py_bench.py` and `py_loop.py` mirror each case with `time.perf_counter` and the same warm-up and median rules.
- **Data.**
  - Keys are `Range[n]`. Values are `vs = Mod[k·7919, 1000]`, which gives 1000 distinct values.
  - The string variants use keys `"k"<>ToString[i]` and values `"v"<>ToString[v]`.
  - Merge and Join combine `a` with a copy whose keys are shifted by n/2, which gives 50% overlap.
  - Bulk Lookup uses 10⁵ pseudo-random keys, `Mod[i·104729, n] + 1`.
- **Value checks.** Every case prints a checksum next to its timing: Length, Total of values, the sum of first and last key, and so on. **All 64 (operation, n) triples agree across all three systems.** Every system is therefore doing the same work.

### 4.2 Results (milliseconds, warm median)

| Operation | n | Mathilda | Mathematica | Python | Mathilda/Mathematica ratio |
|---|---:|---:|---:|---:|---:|
| Construct `Association[rules]` | 10³ | 0.086 | 0.069 | 0.018 | 1.25 |
| | 10⁵ | 22.4 | 22.5 | 1.70 | 1.00 |
| | 10⁶ | 255 | 363 | 19.5 | **0.70** |
| Construct `AssociationThread` (int keys) | 10³ | 0.093 | 0.071 | 0.021 | 1.31 |
| | 10⁵ | 29.6 | 22.4 | 1.84 | 1.32 |
| | 10⁶ | 337 | 333 | 20.7 | 1.01 |
| Construct `AssociationThread` (string keys) | 10³ | 0.097 | 0.067 | 0.040 | 1.45 |
| | 10⁵ | 29.8 | 20.2 | 4.10 | 1.48 |
| | 10⁶ | 319 | 288 | 130 | 1.11 |
| Single-key `a[k]`, 10⁵ probes in a `Do` loop (int) | 10³ | 41.7 | 30.5 | 5.35 | 1.37 |
| | 10⁵ | 55.6 | 55.0 | 4.39 | 1.01 |
| | 10⁶ | 60.8 | 51.0 | 4.56 | 1.19 |
| Single-key `a[k]`, 10⁵ probes (string) | 10³ | 45.6 | 40.0 | 5.17 | 1.14 |
| | 10⁵ | 65.5 | 68.6 | 9.82 | 0.96 |
| | 10⁶ | 59.9 | 62.5 | 11.4 | 0.96 |
| `Lookup[a, 10⁵ keys]` | 10³ | 2.19 | 2.80 | 2.52 | 0.78 |
| | 10⁵ | 14.2 | 22.9 | 5.96 | 0.62 |
| | 10⁶ | 19.2 | 31.7 | 10.4 | **0.61** |
| `AssociateTo[b, i->2i]` loop, n inserts | 10³ | 31.7 | 0.294 | 0.060 | **108** |
| | 3·10³ | 322 | 1.03 | 0.182 | **311** |
| | 10⁴ | 5 481 | 4.00 | 0.598 | **1 371** |
| | 2·10⁴ | 23 688 ¹ | — | — | — |
| | 10⁵ | not feasible ² | 39.7 | 6.52 | — |
| `b[i] = 2i` loop, new keys | 10³ | 13.5 | 0.518 | 0.024 | **26** |
| | 3·10³ | 155 | 1.96 | 0.073 | **79** |
| | 10⁴ | 2 571 | 8.64 | 0.238 | **298** |
| | 2·10⁴ | 12 079 ¹ | — | — | — |
| | 10⁵ | not feasible ² | 106 | 2.90 | — |
| `b[i] = 2i` loop, overwrite existing keys | 10³ | 25.7 | 0.461 | 0.036 | **56** |
| | 10⁴ | 5 204 | 6.51 | 0.376 | **800** |
| | 10⁵ | not feasible ² | 91.4 | 4.26 | — |
| `KeyDrop[a, n/2 keys]` | 10³ | 0.017 | 0.065 | 0.034 | 0.26 |
| | 10⁵ | 7.50 | 35.0 | 3.03 | 0.21 |
| | 10⁶ | 83.9 | 453 | 39.7 | **0.19** |
| `Keys[a]` | 10³ | 0.008 | 0.024 | 0.003 | 0.33 |
| | 10⁵ | 1.52 | 10.1 | 0.221 | 0.15 |
| | 10⁶ | 12.5 | 199 | 2.59 | **0.06** |
| `Values[a]` | 10³ | 0.008 | 0.025 | 0.003 | 0.32 |
| | 10⁵ | 1.47 | 10.9 | 0.258 | 0.13 |
| | 10⁶ | 12.8 | 109 | 3.00 | **0.12** |
| `Map[# + 1 &, a]` | 10³ | 0.231 | 0.219 | 0.039 | 1.05 |
| | 10⁵ | 70.8 | 59.3 | 2.90 | 1.19 |
| | 10⁶ | 732 | 597 | 32.1 | 1.23 |
| `KeySort` (reverse-ordered input) | 10³ | 0.053 | 0.123 | 0.032 | 0.43 |
| | 10⁵ | 19.7 | 42.4 | 5.81 | 0.46 |
| | 10⁶ | 327 | 671 | 65.0 | **0.49** |
| `Merge[{a, a2}, Total]` (2n entries, 50% overlap) | 10³ | 0.491 | 1.08 | 0.131 | 0.45 |
| | 10⁵ | 116 | 273 | 13.4 | 0.42 |
| | 10⁶ | 1 352 | 3 313 | 135 | **0.41** |
| `Join[a, a2]` | 10³ | 0.158 | 0.124 | 0.015 | 1.27 |
| | 10⁵ | 48.9 | 44.4 | 1.37 | 1.10 |
| | 10⁶ | 591 | 701 | 16.4 | 0.84 |
| `GroupBy[vs, Mod[#,100]&]` | 10³ | 0.176 | 0.168 | 0.040 | 1.05 |
| | 10⁵ | 23.3 | 20.2 | 2.92 | 1.15 |
| | 10⁶ | 287 | 166 | 29.6 (153 ³) | **1.72** |
| `Counts[vs]` (int, packed) | 10³ | 0.084 | 0.171 | 0.021 | 0.49 |
| | 10⁵ | 0.286 | 0.219 | 2.25 | 1.31 |
| | 10⁶ | 1.28 | 0.815 | 29.1 | 1.57 |
| `Counts` (string values) | 10³ | 0.080 | 0.151 | 0.028 | 0.53 |
| | 10⁵ | 0.956 | 2.86 | 2.71 | 0.33 |
| | 10⁶ | 7.39 | 16.6 | 25.4 | **0.45** |
| `PositionIndex[vs]` | 10³ | 0.139 | 0.212 | 0.128 | 0.66 |
| | 10⁵ | 7.94 | 1.02 | 5.27 | **7.8** |
| | 10⁶ | 64.6 | 3.65 | 48.3 | **17.7** |
| `Normal[a]` | 10³ | 0.011 | 0.076 | 0.015 | 0.14 |
| | 10⁵ | 2.68 | 20.0 | 2.76 | 0.13 |
| | 10⁶ | 48.9 | 187 | 28.6 | **0.26** |

¹ Single run rather than a median (`one.m`).
² Extrapolated from the measured n² scaling: about 10 minutes per run, and a median needs four runs.
³ With `itertools.groupby(sorted(...))` in place of `defaultdict(list)`. The sort makes it 5× slower than the dict idiom.

Python's idiomatic equivalents for each row:

| Mathilda operation | Python equivalent |
|---|---|
| `Association`, `AssociationThread` | `dict(zip(...))` |
| `Lookup` of many keys | a list comprehension |
| `KeyDrop` | a dict comprehension with a `set` |
| `Merge[..., Total]` | `defaultdict(int)` |
| `Join` | `a \| a2` |
| `GroupBy`, `PositionIndex` | `defaultdict(list)` |
| `Counts` | `Counter` |
| `AssociateTo` | `dict.update` |
| `Normal` | `list(a.items())` |

### 4.3 Scaling observations

We computed the scaling exponent as log(t₂/t₁) / log(n₂/n₁). The 10⁵→10⁶ step is the reliable one, because at n = 10³ fixed call overhead dominates.

- **Every bulk operation is linear.**
  - Mathilda's exponents from 10⁵ to 10⁶ are 1.01–1.08 for construction, Map, Merge, Join, KeyDrop and GroupBy.
  - Keys and Values come out at 0.91–0.94.
  - KeySort is 1.22 and Normal 1.26. Both are consistent with n log n and with allocator and cache pressure at 10⁶.
  - Mathematica's exponents on the same steps are 0.92–1.29. Python's are 1.00–1.12.
  - The hash-indexed design stated at `src/assoc.c:1-9` holds up: no bulk path has quadratic behaviour.
- **Single-key lookup is O(1).** The time for 10⁵ probes is flat in n: 55.6 ms at n = 10⁵ and 60.8 ms at 10⁶, an exponent of 0.04. Mathematica shows the same flat profile (−0.03). The lazy persistent index (`src/assoc.c:297-322`, `src/assoc_index.c:38-77`) survives across `Do` iterations, as `tests/bench_assoc.c` intends. Most of those ~55 ms is interpreter `Do`-loop overhead in both CAS systems, about 0.55 µs per iteration. Python runs the same loop 12× faster.
- **Bulk Lookup is O(n + m).** The exponent is 0.13 at a fixed m = 10⁵, which is simply the one-off index build (`src/assoc.c:408-431`).
- **Counts on a packed integer list is sub-linear in n** (0.65), because it only walks the 1000 distinct values after a machine-word Tally (`src/assoc.c:782`, and `Counts` is on the `AWARE` list in `src/pack.c:812`). That makes it 23× faster than Python's `Counter` at 10⁶.
- **Super-linear: every incremental-update loop.**
  - Mathilda exponents are 2.11–2.35 across 10³→3·10³→10⁴, and 2.11 from 10⁴ to 2·10⁴ for `AssociateTo`.
  - `a[k] = v` has the same profile, whether it appends new keys or overwrites existing ones.
  - Mathematica stays at 0.92–1.32 up to 10⁵, and Python at about 1.0.
  - This is exactly the "O(n²) behaviour from repeated copying" failure mode. Each insert costs O(n), about 0.5 µs per existing entry at n = 10⁴, so the gap widens without bound.

**Memory (peak RSS from `/usr/bin/time -l`).**

| Measurement | Result |
|---|---|
| Mathilda, empty script | 17.6 MB |
| Mathilda, `Range[10⁶]` plus a values list | 41.8 MB |
| Mathilda, the same plus one 10⁶-entry `AssociationThread[ks, vs]` | **524 MB** |
| Mathilda, marginal cost of that association | ≈ **482 MB** |
| Mathilda, string keys | 527 MB |
| Mathilda, three coexisting 10⁶-entry associations | 1 087 MB |
| Mathilda, steady-state cost per entry | ≈ **350 B** |
| Mathematica, same association: `ByteCount` | 126 MB (126 B/entry) |
| Mathematica, same association: `MemoryInUse` growth | 110 MB |
| Python, same dict: RSS growth | 82 MB (`dict` object 42 MB) |
| Mathilda, full 10⁶ benchmark script (all cases) | 3.2 GB |

Mathilda's association is therefore about **2.8× Mathematica's and 4× Python's** footprint. The cause is that each entry holds four separately malloc'd nodes:
- a `Rule` node;
- its 2-slot argument array;
- a *fresh* `Rule` symbol node from `make_rule` (`src/assoc.c:90-93`);
- the key.

Each `Expr` carries a 16-byte header of `last_evaluated_at` and `hash_cache` (`src/expr.h`).

### 4.4 The three biggest losses

**1. `AssociateTo` / `a[k] = v` loops are O(n²): 300–1400× slower than Mathematica at n = 10⁴, and worsening linearly with n.**

Every single-key update rebuilds the entire association:
- **`AssociateTo`** (`src/assoc.c:1057-1089`):
  - It evaluates the symbol (`:1063`).
  - It gathers all `base` existing rules plus the new one into a fresh array (`:1067-1080`).
  - It re-canonicalises through `assoc_from_rules` (`:1082`), which allocates a new `Rule` node, a new `Rule` symbol and new copies for **every** entry (`:151`, `make_rule` at `:90-93`) and builds a full hash table (`:135`).
  - It writes the result back through `symtab_add_own_value` (`:1087`).
- **`a[k] = v`** is rerouted by `apply_assignment` into Part assignment (`src/eval.c:1107-1133`). `expr_part_assign` then:
  - calls `evaluate(sym)` (`src/part.c:413`);
  - copies all n argument pointers into a new array (`src/part.c:220-223`);
  - finds the key by **linear scan** (`src/part.c:283-287`), ignoring the persistent index;
  - appends via `realloc` (`:295`) and builds a new node.
- **Profile.** A `sample` profile of the 10⁴-key `a[i] = 2i` loop attributes about 87% of samples to `expr_part_assign → evaluate → builtin_association → assoc_from_rules`, with `malloc`/`free` at the top of the stack.
  - The Association node that has just been reassigned is re-evaluated, and `builtin_association` rebuilds a complete canonical copy only to discover it was already canonical and throw it away (`src/assoc.c:219-226`).
  - Because each update yields a new node, the lazily built key index (`src/assoc.c:299-317`) is also discarded and rebuilt on the next read.

Mathematica mutates in place when the value is uniquely referenced, which gives amortised O(1) per update.

**Hypothesis.** An in-place fast path when `refcount == 1` would remove the quadratic cost:
- append, or overwrite `args[i]`, and update the attached `AssocIndex` incrementally, which needs a growable variant of `assoc_index.c`;
- skip the canonicality re-check on a node the assoc code itself produced, for example by marking it canonical or ground.

None of the regression gates in `tests/bench_assoc.c` exercises an update loop; its 9 bulk operations plus the lookup gates at `:147-157` include no update. That is why this went unnoticed.

**2. `PositionIndex` is 18× slower than Mathematica at 10⁶ (64.6 ms vs 3.65 ms).**

`builtin_positionindex` (`src/assoc.c:1416-1449`) accepts only a boxed `List` (`:1419`), and `PositionIndex` is not on the `AWARE` list in `src/pack.c`. As a result:
- **A packed input is first materialised** by the transparency gate into 10⁶ boxed `Expr`s. Measured separately (`pi2.m`), this costs 52 ms packed against 31 ms on an already-unpacked list, so about 20 ms is pure unpacking.
- **Every element pays** a `ki_lookup` over a boxed `Expr`, which costs an `expr_hash` and an `expr_eq`.
- **Every position is a separate heap `Integer`** (`expr_new_integer`, `:1438`).
- **Every group grows** through its own `realloc`'d array (`:1437`).

Mathematica evidently returns packed integer position vectors: 5 ms packed, 15 ms unpacked.

**Hypothesis.** An int64 buffer kernel in the style of `counts_from_ndarray` (`src/assoc.c:704-757`), which already makes `Counts` 23× faster than Python, would close this gap. It would direct-index or hash machine words and emit packed position lists. That also satisfies the project's own packed-array rule for numeric heads.

**3. `GroupBy` is 1.7× slower than Mathematica at 10⁶ (287 ms vs 166 ms).**

`Map` over values shares the same root cause, at 1.2× slower.

`builtin_groupby` (`src/assoc.c:820-903`) builds and fully evaluates a fresh `f[x]` call for every element (`:854-857`). That is one `expr_new_function`, one evaluator dispatch and one `Function` application per element, with no attempt to compile or vectorise `f` over the (materialised) packed input.

In isolation, `Map[Mod[#,100]&, list]` costs 209 ms in Mathilda against 9 ms in Mathematica, which compiles pure functions over packed data. Mathematica's own GroupBy does not fully exploit that either, which is why the gap is only 1.7×.

`assoc_map_values` (`src/assoc.c:1105-1122`) likewise emits 10⁶ unevaluated `f[v]` nodes for the evaluator to reduce one at a time. The resulting association is then re-canonicalised through `builtin_association`/`assoc_from_rules` (`src/assoc.c:201-235`). The same double canonicalisation (a rebuild, then a discarded rebuild) is visible in the profile of `Association[rules]`, where `expr_free` of the discarded copy alone takes about 14% of samples. It inflates the constant factor of every operation that produces an association, and helps explain why Join and construction trail Python by 16–36×.

**Hypothesis.** Evaluating `f` over the key-input vector once, either through the auto-compiler (`src/compile/autocompile.c`) or through a packed `Map`, and then hashing the resulting machine keys would turn GroupBy into an O(n) hash pass. Separately, having `builtin_association` trust nodes built by `assoc_*` constructors would remove one full O(n) allocate/free cycle from every association-producing builtin.

### Key findings

- **[strength]** Every bulk association operation scales linearly from 10⁵ to 10⁶ (exponents 0.9–1.26), and all 64 value checks agree with Mathematica and Python. Mathilda beats Mathematica 15 on KeyDrop (5×), Keys (16×), Values (8.5×), Normal (3.8×), Merge (2.4×), KeySort (2×), string Counts (2.2×) and bulk Lookup (1.6×) at n = 10⁶.
- **[strength]** Single-key lookup is genuinely O(1): 10⁵ probes cost about 55–60 ms from n = 10³ to 10⁶, on par with Mathematica. Packed-integer `Counts` (1.3 ms at 10⁶) is 23× faster than Python's `Counter`.
- **[risk]** `AssociateTo` and `a[k] = v` in a loop are O(n²): about 5.5 s for 10⁴ inserts against 4 ms in Mathematica, 23.7 s at 2·10⁴, and an estimated 10 minutes at 10⁵. Each update copies and re-canonicalises the whole association (`src/assoc.c:1067-1087`; `src/part.c:413`, `:220-223`, `:283-287`). Incremental dictionary building is the most common association idiom in user code.
- **[gap]** `tests/bench_assoc.c` gates 9 bulk operations plus lookup but no update loop, `Join`, `PositionIndex` or `Normal`, so the quadratic update path is invisible to CI. `benchmarks/` has no association experiment at all.
- **[gap]** `PositionIndex` (18× slower than Mathematica) and `GroupBy` (1.7×) have no packed or compiled path. `PositionIndex` is absent from `AWARE` in `src/pack.c` and boxes every position (`src/assoc.c:1419`, `:1438`).
- **[risk]** Memory: a 10⁶-entry association costs about 480 MB peak and about 350 B per entry at steady state. That is 2.8× Mathematica's 126 B and 4× Python's, driven by four mallocs per entry, including a fresh `Rule` symbol node per rule (`src/assoc.c:90-93`). The full 10⁶ benchmark peaked at 3.2 GB RSS.
- **[risk]** Every association-producing builtin pays a hidden second O(n) pass: `builtin_association` rebuilds an already-canonical node and then frees it (`src/assoc.c:219-226`). This inflates constant factors, especially for Map (1.2× slower than Mathematica) and Join and construction (16–36× slower than Python).
- **[gap]** Python's `dict` remains 5–40× faster than both CAS systems on construction, Join, Map and GroupBy. The single-key loop gap (12×) is interpreter overhead shared with Mathematica, not an association defect.

*Reproduction:* scripts and raw outputs are in `/tmp/assocreport/perf-bench/`:
- `body.m` and `loopbody.m`: prefix with `n = …;`;
- `py_bench.py` and `py_loop.py`;
- `mem*.m` for the RSS runs;
- `pi2.m` for the packed/unpacked split.
