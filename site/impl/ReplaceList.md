---
source: src/replace.c
---
**Algorithm.** `builtin_replacelist` returns *all* ways a rule (or rule list) can match `expr`, rather than the first. It flattens the rule argument into a `ReplaceRule[]` array (recording the `delayed` flag for `RuleDelayed`), then for each rule runs the pattern matcher with a callback installed on the `MatchEnv`. Instead of stopping at the first successful binding, the matcher's backtracking enumerates every distinct binding (notably every partition of a `BlankSequence`/`BlankNullSequence`), and `replacelist_callback` materialises each one: it builds the replacement via `replace_bindings`, additionally `eval_and_free`-ing it for delayed rules, and appends it to a growable `ReplaceListState.results` buffer. An optional third argument caps the number of results (`state.limit`), checked both before invoking the matcher's callback and between rules. The accumulated results are wrapped in a `List`.

**Data structures.** `ReplaceRule { pattern, replacement, delayed }` and `ReplaceListState { results, count, cap, limit, replacement, delayed }` — the matcher signals each match through `env->callback`/`env->callback_data`.
