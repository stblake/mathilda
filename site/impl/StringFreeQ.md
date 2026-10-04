---
source: src/strings/regex/stringcontainsq.c
---
**Algorithm.** `StringFreeQ` is the `SQ_FREE` face of the shared `sq_dispatch` core (`src/strings/regex/stringcontainsq.c`). The matching is identical to `StringContainsQ` — one `regex_match` per rule with an early exit — but `sq_answer` inverts the result, so `StringFreeQ` is `True` exactly when no rule matches anywhere. Option seeding, the curried operator form, and list-threading are all shared.

**Data structures.** As `StringContainsQ`: a `RegexRule` array, whole-match pair only, no capture pool.

**Complexity / limits.** `StringFreeQ[s, p] == !StringContainsQ[s, p]` by construction. An `IgnoreCase` option is honoured; a list of patterns means "free of any of them" (an empty list matches nothing). A non-string subject leaves the call unevaluated; byte semantics.
