/-
  Mathilda: Lean 4 tactics backed by the Mathilda computer algebra system.

  Core-only (no Mathlib).  Mathilda is an untrusted oracle; every certificate
  it returns is re-checked by the Lean kernel.  See README.md.
-/
import Mathilda.Prelude
import Mathilda.Serialize
import Mathilda.Oracle
import Mathilda.Tactics
