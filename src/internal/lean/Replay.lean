/-
  Mathilda kernel-replay checker (run as `lean --run Replay.lean <file>.olean`).

  `Certify` and `LeanVerify` compile the Lean source to an .olean in a private
  temp directory and then run this script on it, in a fresh `lean` process
  that executes none of the checked file's code.  It

  1. imports the module's imports (core, Mathlib, ... -- the compiled
     libraries are trusted as built);
  2. replays every declaration of the module into a fresh kernel environment
     (`Lean.Kernel.Environment.replay`, the per-module mode of `leanchecker`),
     so a declaration that the elaborator added without a kernel check
     (`debug.skipKernelTC`, `addDecl` from a user tactic, an environment hack)
     is re-checked and rejected;
  3. collects, from the replayed environment, the axioms every declaration of
     the module depends on.

  The verdict therefore never comes from the checked file's own output (a
  user elaborator can print anything, including a fake `#print axioms`).

  Output (stdout), parsed by src/lean_run.c:

      MATHILDA-REPLAY-BEGIN
      C <name>          one line per declaration of the module, followed by
      A <axiom>         one line per axiom it depends on
      MATHILDA-REPLAY-OK            (or: MATHILDA-REPLAY-FAIL <reason>)

  Names are printed with `Name.toString`; a newline in a name is replaced by
  `?` (such a name then matches nothing, which fails closed).  Exit code 0
  only on success.  Unsafe and `partial` declarations are not replayed (the
  kernel does not accept them); anything that uses one fails the replay.
-/
import Lean.Replay
import Lean.Util.CollectAxioms
import Lean.Util.Path

open Lean

namespace MathildaReplay

/-- `MonadEnv` for collecting the axioms of imported constants. -/
instance : MonadEnv (StateM Environment) := { getEnv := get, modifyEnv := modify }

/-- Walk state: constants visited for the current root, the axioms found so
    far, and a cache of the (precomputed) axioms of imported constants that
    is shared by all roots. -/
structure St where
  visited : NameSet := {}
  axioms : NameSet := {}
  imported : NameMap (Array Name) := {}

/-- Axioms reachable from `n`: declarations of the module are walked in the
    replayed kernel environment `k`; imported ones use Lean's precomputed
    axiom table of the imported environment `env`. -/
partial def axiomsOf (k : Kernel.Environment) (env : Environment)
    (isNew : Name → Bool) (n : Name) : StateM St Unit := do
  if (← get).visited.contains n then return
  modify fun s => { s with visited := s.visited.insert n }
  if isNew n then
    match k.find? n with
    | some ci =>
      if let .axiomInfo _ := ci then modify fun s => { s with axioms := s.axioms.insert n }
      for u in ci.type.getUsedConstants do axiomsOf k env isNew u
      if let some v := ci.value? (allowOpaque := true) then
        for u in v.getUsedConstants do axiomsOf k env isNew u
      if let .inductInfo v := ci then
        for c in v.ctors do axiomsOf k env isNew c
    | none =>
      -- Not in the replayed environment (an unsafe / partial declaration):
      -- report a pseudo-axiom, which is never standard.
      modify fun s => { s with axioms := s.axioms.insert (`mathilda.notReplayed ++ n) }
  else
    let axs ← match (← get).imported.find? n with
      | some axs => pure axs
      | none => do
        let axs := (collectAxioms (m := StateM Environment) n).run' env
        modify fun s => { s with imported := s.imported.insert n axs }
        pure axs
    modify fun s => { s with axioms := axs.foldl NameSet.insert s.axioms }

def clean (s : String) : String := s.map fun c => if c == '\n' || c == '\r' then '?' else c

end MathildaReplay

open MathildaReplay in
unsafe def main (args : List String) : IO UInt32 := do
  let [olean] := args
    | IO.println "MATHILDA-REPLAY-FAIL usage: lean --run Replay.lean <file>.olean"; return 2
  try
    -- The sysroot of the `lean` running this script, plus LEAN_PATH (set by
    -- `lake env` for a Mathlib project).
    initSearchPath ((← IO.appDir) / "..")
    let mFile : System.FilePath := olean
    let mut fnames := #[mFile]
    let sFile := OLeanLevel.server.adjustFileName mFile
    if ← sFile.pathExists then
      fnames := fnames.push sFile
      let pFile := OLeanLevel.private.adjustFileName mFile
      if ← pFile.pathExists then fnames := fnames.push pFile
    let parts ← readModuleDataParts fnames
    let some (mod, _) := parts[0]? | throw <| IO.userError "cannot read the module data"
    let some (last, _) := parts[parts.size - 1]? | throw <| IO.userError "no module data"
    let (_, s) ← importModulesCore mod.imports |>.run
    let env ← finalizeImport s mod.imports {} 0 false false (isModule := true)
    let mut newConstants : Std.HashMap Name ConstantInfo := {}
    for name in last.constNames, ci in last.constants do
      newConstants := newConstants.insert name ci
    let kenv ← env.toKernelEnv.replay newConstants
    IO.println "MATHILDA-REPLAY-BEGIN"
    let mut cache : NameMap (Array Name) := {}
    for n in last.constNames do
      let (_, st) := (axiomsOf kenv env newConstants.contains n).run { imported := cache }
      cache := st.imported
      IO.println s!"C {clean n.toString}"
      for a in st.axioms.toList do IO.println s!"A {clean a.toString}"
    IO.println "MATHILDA-REPLAY-OK"
    return 0
  catch e =>
    IO.println s!"MATHILDA-REPLAY-FAIL {clean (toString e)}"
    return 1
