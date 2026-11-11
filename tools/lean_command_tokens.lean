import Lean
/-!
# lean_command_tokens.lean

Trusted helper of Mathilda's test suite (never given to LeanVerify).  Prints
the facts about Lean's parser that LeanVerify's source check relies on:

* `lead TOKEN`  every first token of every leading parser of the `command`
                category (`Lean.Parser.parserExtension`, the leading table of
                the category plus its non-indexed leading parsers);
* `nonindexed N` the number of leading command parsers that are not indexed
                by a token (they could start at an identifier; must be 0);
* `token TOKEN` every entry of the token table.

`tests/scripts/lean_check_tokens.sh` runs it (with `import Lean`, or with
`import Mathlib` in a Lake project) and checks every token against
`Lean`Internal`CheckSource`; `src/lean_check_tokens.h` is generated from the
`lead` words.  Run it only through that script.
-/
open Lean Parser

def firstToks (p : Parser) : List String :=
  match p.info.firstTokens with
  | .tokens ts => ts
  | .optTokens ts => ts
  | _ => []

#eval show CoreM Unit from do
  let env ← getEnv
  let st := parserExtension.getState env
  let some cat := st.categories.find? `command | throwError "no command category"
  let mut out : Array String := #[]
  for (_, ps) in cat.tables.leadingTable.toList do
    for (p, _) in ps do
      for t in firstToks p do out := out.push ("lead " ++ t)
  for (p, _) in cat.tables.leadingParsers do
    for t in firstToks p do out := out.push ("lead " ++ t)
  out := out.push ("nonindexed " ++ toString cat.tables.leadingParsers.length)
  for t in (getTokenTable env).values do out := out.push ("token " ++ t)
  for l in (out.qsort (· < ·)).toList.eraseDups do IO.println l
