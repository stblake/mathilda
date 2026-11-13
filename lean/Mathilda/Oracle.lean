/-
  Transport to the Mathilda CAS (the untrusted oracle).

  `callMathilda q` writes a temporary Mathilda script

      Print["@@MATHILDA-BEGIN@@"]; Print[q]; Print["@@MATHILDA-END@@"]

  runs `Mathilda -file <script>` (stdin is null, so Mathilda never enters its
  interactive or NDJSON pipe mode), returns the text between the sentinels and
  deletes the script.  The sentinels keep the answer apart from anything else
  that may reach stdout (a user init file, a future module's `Print`).

  Nothing returned here is trusted: every tactic in `Mathilda.Tactics`
  re-checks the answer in the Lean kernel.

  The binary is, in order:
  1. the environment variable `MATHILDA_BIN` (when non-empty);
  2. the option `set_option mathilda.path "/path/to/Mathilda"` (when non-empty);
  3. a `Mathilda` executable file in the working directory or its parent (the
     repository layout: the package lives in `<repo>/lean`, the binary is
     `<repo>/Mathilda`);
  4. `Mathilda` on `PATH`.
  A binary that cannot be started gives an error naming all of these.
-/
import Lean

namespace Mathilda

open Lean

register_option mathilda.path : String := {
  defValue := ""
  descr := "path of the Mathilda executable used by the mathilda_* tactics and #mathilda \
    (the environment variable MATHILDA_BIN takes precedence)"
}

def beginSentinel : String := "@@MATHILDA-BEGIN@@"
def endSentinel : String := "@@MATHILDA-END@@"

/-- Is `p` an existing regular file (not a directory such as `lean/Mathilda/`)? -/
def isFile (p : System.FilePath) : IO Bool := do
  try return (← p.metadata).type == .file
  catch _ => return false

/-- The Mathilda executable to run (see the module docstring for the order). -/
def mathildaPath (opts : Options) : IO String := do
  if let some p ← IO.getEnv "MATHILDA_BIN" then
    if p != "" then return p
  let o := mathilda.path.get opts
  if o != "" then return o
  for c in ["Mathilda", "../Mathilda"] do
    let p : System.FilePath := c
    if ← isFile p then return (← IO.FS.realPath p).toString
  return "Mathilda"

/-- A Mathilda string literal for `s` (`\` and `"` escaped). -/
def mathildaString (s : String) : String :=
  "\"" ++ (s.replace "\\" "\\\\" |>.replace "\"" "\\\"") ++ "\""

/-- Run the Mathilda expression `query` and return what it prints (its value as
    text).  An answer starting with `MATHILDA-FAIL:` becomes an error. -/
def callMathilda (query : String) (opts : Options := {}) : IO String := do
  let bin ← mathildaPath opts
  let notFound (why : String) : IO.Error := IO.userError s!"mathilda: cannot run the \
    Mathilda executable '{bin}' ({why}). Set the environment variable MATHILDA_BIN, or \
    `set_option mathilda.path \"/path/to/Mathilda\"`, or put Mathilda on PATH."
  if bin.contains '/' && !(← isFile bin) then throw (notFound "no such file")
  let (h, tmp) ← IO.FS.createTempFile
  h.putStr s!"Print[\"{beginSentinel}\"];\nPrint[{query}];\nPrint[\"{endSentinel}\"];\n"
  h.flush
  try
    let out ← try IO.Process.output { cmd := bin, args := #["-file", tmp.toString] }
      catch e => throw (notFound (toString e))
    -- A spawn failure surfaces as exit code 255 from the forked child on POSIX.
    if out.exitCode == 255 && out.stderr.startsWith "could not execute external process" then
      throw (notFound "not found on PATH")
    if out.exitCode != 0 then
      throw <| IO.userError s!"mathilda: '{bin}' exited with code {out.exitCode}: \
        {out.stderr.trimAscii.toString}"
    match out.stdout.splitOn (beginSentinel ++ "\n") with
    | [_, rest] =>
      match rest.splitOn ("\n" ++ endSentinel) with
      | ans :: _ :: _ =>
        let ans := ans.trimAscii.toString
        if ans.startsWith "MATHILDA-FAIL:" then
          throw <| IO.userError s!"mathilda: Mathilda found no certificate: \
            {(ans.drop 14).trimAscii.toString}"
        return ans
      | _ => throw <| IO.userError s!"mathilda: no end sentinel in Mathilda's output \
          (stderr: {out.stderr.trimAscii.toString})"
    | _ => throw <| IO.userError s!"mathilda: no begin sentinel in Mathilda's output \
        (stderr: {out.stderr.trimAscii.toString})"
  finally
    IO.FS.removeFile tmp

end Mathilda
