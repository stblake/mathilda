(* Standalone PARANOIA runner.  From the repository root:
       ./Mathilda -file tests/paranoia/run_paranoia.m
   Runs every mode (exact comparisons, gated) plus "Machine" with the
   language's own comparisons (informational), prints each report and
   the gate verdict.  tests/test_paranoia.c is the ctest equivalent. *)
Get["tests/paranoia/paranoia.m"];
Paranoia`ParanoiaGate[];
