(* solve_corpus_run.m
 *
 * Standalone driver for the Solve[] stress corpus, for fast local
 * iteration without rebuilding the C runner.  Run from the repository
 * root:
 *
 *     ./Mathilda -file tests/solve_corpus_run.m
 *
 * It loads the verifier prelude and the corpus, runs every case through
 * solveReport (which prints one `SOLVE<tab>PASS|FAIL|UNEVAL<tab>label`
 * line per case), and prints a final SOLVE-SUMMARY tally.  The C runner
 * tests/test_solve_corpus.c is the CI gate; this driver is the human
 * dashboard and mirrors its verdicts (modulo degenerate True/False
 * equations, which the C runner keeps held).
 *)

Get["tests/solve_check_prelude.m"];
cases = Get["tests/solve_corpus.m"];

npass = 0; nfail = 0; nuneval = 0;
Scan[
  Function[rec,
    Module[{code},
      code = solveReport @@ Take[rec, 5];
      Which[code === 0, npass++, code === 2, nuneval++, True, nfail++]]],
  cases];

Print["SOLVE-SUMMARY\ttotal=", Length[cases],
      "\tpass=", npass, "\tfail=", nfail, "\tuneval=", nuneval];
