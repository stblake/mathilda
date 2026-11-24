/* Builtin-audit regressions: file I/O (review/issues/<Head>.md; tests/test_audit.h).
 * Every file lives under $TemporaryDirectory and is deleted by the check that made it. */
#include "test_audit.h"

/* Export#1 (critical, crash): Export of a Graphics to PNG/JPEG segfaulted inside raylib's
 * InitWindow whenever the window or its GL context was unusable -- InitPlatform with a NULL
 * window under a sandbox that denies IOKit, or rlglInit -> rlLoadTexture jumping to a never-loaded
 * GL entry point (address 0) in a tool process that has a GUI session.  Nothing in-process predicts
 * either, so the first raster export now does the whole render once in a spawned child
 * (gl_probe.c) and declines with Export::nogl and $Failed if the child does not succeed.  In a
 * usable session the file is written.  Either outcome passes; a crash kills the suite.  (Checked by
 * hand under `sandbox-exec (deny iokit-open)`: SIGSEGV before, Export::nogl + $Failed after.) */
static void test_export_graphics_raster_no_crash(void) {
    audit_true("Export#1",
               "Module[{p = FileNameJoin[{$TemporaryDirectory, \"mth_audit_export1.png\"}], r, c},"
               " c = Quiet[Check[r = Export[p, Graphics[Disk[]]]; \"ok\", \"nogl\", Export::nogl]];"
               " If[FileExistsQ[p], DeleteFile[p]];"
               " (r === p && c === \"ok\") || (r === $Failed && c === \"nogl\")]");
    audit_true("Export#1",
               "Module[{p = FileNameJoin[{$TemporaryDirectory, \"mth_audit_export1.jpg\"}], r, c},"
               " c = Quiet[Check[r = Export[p, Graphics[Disk[]]]; \"ok\", \"nogl\", Export::nogl]];"
               " If[FileExistsQ[p], DeleteFile[p]];"
               " (r === p && c === \"ok\") || (r === $Failed && c === \"nogl\")]");
    /* sanity: a raster Image export (no renderer involved) still writes and round-trips */
    audit_true("Export#1",
               "Module[{p = FileNameJoin[{$TemporaryDirectory, \"mth_audit_export1b.png\"}], r, d},"
               " r = Export[p, Image[{{0., 1.}, {1., 0.}}]]; d = ImageData[Import[p]];"
               " DeleteFile[p]; r === p && d == {{0., 1.}, {1., 0.}}]");
}

/* Read#1 (high): Read[s, Number] parsed each token as Mathilda code and EVALUATED it --
 * Print["PWN"] in the file ran, 1/2 became a Rational and 2^3 became 8.  A Number is a
 * literal (sign, digits, point, e/E or *^ exponent); anything else is Read::readn and
 * $Failed, and is never evaluated. */
static void test_read_number_literal_only(void) {
    audit_run("audit$f = FileNameJoin[{$TemporaryDirectory, \"mth_audit_read1.txt\"}];"
              " audit$o = OpenWrite[audit$f];"
              " WriteString[audit$o, \"1/2 audit$side[] 2^3 -4 2.5 1.5e3 -2.*^-2 .5 7\\n\"];"
              " Close[audit$o]; audit$hit = 0; audit$side[] := (audit$hit = 1; 99)");
    audit_chk("Read#1",
              "Quiet[Read[audit$f, {Number, Number, Number, Number, Number, Number, Number, Number, Number}]]",
              "{$Failed, $Failed, $Failed, -4, 2.5, 1500.0, -0.02, 0.5, 7}");
    audit_chk("Read#1", "audit$hit", "0");                 /* the token was not run */
    audit_run("Close[audit$f]");
    audit_chk("Read#1", "Quiet[Check[Read[audit$f, Number], fired, Read::readn]]", "fired");
    audit_run("Close[audit$f]");
    audit_chk("Read#1", "Quiet[Read[audit$f, Real]]", "$Failed");
    audit_run("Close[audit$f]; DeleteFile[audit$f]");
    /* sanity: an ordinary numeric file reads as before */
    audit_chk("Read#1",
              "Module[{f = FileNameJoin[{$TemporaryDirectory, \"mth_audit_read1b.txt\"}], o, r},"
              " o = OpenWrite[f]; WriteString[o, \"10 -3.25\\n+6 1E2\\n\"]; Close[o];"
              " r = ReadList[f, Number]; DeleteFile[f]; r]",
              "{10, -3.25, 6, 100.0}");
}

int main(void) {
    audit_init();
    TEST(test_export_graphics_raster_no_crash);
    TEST(test_read_number_literal_only);
    return audit_finish("audit_file_io");
}
