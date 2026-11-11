/*
 * lean_run.h -- running the Lean 4 checker from Mathilda.
 *
 * The process layer is ISO C system(): no fork / popen (POSIX <unistd.h> /
 * <fcntl.h> are used only for the private temp directory: mkdtemp, open with
 * O_EXCL, rmdir).  The
 * command line is assembled from an argv whose every element is quoted for
 * the host shell, and holds only executable paths, flags written by Mathilda
 * and temp-file paths: the Lean source, tactic strings and claims are always
 * written to a .lean FILE, never to the command line.
 *
 *   POSIX:  set -C; [cd 'CWD' &&] 'ARG0' 'ARG1' ... > 'OUT' 2> 'ERR'
 *           Single quotes make every character literal (an embedded ' is
 *           written '\''), so the quoting is complete; OUT and ERR lie in
 *           a fresh mkdtemp directory (mode 0700), and set -C (noclobber)
 *           makes the redirections fail on an existing file anyway.  The
 *           exit status is decoded with WEXITSTATUS.
 *   _WIN32: [cd /d "CWD" &&] "ARG0" "ARG1" ... > "OUT" 2> "ERR"; an argument
 *           containing any of " % ! ^ & | < > is refused (cmd.exe has no
 *           complete quoting).  This branch compiles under MinGW but is not
 *           exercised at run time.
 *
 * system() blocks.  TimeConstrained's profiling timer measures this
 * process's CPU, which a parent waiting in system() does not use, so a
 * timeout never fires mid-run and no child is orphaned by a longjmp; the
 * resource bound on a Lean run is the deterministic maxHeartbeats / -T set
 * by the generator.
 */
#ifndef MATHILDA_LEAN_RUN_H
#define MATHILDA_LEAN_RUN_H

#include "expr.h"

typedef struct {
    int exit_code;      /* the process's exit status (127: command not found) */
    char* out;          /* captured stdout (owned, never NULL after a run) */
    char* err;          /* captured stderr (owned, never NULL after a run) */
    double seconds;     /* wall time of the run */
} LeanProcResult;

/* Run argv[0..argc-1] (quoted for the host shell) in `cwd` (NULL: the current
 * directory), capturing stdout / stderr through temp files.  Returns 0 when the
 * process ran (r filled; free with lean_proc_result_free), -1 when it could not
 * be started (an argument cannot be quoted safely, or the temp files could not
 * be created); *why then names the reason. */
int lean_run_process(const char* const* argv, int argc, const char* cwd, LeanProcResult* r,
                     const char** why);
void lean_proc_result_free(LeanProcResult* r);

/* Registers LeanForm's strict translator (Lean`Internal`ToLean) and the process
 * builtins Lean`Internal`RunProcess and Lean`Internal`VerifySource.  Called from core_init(). */
void lean_init(void);

#endif /* MATHILDA_LEAN_RUN_H */
