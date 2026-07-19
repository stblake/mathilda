#ifndef FILES_H
#define FILES_H

#include "expr.h"

/* Filesystem predicates and path-string utilities.
 *
 *   FileExistsQ["name"]    True if a filesystem object named "name"
 *                          exists (file, directory, symlink, FIFO,
 *                          socket, device, ...).  False otherwise.
 *                          Interpreted relative to the current working
 *                          directory; $Path is not searched.
 *
 *   FileExtension["name"]  The substring after the last `.` in the
 *                          file-name component of "name", excluding
 *                          the dot.  "" when there is no extension,
 *                          when the name ends with `.`, when the name
 *                          is a pure directory (ends with `/`), or
 *                          when the leading dots are leading (".foo"
 *                          has no extension; "a.b.c" has extension
 *                          "c").  Pure string manipulation; no
 *                          filesystem access.
 *
 *   FileBaseName["name"]   The file-name component of "name" with
 *                          its final extension (if any) removed.
 *                          Mirrors FileExtension's notion of what
 *                          counts as an extension.  Pure string
 *                          manipulation; no filesystem access.
 *
 *   FilePrint["name"]      Prints the raw textual contents of "name"
 *   FilePrint["name", n]   to stdout.  With a positive integer second
 *   FilePrint["name", -n]  argument prints the first n lines; with a
 *   FilePrint["name",      negative integer prints the last |n|.  A
 *       m;;n]              Span argument selects a half-open line
 *   FilePrint["name",      range (1-based, inclusive) with optional
 *       m;;n;;s]           signed step `s`.  Returns Null on success
 *                          and $Failed if the file cannot be opened.
 *
 *   FileNameJoin[{...}]    Joins path components into a single file name
 *   FileNameJoin["name"]   suitable for the current operating system, or
 *                          canonicalizes a lone name.  An empty leading
 *                          component yields an absolute path.  Accepts an
 *                          OperatingSystem->"Windows"|"MacOSX"|"Unix"
 *                          option to select the separator convention.
 *                          Pure string manipulation; no filesystem access.
 */

void files_init(void);

Expr* builtin_fileexistsq(Expr* res);
Expr* builtin_fileextension(Expr* res);
Expr* builtin_filebasename(Expr* res);
Expr* builtin_fileprint(Expr* res);
Expr* builtin_filenamejoin(Expr* res);

#endif /* FILES_H */
