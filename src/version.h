#ifndef MATHILDA_VERSION_H
#define MATHILDA_VERSION_H

/*
 * version.h — Mathilda release version.
 *
 * MATHILDA_VERSION_NUMBER is the single source of truth for the release and
 * backs the $VersionNumber system variable. MATHILDA_VERSION_STRING is its
 * textual form (kept in sync by hand) and is embedded in the $Version string.
 *
 * The full $Version string — Mathilda's version plus the versions of the
 * libraries it was compiled against — is assembled entirely at compile time
 * from preprocessor macros; see version.c. mathilda_version() returns it.
 */

#define MATHILDA_VERSION_NUMBER 0.017    /* keep in sync with the string below */
#define MATHILDA_VERSION_STRING "0.017"

/* Full descriptive version string, e.g.
 *   "Mathilda 0.01 (Apple LLVM 17.0.0, GMP 6.3.0, MPFR 4.2.2, FLINT 3.6.0, ...)"
 * The pointer refers to static storage and must not be freed. */
const char* mathilda_version(void);

#endif /* MATHILDA_VERSION_H */
