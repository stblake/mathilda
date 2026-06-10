#!/bin/sh
set -e

# Auto-initialize the ECM git submodule if the user cloned without
# --recurse-submodules.  Without this, src/external/ecm/ is an empty
# directory and `autoreconf -i` below fails with
# "'configure.ac' is required" (see GitHub issue #3).
if [ ! -f src/external/ecm/configure.ac ]; then
    if [ -d .git ] || [ -f .git ]; then
        echo "build_ecm.sh: ECM submodule not initialized, fetching..."
        git submodule update --init --recursive src/external/ecm
    else
        echo "build_ecm.sh: ERROR: src/external/ecm/configure.ac is missing" >&2
        echo "  and this is not a git checkout, so the submodule cannot be" >&2
        echo "  fetched automatically.  Re-clone with" >&2
        echo "    git clone --recurse-submodules https://github.com/stblake/Mathilda.git" >&2
        exit 1
    fi
fi

cd src/external/ecm

# Generate ECM's build system (configure, Makefile.in, libtool, ...) if it is
# not already present.  This is the ONLY step in the entire ECM build that
# needs Perl: autoreconf drives autoconf/automake/libtoolize, which are Perl
# programs.  Everything afterwards -- the generated `configure` script, libtool,
# and `make` -- is pure POSIX shell / m4 and needs no Perl.
#
# Therefore, if a pre-generated `configure` is already shipped (e.g. copied from
# a machine where this build has succeeded, or unpacked from a GMP-ECM release
# tarball), the whole ECM build runs Perl-free and this block is skipped.
if [ ! -f configure ]; then
    # `autoreconf` may be installed but non-functional -- the classic symptom is
    # a stale `#!/usr/bin/perlX.Y` shebang pointing at a Perl that no longer
    # exists, which makes the kernel reject it with "bad interpreter" (exit 126).
    # Probe it for real (not just `command -v`) so we can fail with an
    # actionable message instead of that cryptic error.
    if autoreconf --version >/dev/null 2>&1; then
        autoreconf -i
    else
        echo "build_ecm.sh: ERROR: ECM's 'configure' is missing and 'autoreconf'" >&2
        echo "  is unavailable or non-functional (it requires Perl + the GNU" >&2
        echo "  autotools, and one of those is broken on this machine)." >&2
        echo >&2
        echo "  Perl is only needed for this one-time bootstrap that generates" >&2
        echo "  'configure'.  The rest of the ECM build is Perl-free, so you have" >&2
        echo "  two options:" >&2
        echo >&2
        echo "    1. Copy a pre-generated build system into src/external/ecm/" >&2
        echo "       from a machine where 'make' has already succeeded -- at" >&2
        echo "       minimum 'configure', plus 'config.h.in', 'aclocal.m4', the" >&2
        echo "       'Makefile.in' files, the build-aux scripts (compile," >&2
        echo "       config.guess, config.sub, depcomp, install-sh, ltmain.sh," >&2
        echo "       missing, test-driver) and 'm4/lt*.m4'.  Then re-run 'make'." >&2
        echo >&2
        echo "    2. Install a working Perl plus autoconf, automake and libtool" >&2
        echo "       (e.g. 'brew install autoconf automake libtool', or your" >&2
        echo "       distro's package manager), then re-run 'make'." >&2
        echo >&2
        echo "  Alternatively, build Mathilda without ECM:  make USE_ECM=0" >&2
        exit 1
    fi
fi

# Fix endian functions on macOS
sed -i.bak 's/htole32/ /g' resume.c
sed -i.bak 's/htole64/ /g' resume.c
sed -i.bak 's/le32toh/ /g' resume.c
sed -i.bak 's/le64toh/ /g' resume.c

sed -i.bak 's/htole32/ /g' torsions.c
sed -i.bak 's/htole64/ /g' torsions.c
sed -i.bak 's/le32toh/ /g' torsions.c
sed -i.bak 's/le64toh/ /g' torsions.c

MY_CFLAGS="-O2"

# Force ECM's configure to skip its optional libprimesieve detection.
# configure.ac runs an unconditional AC_CHECK_LIB(primesieve,primesieve_init);
# on hosts where libprimesieve happens to be installed, libecm.a ends up
# with unresolved primesieve_* references that the Mathilda top-level link
# line doesn't satisfy (no -lprimesieve there), and the final link of
# Mathilda fails inside getprime_r.o.  ECM has its own internal prime
# iterator, so we pin the autoconf cache var to "no" to make the bundled
# build deterministic across hosts regardless of what's installed.
PRIMESIEVE_OFF="ac_cv_lib_primesieve_primesieve_init=no"

if [ -d /opt/homebrew/include ]; then
    env $PRIMESIEVE_OFF ./configure --with-gmp=/opt/homebrew CFLAGS="$MY_CFLAGS -I/opt/homebrew/include" LDFLAGS="-L/opt/homebrew/lib"
elif [ -d /usr/local/include ]; then
    env $PRIMESIEVE_OFF ./configure --with-gmp=/usr/local CFLAGS="$MY_CFLAGS -I/usr/local/include" LDFLAGS="-L/usr/local/lib"
else
    env $PRIMESIEVE_OFF ./configure CFLAGS="$MY_CFLAGS"
fi

make -j4
