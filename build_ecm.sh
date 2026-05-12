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
        echo "    git clone --recurse-submodules https://github.com/stblake/picocas.git" >&2
        exit 1
    fi
fi

cd src/external/ecm

if [ ! -f configure ]; then
    autoreconf -i
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

if [ -d /opt/homebrew/include ]; then
    ./configure --with-gmp=/opt/homebrew CFLAGS="$MY_CFLAGS -I/opt/homebrew/include" LDFLAGS="-L/opt/homebrew/lib"
elif [ -d /usr/local/include ]; then
    ./configure --with-gmp=/usr/local CFLAGS="$MY_CFLAGS -I/usr/local/include" LDFLAGS="-L/usr/local/lib"
else
    ./configure CFLAGS="$MY_CFLAGS"
fi

make -j4
