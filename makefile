USE_ECM ?= 1
USE_MPFR ?= 1
USE_LAPACK ?= 1
USE_GRAPHICS ?= 1
USE_FLINT ?= 1
USE_REGEX ?= 1
USE_FFTW ?= 1

# Platform detection — used for readline and other OS-specific choices.
# On Windows under MSYS2/MinGW, uname returns "MINGW64_NT-*" or similar;
# map anything that isn't Darwin/Linux to "Windows" for build logic.
UNAME_S := $(shell uname -s 2>/dev/null || echo Windows)
ifneq ($(filter Darwin Linux,$(UNAME_S)),)
  BUILD_PLATFORM := $(UNAME_S)
else
  BUILD_PLATFORM := Windows
endif

# Compiler. Mathilda must be built with a REAL GCC, never Apple's clang shim:
# on macOS the plain `gcc` is a symlink to Apple clang, so a naive `CC = gcc`
# silently builds with LLVM (the $Version banner then reads "Apple LLVM ...").
# Auto-detect a genuine GCC by trying the newest Homebrew `gcc-NN` first and
# falling back to a plain `gcc` (which IS real GCC on most Linux distros). We
# only pick a name when CC was not set explicitly, so `make CC=clang` (or any
# other override on the command line / in the environment) is still honoured.
# A versioned name is fine here because it is *tried*, not hardcoded — hosts
# without it fall through to the next candidate.
ifeq ($(origin CC),default)
  CC := $(shell for c in gcc-16 gcc-15 gcc-14 gcc-13 gcc; do \
                  command -v $$c >/dev/null 2>&1 && { echo $$c; break; }; \
                done)
  # Warn loudly if the only `gcc` we found is really Apple clang — the build
  # still proceeds, but the operator should install a real GCC (brew install gcc).
  ifneq ($(shell $(CC) --version 2>/dev/null | grep -ci clang),0)
    $(warning Mathilda: '$(CC)' is Apple clang, not real GCC — install Homebrew GCC (brew install gcc) for a supported build.)
  endif
endif
# `-Werror=implicit-function-declaration`: under -std=c99 glibc hides every
# POSIX symbol, so reaching for one without a feature-test macro leaves the
# compiler guessing a signature — `jn` was assumed to return `int`, silently
# truncating Bessel values (issue #37). GCC 14 makes this an error on its own;
# older GCC and clang only warn, and a warning in a -j12 build scrolls past.
# Promote it everywhere so the failure is loud rather than a wrong answer.
# `make check-c99` catches the same class before the compiler ever sees it.
#
# `-Werror=incompatible-pointer-types` / `-int-conversion` / `-implicit-int`:
# the other diagnostics GCC 14 promoted to errors. Every one of them is a real
# type confusion, and each can be latent on macOS and fatal on Linux — issue #40
# was `ci_powi` handed an `int64_t*` where a `long long*` was expected, which is
# the SAME type on Darwin and a different one under glibc. The Linux CI job runs
# whatever GCC the runner ships (13 at the time of writing), where these are
# warnings that pass the build; promoting them here makes that job a real gate
# no matter which compiler version it lands on.
#
# `-Werror=unused-function`: a static function used only inside an `#ifdef
# USE_FLINT` / `USE_MPFR` block but defined unconditionally is dead code in the
# degrade config, where -Wall reports it as "defined but not used" / "declared
# static but never defined". Those warnings scroll past a -j build unnoticed on
# the dev machine (which has the optional deps) and only ever appear in the
# no-FLINT / no-MPFR CI job — the one place nobody reads the warning stream.
# Promote it so that job fails loudly instead: the fix is always to guard the
# definition with the same `#ifdef` as its caller.
CFLAGS = -O3 -std=c99 -Wall -Wextra -Werror=implicit-function-declaration \
         -Werror=incompatible-pointer-types -Werror=int-conversion \
         -Werror=implicit-int -Werror=unused-function -g -I./src -I./src/list -I./src/ml -I./src/linalg -I./src/numbertheory -I./src/poly -I./src/simp -I./src/stats -I./src/calculus -I./src/sum -I./src/product -I./src/special_functions -I./src/numerical_calculus -I./src/numerical_roots -I./src/graphics -I./src/graph -I./src/strings -I./src/strings/regex -I./src/ffi -I/usr/include -I/usr/local/include

# Readline is available on macOS and Linux but not on Windows (MinGW).
# Build with USE_READLINE=0 to disable it explicitly (e.g. for cross-builds
# or when only the pipe-mode sidecar is needed).
USE_READLINE ?= 1
ifeq ($(BUILD_PLATFORM),Windows)
  override USE_READLINE := 0
endif

ifeq ($(USE_READLINE),0)
  CFLAGS      += -DNO_READLINE
  READLINE_LIBS =
else
  READLINE_LIBS = -lreadline
endif

# POSIX threads accelerate the element-wise NDArray kernels (Erf, Sin, Exp, ...)
# by splitting large-array maps across cores. Available on macOS and Linux; the
# -pthread driver flag sets the right defines and links libpthread on both. Build
# with USE_THREADS=0 to force the serial path (e.g. thread-less platforms).
USE_THREADS ?= 1
ifeq ($(BUILD_PLATFORM),Windows)
  override USE_THREADS := 0
endif
ifeq ($(USE_THREADS),1)
  CFLAGS  += -DMATHILDA_THREADS -pthread
  THREAD_LIBS = -pthread
else
  THREAD_LIBS =
endif

LDFLAGS = $(READLINE_LIBS) $(THREAD_LIBS) -L/usr/local/lib -lgmp -lm

# Site-specific link libraries, appended verbatim to the link line. Some
# distributions need extra libraries the autodetection can't infer — e.g. on
# certain Ubuntu setups a statically-linked raylib pulls in `-lX11`, or a
# minimal LAPACKE package needs `-llapack` spelled out. Pass them on the
# command line rather than editing this file:
#   make EXTRA_LIBS="-llapack -lX11"
EXTRA_LIBS ?=

# Optional compile-time install prefix. When set, the kernel also looks for its
# bundled src/internal tree under $(PREFIX)/share/mathilda/internal, so a binary
# installed to $(PREFIX)/bin finds its modules with no MATHILDA_HOME needed:
#   make PREFIX=/usr/local && cp Mathilda /usr/local/bin
ifdef PREFIX
CFLAGS += -DMATHILDA_PREFIX=\"$(PREFIX)\"
endif

# GMP-ECM for advanced integer factorisation (facint.c: ecm_init/ecm_factor via
# the public ecm.h). System library only — no longer vendored as a submodule.
# GMP-ECM ships no pkg-config .pc file, so detection is a compile+link probe
# against the shared libecm (its transitive deps — primesieve, libomp — resolve
# automatically). When absent the build still succeeds with a runtime-degraded
# factoriser, matching the USE_FLINT=0 / USE_MPFR=0 graceful-degrade policy.
#   macOS (Homebrew): brew install gmp-ecm
#   Ubuntu/Debian:    sudo apt install libecm-dev
# NOTE: the probe MUST include <stdio.h> before <ecm.h>. GMP-ECM's ecm.h uses
# the FILE type (ecm_params has `FILE *os, *es;`) but on Debian/Ubuntu's 7.0.5
# it does not include <stdio.h> itself — it relies on the includer. macOS system
# headers pull in <stdio.h> transitively (masking the bug), Linux glibc does not,
# so without this the probe fails to COMPILE on Linux and ECM is wrongly reported
# "not detected" even when libecm-dev is installed. The multiarch lib dir is on
# the default linker path, so -lecm resolves without an explicit -L there.
ifeq ($(USE_ECM), 1)
  ECM_PROBE := $(shell printf '\#include <stdio.h>\n\#include <ecm.h>\nint main(void){ecm_params p;ecm_init(p);return 0;}\n' > /tmp/mathilda_ecmprobe.c 2>/dev/null && \
    $(CC) /tmp/mathilda_ecmprobe.c -o /tmp/mathilda_ecmprobe -I/usr/include -I/usr/local/include -I/opt/homebrew/include \
      -L/usr/local/lib -L/opt/homebrew/lib -lecm -lgmp 2>/dev/null && echo y; \
    rm -f /tmp/mathilda_ecmprobe.c /tmp/mathilda_ecmprobe)
  ifeq ($(ECM_PROBE), y)
    LDFLAGS += -lecm
  else
    $(warning GMP-ECM not detected; building with USE_ECM=0 (advanced factorisation disabled))
    $(warning   macOS (Homebrew): brew install gmp-ecm)
    $(warning   Ubuntu/Debian:    sudo apt install libecm-dev)
    override USE_ECM := 0
  endif
endif
ifneq ($(USE_ECM), 1)
CFLAGS += -DNO_ECM
endif

# Arbitrary-precision reals (MPFR) — enables N[expr, prec], Precision/
# Accuracy/SetPrecision/SetAccuracy, and `3.98`50` precision literals.
# Disable with `USE_MPFR=0` to build without the MPFR dependency; the
# machine-precision path continues to work with a runtime warning on the
# unsupported operations.
ifeq ($(USE_MPFR), 1)
CFLAGS  += -DUSE_MPFR
LDFLAGS += -lmpfr
endif

ifneq ($(wildcard /opt/homebrew/include),)
CFLAGS += -I/opt/homebrew/include
endif

ifneq ($(wildcard /opt/homebrew/lib),)
LDFLAGS += -L/opt/homebrew/lib
endif

# pkg-config wrapper. On macOS, Homebrew installs its .pc files outside the
# search path used by a MacPorts/system pkg-config, so `pkg-config --exists
# raylib` (below) fails on a stock shell even when raylib is installed. Prepend
# the Homebrew prefixes (Apple-Silicon + Intel) on Darwin so plain `make`
# detects raylib — and graphics is built by default — without a manual
# PKG_CONFIG_PATH override. Harmless on Linux (the directories simply don't
# exist) and respects any PKG_CONFIG_PATH the user already exported.
ifeq ($(BUILD_PLATFORM),Darwin)
  PKG_CONFIG = PKG_CONFIG_PATH="/opt/homebrew/lib/pkgconfig:/usr/local/lib/pkgconfig:$$PKG_CONFIG_PATH" pkg-config
else
  PKG_CONFIG = pkg-config
endif

# BLAS/LAPACK for fast machine-precision linear-algebra kernels
# (machine-precision QRDecomposition; later: Inverse, LinearSolve, Det,
# Eigenvalues, LeastSquares, SVD).  Four-tier autodetection:
#   1. Darwin              -> Apple Accelerate framework (zero install).
#   2. pkg-config lapacke  -> use pkg-config flags (OpenBLAS / MKL / etc).
#   3. /usr/include/lapacke.h or /usr/local/include/lapacke.h
#                          -> link -llapacke -llapack -lblas.
#   4. nothing found       -> warn and override USE_LAPACK := 0.
# When USE_LAPACK is off, machine-precision linalg falls back to the
# MPFR / symbolic path with a one-time runtime warning.  This matches
# the existing USE_MPFR=0 / USE_ECM=0 graceful-degrade policy so that
# `git clone && make` always succeeds, no matter the host environment.
ifeq ($(USE_LAPACK), 1)
  ifeq ($(BUILD_PLATFORM),Darwin)
    # Apple's vecLib/vBasicOps.h headers (pulled in by Accelerate.h) pass typed
    # SSE vectors (vUInt16, vUInt32, ...) into the compiler's _mm_* intrinsics,
    # which under GCC are declared with strict __m128i parameters. Clang permits
    # the implicit vector conversion; GCC needs -flax-vector-conversions or the
    # Accelerate include fails to compile. Harmless to our own code (it uses no
    # vector types) and scoped to the Darwin/Accelerate path only.
    CFLAGS  += -DUSE_LAPACK -DMATHILDA_USE_ACCELERATE -flax-vector-conversions
    LDFLAGS += -framework Accelerate
  else ifneq ($(shell pkg-config --exists lapacke 2>/dev/null && echo y),)
    CFLAGS  += -DUSE_LAPACK $(shell pkg-config --cflags lapacke)
    LDFLAGS += $(shell pkg-config --libs lapacke)
  else ifneq ($(wildcard /usr/include/lapacke.h)$(wildcard /usr/local/include/lapacke.h),)
    CFLAGS  += -DUSE_LAPACK
    LDFLAGS += -llapacke -llapack -lblas
  else
    $(warning LAPACK/LAPACKE not detected; building with USE_LAPACK=0)
    $(warning   Ubuntu/Debian:  sudo apt install liblapacke-dev libopenblas-dev)
    $(warning   Fedora/RHEL:    sudo dnf install lapack-devel lapacke-devel openblas-devel)
    $(warning   Arch:           sudo pacman -S openblas lapacke)
    override USE_LAPACK := 0
  endif
endif

# Raylib for the 2D graphics engine (Graphics[]/Show[]/Plot[] rendering).
# System library only (like GMP/Readline) — never vendored, since it's a
# full windowing/OpenGL framework. Autodetected via pkg-config; when absent
# the build still succeeds and Show/Plot fall back to a text placeholder
# at runtime, matching the USE_MPFR=0 / USE_LAPACK=0 graceful-degrade policy.
ifeq ($(USE_GRAPHICS), 1)
  ifneq ($(shell $(PKG_CONFIG) --exists raylib 2>/dev/null && echo y),)
    CFLAGS  += -DUSE_GRAPHICS $(shell $(PKG_CONFIG) --cflags raylib)
    # `--static` appends the `Libs.private:` transitive deps (on Linux: -lX11
    # -lGL -lpthread -ldl -lrt ...) that a static libraylib.a needs but `--libs`
    # alone omits — otherwise the link fails with undefined XInternAtom/XSync
    # (issue #18). Identical to `--libs` on macOS, where Libs.private is empty
    # (frameworks are recorded in the dylib). EXTRA_LIBS below remains as a
    # manual escape hatch for anything pkg-config still can't infer.
    LDFLAGS += $(shell $(PKG_CONFIG) --static --libs raylib)
  else
    $(warning Raylib not detected; building with USE_GRAPHICS=0 (Show/Plot will print a text placeholder))
    $(warning   macOS (Homebrew): brew install raylib)
    $(warning   Ubuntu/Debian:    sudo apt install libraylib-dev)
    override USE_GRAPHICS := 0
  endif
endif

# FLINT (>= 3.0) for fast, rigorous polynomial arithmetic over algebraic
# extensions — multivariate GCD/factoring over Q (fmpq_mpoly), univariate
# GCD/factoring over a number field Q(alpha) (the generic-ring `gr` layer +
# ANTIC, merged into FLINT at 3.0), and the finite-field multivariate workhorse
# (fq_nmod_mpoly) used by the parametric Q(t)(alpha) outer loop. See
# ALGEBRAIC_EXTENSION_ARITHMETIC_PLAN.md. Hard deps (GMP, MPFR) are already
# linked. System library only (LGPL-3, GPLv3-compatible) — never vendored.
# Autodetected via pkg-config with a >= 3.0 version floor (older packages lack
# ANTIC); when absent the build still succeeds and the algebraic-extension
# Cancel/Together/Apart/Factor paths fall back to the classical (slower but
# rigorous) code, matching the USE_MPFR=0 / USE_LAPACK=0 graceful-degrade policy.
ifeq ($(USE_FLINT), 1)
  ifneq ($(shell $(PKG_CONFIG) --exists 'flint >= 3.0' 2>/dev/null && echo y),)
    CFLAGS    += -DUSE_FLINT $(shell $(PKG_CONFIG) --cflags flint)
    FLINT_LIBS := $(shell $(PKG_CONFIG) --libs flint)
    # Debian/Ubuntu's flint.pc lists only the transitive `-lgmp -lmpfr` and
    # omits -lflint itself, so a pkg-config-only link fails with thousands of
    # undefined references to fmpz_*/arb_*/gr_* while the headers resolve fine.
    # Homebrew's .pc gets this right, which is why it was never seen locally.
    # Prepended, not appended: ld resolves left to right, so libflint has to
    # precede the -lgmp/-lmpfr it depends on.
    ifeq ($(filter -lflint,$(FLINT_LIBS)),)
      FLINT_LIBS := -lflint $(FLINT_LIBS)
    endif
    LDFLAGS   += $(FLINT_LIBS)
  else
    $(warning FLINT >= 3.0 not detected; building with USE_FLINT=0 (algebraic-extension GCD/Factor use the classical fallback))
    $(warning   macOS (Homebrew):      brew install flint)
    $(warning   Ubuntu 24.04+/Debian:  sudo apt install libflint-dev   (needs >= 3.0 for ANTIC))
    $(warning   Ubuntu 22.04 (apt has 2.x): ./tools/install-flint.sh   (builds 3.x from source; see docs/building.md))
    override USE_FLINT := 0
  endif
endif

# PCRE2 (Perl-Compatible Regular Expressions, 8-bit code units) backs
# RegularExpression[] and the regex-aware string functions (StringMatchQ,
# StringCases, StringReplace, StringSplit). This is the same engine the
# Wolfram Language uses, so RegularExpression syntax is faithful. System
# library only (BSD-licensed) — never vendored. Autodetected via pkg-config;
# when absent the build still succeeds and those builtins warn and stay
# unevaluated at runtime, matching the USE_MPFR=0 / USE_LAPACK=0 graceful-
# degrade policy.
ifeq ($(USE_REGEX), 1)
  ifneq ($(shell $(PKG_CONFIG) --exists libpcre2-8 2>/dev/null && echo y),)
    CFLAGS  += -DUSE_REGEX -DPCRE2_CODE_UNIT_WIDTH=8 $(shell $(PKG_CONFIG) --cflags libpcre2-8)
    LDFLAGS += $(shell $(PKG_CONFIG) --libs libpcre2-8)
  else
    $(warning PCRE2 not detected; building with USE_REGEX=0 (RegularExpression/StringMatchQ/StringCases/StringReplace/StringSplit warn and stay unevaluated))
    $(warning   macOS (Homebrew): brew install pcre2)
    $(warning   Ubuntu/Debian:    sudo apt install libpcre2-dev)
    override USE_REGEX := 0
  endif
endif

# FFTW (Fastest Fourier Transform in the West) backs the machine-precision path
# of Fourier[]/InverseFourier[] with an O(n log n) discrete Fourier transform.
# System library only (GPL) — never vendored. Only the double-precision `fftw3`
# module is used; the arbitrary-precision path is Mathilda's own MPFR-complex
# FFT. Autodetected via pkg-config; when absent the build still succeeds and the
# machine path falls back to a naive O(n^2) DFT, matching the USE_MPFR=0 /
# USE_LAPACK=0 graceful-degrade policy.
ifeq ($(USE_FFTW), 1)
  ifneq ($(shell $(PKG_CONFIG) --exists fftw3 2>/dev/null && echo y),)
    CFLAGS  += -DUSE_FFTW $(shell $(PKG_CONFIG) --cflags fftw3)
    LDFLAGS += $(shell $(PKG_CONFIG) --libs fftw3)
  else
    $(warning FFTW not detected; building with USE_FFTW=0 (Fourier uses a naive O(n^2) fallback))
    $(warning   macOS (Homebrew): brew install fftw)
    $(warning   Ubuntu/Debian:    sudo apt install libfftw3-dev)
    override USE_FFTW := 0
  endif
endif

SRC_DIR = src
SRC = $(wildcard $(SRC_DIR)/*.c) $(wildcard $(SRC_DIR)/list/*.c) $(wildcard $(SRC_DIR)/ml/*.c) $(wildcard $(SRC_DIR)/linalg/*.c) $(wildcard $(SRC_DIR)/numbertheory/*.c) $(wildcard $(SRC_DIR)/poly/*.c) $(wildcard $(SRC_DIR)/simp/*.c) $(wildcard $(SRC_DIR)/stats/*.c) $(wildcard $(SRC_DIR)/calculus/*.c) $(wildcard $(SRC_DIR)/sum/*.c) $(wildcard $(SRC_DIR)/product/*.c) $(wildcard $(SRC_DIR)/special_functions/*.c) $(wildcard $(SRC_DIR)/compile/*.c) $(wildcard $(SRC_DIR)/numerical_calculus/*.c) $(wildcard $(SRC_DIR)/numerical_roots/*.c) $(wildcard $(SRC_DIR)/graphics/*.c) $(wildcard $(SRC_DIR)/graph/*.c) $(wildcard $(SRC_DIR)/strings/*.c) $(wildcard $(SRC_DIR)/strings/regex/*.c)
ifneq ($(USE_GRAPHICS), 1)
SRC := $(filter-out $(SRC_DIR)/graphics/render.c $(SRC_DIR)/graphics/render3d.c $(SRC_DIR)/graphics/label_font.c, $(SRC))
endif
OBJ = $(SRC:.c=.o)

# FFI objects (src/ffi/*.c) are the embedding entry points. They belong in
# libmathilda.a REGARDLESS of USE_GRAPHICS, but NOT in the desktop Mathilda
# binary (which uses repl.c's main()); keep them out of SRC/OBJ and add them to
# the archive explicitly. (Tying them to USE_GRAPHICS=0 previously left the
# default host archive without any mathilda_ffi_* symbols.)
FFI_OBJ = $(patsubst %.c,%.o,$(wildcard $(SRC_DIR)/ffi/*.c))
# Per-object dependency files (.d), generated by the compiler's -MMD -MP. Each
# .d lists the headers its .o #includes, so editing a header (e.g. version.h)
# rebuilds every object that includes it. Pulled in via `-include` at the bottom.
DEP = $(OBJ:.o=.d)
TARGET = Mathilda

TEST_BINARIES = eval_tests expr_tests parse_tests test_ld test_ops test_pattern list_tests stats_tests expand_tests
TEST_DIR = tests
CMAKE_TEST_BINARIES = comparisons_tests eval_tests expr_tests match_tests match_extensive_tests parse_tests regression_tests symtab_tests list_tests trig_tests hyperbolic_tests logexp_tests piecewise_tests purefunc_tests stats_tests expand_tests numeric_tests numeric_largearg_tests numeric_stress_tests

all: $(TARGET)

# --- Link-flag de-duplication -------------------------------------------------
# pkg-config for FLINT (and GMP-ECM) re-lists -lgmp/-lmpfr that the base LDFLAGS
# already carry, so the raw link line repeats them and macOS ld warns "ignoring
# duplicate libraries". Collapse each -l flag to a single copy while keeping its
# LAST occurrence, so a provider (e.g. -lgmp) still follows every consumer
# (-lflint, -lmpfr, -lecm) — the right-to-left order a static archive link needs,
# and exactly why FLINT_LIBS is appended after the base -lgmp above. Only -l
# library words are de-duplicated: -L search paths, -pthread, and two-word
# "-framework X" tokens pass through untouched, ahead of the libraries.
ld_reverse    = $(if $(1),$(call ld_reverse,$(wordlist 2,$(words $(1)),$(1))) $(firstword $(1)))
ld_uniq_first = $(if $(1),$(firstword $(1)) $(call ld_uniq_first,$(filter-out $(firstword $(1)),$(1))))
ld_uniq_last  = $(call ld_reverse,$(call ld_uniq_first,$(call ld_reverse,$(1))))
LDFLAGS_DEDUP  = $(filter-out -l%,$(LDFLAGS)) $(call ld_uniq_last,$(filter -l%,$(LDFLAGS)))

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS_DEDUP) $(EXTRA_LIBS)

# ---------------------------------------------------------------------------
# Static library for embedding the kernel in-process (mobile hosts, FFI tests).
#
# iOS/Android sandboxes forbid spawning the `mathilda` sidecar, so the kernel is
# linked directly into the host app via src/ffi/mathilda_ffi.c. The archive is
# every kernel object EXCEPT repl.o (which carries main() + readline + the stdio
# pipe loop — none of which belong in an embedded library). The host links this
# .a plus -lgmp and provides its own entry point.
#
# The default host build works as-is. For a minimal, dependency-light kernel
# (what the iOS cross-build uses) pass the USE_* toggles off, e.g.:
#   make libmathilda.a USE_READLINE=0 USE_THREADS=0 USE_ECM=0 USE_MPFR=0 \
#        USE_LAPACK=0 USE_GRAPHICS=0 USE_FLINT=0 USE_REGEX=0 USE_FFTW=0
AR ?= ar
LIB_OBJ = $(filter-out $(SRC_DIR)/repl.o,$(OBJ)) $(FFI_OBJ)
# Remove any existing archive first: `ar rcs` MERGES into an existing .a rather
# than replacing it, so without this a cross-arch rebuild (iOS/Android) would
# leave stale host-arch members (render.o, label_font.o, …) behind, producing
# "neither ET_REL nor LLVM bitcode" linker warnings and a polluted archive.
libmathilda.a: $(LIB_OBJ)
	rm -f $@
	$(AR) rcs $@ $(LIB_OBJ)

# -MMD writes a .d file next to each .o listing the (non-system) headers it
# includes; -MP adds a phony target for each header so deleting a header does
# not break the build with a "no rule to make target" error. The .d files are
# consumed by the `-include $(DEP)` line below.
$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -f $(OBJ) $(DEP) $(TARGET)
	rm -rf *.dSYM
	# Sweep the whole src tree, not just the current source list, so objects and
	# dep files orphaned by a moved/renamed/deleted source (e.g. src/foo.c ->
	# src/bar/foo.c) don't linger as dead, unlinked artifacts.
	find $(SRC_DIR) \( -name '*.o' -o -name '*.d' \) -delete
	rm -f $(TEST_BINARIES)
	rm -f *~
	rm -f $(SRC_DIR)/*~
	rm -f $(TEST_DIR)/*~
	rm -f *.o
	if [ -f $(TEST_DIR)/Makefile ]; then $(MAKE) -C $(TEST_DIR) clean; fi
	rm -f $(addprefix $(TEST_DIR)/, $(CMAKE_TEST_BINARIES))

# Regenerate the documentation website's per-builtin pages from the docstrings,
# attributes, and spec examples. Requires the built ./Mathilda binary (examples
# are verified against it). The generated Markdown is committed; CI only builds
# the MkDocs site from it.
docs: $(TARGET)
	python3 site/generate.py

# Build the static site locally (needs `pip install -r site/requirements.txt`).
docs-build:
	mkdocs build --strict -f site/mkdocs.yml

# Serve the site with live reload at http://127.0.0.1:8000
docs-serve:
	mkdocs serve -f site/mkdocs.yml

# Portability gate: catch POSIX-only symbols that glibc hides under -std=c99 —
# both <math.h> constants (M_PI, M_E, ...) used without a C99 fallback and
# POSIX functions (jn, yn, strdup, fileno, ...) used without a feature-test
# macro. macOS exposes all of them implicitly, so these break only on Linux.
# Not part of `make all` (it needs python3); CI runs it on every push, next to
# a real Linux build — see .github/workflows/build.yml.
check-c99:
	python3 tools/check_c99_portability.py

# `make check-interval` — randomised stress test of Interval[] arithmetic against
# the inclusion (containment) guarantee, plus exactness and determinism. Needs a
# built ./Mathilda and python3, so it is not part of `all` (same status as
# check-c99). See tools/interval_fuzz.py.
check-interval:
	python3 tools/interval_fuzz.py

# `make check-packed-aware` — does every head with an NDArray fast path opt in
# to it? The packing gate materialises for any head NOT on src/pack.c's AWARE
# list, so a missing opt-in is correct, silent, and 30x-658x slow. That has
# happened four times; this reads the dispatch sites out of the source and
# diffs them against the list. Same "needs python3, so not part of `all`"
# status as check-c99.
# `make check-tests` — build and run EVERY test binary, once, in a pinned configuration.
#
# Until this existed, "the tests pass" was a claim no tool backed: CI compiled the tree and ran two
# source-level gates, and running the 400+ binaries was left to whoever remembered — so a PR once
# reported them as passing on the strength of 40 of 426. The script removes the choice of subset and
# pins the configuration, because a run in a tree left at USE_FLINT=OFF reported 36 failures of
# which 32 were the missing library. tests/known_failures.txt carries the standing failures with a
# reason each, so a NEW failure is loud and the gate is not red on arrival.
check-tests:
	tools/run_test_suite.sh

check-packed-aware:
	python3 tools/check_packed_aware.py

# `make check-array-exactness` — does any routine hand back a TWO-HEADED array
# from a machine input? A routine given a packed array must answer with a scalar
# or an array of one element head; an exact 0 invented inside a machine-real
# result is both wrong against Mathematica's numeric tower and unpackable, so
# every consumer downstream falls off the fast path. Six heads were doing it.
#
# Unlike the two checks above this one RUNS the binary (342 probes, ~5 min), so
# it is a release/pre-merge gate rather than a per-push one. Needs ./Mathilda
# built.
check-array-exactness:
	python3 tools/check_array_exactness.py

# `make check-image-packing` — does every image head hand back a PACKED buffer,
# at BOTH ranks? Three times an image operation was 4x to 23x slower than its
# equivalent elsewhere with entirely correct answers, because the marshalling and
# not the algorithm was the cost: image_load walking an NDArray element by
# element, image3d_load still walking after image_load was fixed, and
# bit_image_from_mask building 262144 Expr integers in nested Lists. No test in
# the suite could catch any of them; a benchmark caught each one by accident, one
# at a time. This asks the question mechanically instead.
check-image-packing:
	python3 tools/check_image_packing.py

# `make check-menu-ids` — does every native menu item actually do something?
#
# The menu bar is built in Rust and handled in TypeScript, joined only by a string id travelling
# through a `menu:<id>` event, and BOTH failure directions are silent: an id with no handler is a
# dead command, and a handler nothing emits is dead code that reads like wiring. Both were present
# the first time this was checked by hand.
check-menu-ids:
	python3 tools/check_menu_ids.py

# `make check-nd-surfaces` — does every head reach the SAME fast path from a
# packed List and from a visible NDArray, and agree on the answer?
#
# The two checks above are static: one reads dispatch sites, the other reads the
# registries. Neither can see a head that has no fast path on EITHER surface
# (which is how DeleteDuplicates stayed at 72x NumPy through four sweeps), and
# neither can see the two surfaces DISAGREEING — which they did, silently, until
# 2026-08-01: every real kernel truncated a visible int64 NDArray, so
# Sin[NDArray[{1,2,3}, DataType -> "int64"]] was {0, 0, 0}. The gate that keeps
# the packed surface safe is exactly what leaves the visible one unguarded.
#
# This runs the binary over numeric_sweep.py's 284 probes on all three
# representations, so it is a release/pre-merge gate rather than a per-push one.
# `--survival` is the cheap half and answers a different question: which
# producers hand back a plain List that would have packed, making their
# CONSUMERS slow. Needs ./Mathilda built.
check-nd-surfaces:
	python3 tools/nd_surface_audit.py --survival

# `make check-compile-coverage` — does every head with a numeric fast path also
# COMPILE?
#
# A head earns a numeric fast path by being something a numeric workload runs
# over machine numbers, and Compile[] exists to run those workloads: a head fast
# at the REPL and unlowerable inside Compile[] is a contradiction. It is not
# merely a missed speedup, because the compilable subset is a CLIFF — one
# unlowered head sends the WHOLE body to the interpreter, silently, so
# Compile[{{v,_Real,1}}, Total[v]/Mean[v]] used to lose the compiled Total too.
#
# Joins the kernel registry and pack.c's AWARE list against the binary's own
# CompileDiagnostics over every typed argument shape. Exits non-zero while any
# head is unlowered without a recorded exemption, so the gap list is a work
# queue that cannot quietly grow. Needs ./Mathilda built.
check-compile-coverage:
	python3 tools/compile_coverage.py

# `make check-fastpath-sweep` — is every head that CAN take a machine array
# actually reaching the buffer, measured rather than declared?
#
# The three audits above all either read the source or read a curated list of
# expressions, and share one blind spot: a head nobody thought to name. All
# three were green on the day Commonest was found costing 880 ms where Tally of
# the identical 10^7 buffer cost 21.5 ms. This names nobody — it enumerates
# every registered builtin, discovers by TRIAL which call shapes each accepts,
# and times the survivors packed and under MATHILDA_NO_PACK=1. A head that is
# both expensive per element AND indifferent to whether its input was packed is
# not on the buffer.
#
# The gate pass is the DETECTOR and is what this target runs: it counts, per
# head, how many elements the transparency gate materialised — a count, not a
# duration, so it is deterministic and immune to load. Drop --gate-only to add
# the timing half, which ranks severity and must have an idle machine (a timing
# tool sharing a machine with anything else is measuring the other thing too).
#
# Ratchets against a checked-in OFF_BUFFER list, so it fails on a head that has
# NEWLY fallen off the buffer and merely reports the standing backlog. A gate
# that fails on everything from the day it lands is noise within a week.
#
# Minutes, not seconds: a release gate, not a per-push one. Needs ./Mathilda.
check-fastpath-sweep:
	python3 tools/nd_fastpath_sweep.py --gate-only

# `make bench-gap` — the weekly gap-driven benchmark job.
#
# Runs all 31 `benchmarks/NN-slug/` experiments in Mathilda, in Python
# (numpy/scipy/sympy/networkx) and — when wolframscript is installed — in
# Mathematica, joins the rows by label, and writes a ranked report naming the
# week's work: benchmarks/REPORT.md, benchmarks/ABSENT.md, and the raw rows to
# benchmarks/results/<date>.json.
#
# Distinct from the check-* gates above: those answer "did anything regress",
# this answers "where is Mathilda behind, and is it behind because a kernel is
# slow or because a function does not exist". Those two are reported separately
# and never pooled — a `SLOWER` row carries a ratio, an `ABSENT` row never does.
#
# Minutes, not seconds. Needs ./Mathilda; run `make` first.
#
#   make bench-gap                                   # everything
#   python3 benchmarks/run_all.py --only 31           # one experiment
#   python3 benchmarks/run_all.py --system mathilda,python
bench-gap:
	python3 benchmarks/run_all.py

# `make check-diophantine-heldout` — the held-out Solve[..., Integers] gate.
#
# Runs Mathilda COLD on equations drawn from standard references (NOT the
# co-designed benchmarks/87 cases.py) and cross-checks every answer against an
# independent Python brute-force oracle over the same box. Fails (nonzero exit)
# on any SILENT WRONG ANSWER -- a {}/finite/parametric result the oracle
# contradicts, the one class the developed-against benchmark cannot see. Needs
# only the Mathilda binary (no sympy); writes HELDOUT_REPORT.md.
check-diophantine-heldout:
	python3 benchmarks/87-diophantine-integers/validate.py

# Report the compiler the build will ACTUALLY use. `gcc --version` does not
# answer that: the autodetection above prefers a versioned `gcc-NN` over the
# plain name, so on a host with both, a bare `gcc --version` names one compiler
# while every compile line names another. The Linux CI job printed 13.3.0 that
# way while building with gcc-14, which is how issue #40's diagnosis started
# off pointing at the wrong toolchain.
print-cc:
	@echo "CC = $(CC)"
	@$(CC) --version 2>/dev/null | head -1

.PHONY: all clean docs docs-build docs-serve check-c99 check-interval check-packed-aware \
        check-array-exactness check-nd-surfaces check-compile-coverage \
        check-fastpath-sweep check-menu-ids check-tests bench-gap \
        check-diophantine-heldout print-cc

# Pull in the auto-generated header dependencies. The leading `-` silences the
# "no such file" notice on a fresh tree (no .d files exist until the first
# compile). Placed last so it never becomes the default goal.
-include $(DEP)
