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
CFLAGS = -O3 -std=c99 -Wall -Wextra -g -I./src -I./src/list -I./src/linalg -I./src/numbertheory -I./src/poly -I./src/simp -I./src/calculus -I./src/sum -I./src/product -I./src/special_functions -I./src/numerical_calculus -I./src/numerical_roots -I./src/graphics -I./src/graph -I./src/strings -I./src/strings/regex -I./src/ffi -I/usr/include -I/usr/local/include

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
    CFLAGS  += -DUSE_FLINT $(shell $(PKG_CONFIG) --cflags flint)
    LDFLAGS += $(shell $(PKG_CONFIG) --libs flint)
  else
    $(warning FLINT >= 3.0 not detected; building with USE_FLINT=0 (algebraic-extension GCD/Factor use the classical fallback))
    $(warning   macOS (Homebrew): brew install flint)
    $(warning   Ubuntu/Debian:    sudo apt install libflint-dev   (needs >= 3.0 for ANTIC))
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
SRC = $(wildcard $(SRC_DIR)/*.c) $(wildcard $(SRC_DIR)/list/*.c) $(wildcard $(SRC_DIR)/linalg/*.c) $(wildcard $(SRC_DIR)/numbertheory/*.c) $(wildcard $(SRC_DIR)/poly/*.c) $(wildcard $(SRC_DIR)/simp/*.c) $(wildcard $(SRC_DIR)/calculus/*.c) $(wildcard $(SRC_DIR)/sum/*.c) $(wildcard $(SRC_DIR)/product/*.c) $(wildcard $(SRC_DIR)/special_functions/*.c) $(wildcard $(SRC_DIR)/numerical_calculus/*.c) $(wildcard $(SRC_DIR)/numerical_roots/*.c) $(wildcard $(SRC_DIR)/graphics/*.c) $(wildcard $(SRC_DIR)/graph/*.c) $(wildcard $(SRC_DIR)/strings/*.c) $(wildcard $(SRC_DIR)/strings/regex/*.c)
ifneq ($(USE_GRAPHICS), 1)
SRC := $(filter-out $(SRC_DIR)/graphics/render.c $(SRC_DIR)/graphics/render3d.c $(SRC_DIR)/graphics/hershey_font.c, $(SRC))
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
CMAKE_TEST_BINARIES = comparisons_tests eval_tests expr_tests match_tests match_extensive_tests parse_tests regression_tests symtab_tests list_tests trig_tests hyperbolic_tests logexp_tests piecewise_tests purefunc_tests stats_tests expand_tests

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(EXTRA_LIBS)

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
# leave stale host-arch members (render.o, hershey_font.o, …) behind, producing
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

# hershey_font.c includes the generated glyph table. Header-dependency tracking
# (-MMD, above) picks this up automatically once a .d exists, but the .d is only
# written after the first successful compile — name the .inc explicitly so the
# very first build (and a `make clean` build) also rebuilds when the font data
# is regenerated (tools/gen_hershey.py).
$(SRC_DIR)/graphics/hershey_font.o: $(SRC_DIR)/graphics/hershey_glyphs.inc

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

.PHONY: all clean docs docs-build docs-serve

# Pull in the auto-generated header dependencies. The leading `-` silences the
# "no such file" notice on a fresh tree (no .d files exist until the first
# compile). Placed last so it never becomes the default goal.
-include $(DEP)
