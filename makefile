USE_ECM ?= 1
USE_MPFR ?= 1
USE_LAPACK ?= 1
USE_GRAPHICS ?= 1
USE_FLINT ?= 1
USE_REGEX ?= 1

# Platform detection — used for readline and other OS-specific choices.
# On Windows under MSYS2/MinGW, uname returns "MINGW64_NT-*" or similar;
# map anything that isn't Darwin/Linux to "Windows" for build logic.
UNAME_S := $(shell uname -s 2>/dev/null || echo Windows)
ifneq ($(filter Darwin Linux,$(UNAME_S)),)
  BUILD_PLATFORM := $(UNAME_S)
else
  BUILD_PLATFORM := Windows
endif

CC = gcc
CFLAGS = -O3 -std=c99 -Wall -Wextra -g -I./src -I./src/list -I./src/linalg -I./src/numbertheory -I./src/poly -I./src/simp -I./src/calculus -I./src/sum -I./src/product -I./src/special_functions -I./src/numerical_calculus -I./src/numerical_roots -I./src/graphics -I./src/strings -I./src/strings/regex -I/usr/local/include

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

LDFLAGS = $(READLINE_LIBS) -L/usr/local/lib -lgmp -lm

# GMP-ECM for advanced integer factorisation (facint.c: ecm_init/ecm_factor via
# the public ecm.h). System library only — no longer vendored as a submodule.
# GMP-ECM ships no pkg-config .pc file, so detection is a compile+link probe
# against the shared libecm (its transitive deps — primesieve, libomp — resolve
# automatically). When absent the build still succeeds with a runtime-degraded
# factoriser, matching the USE_FLINT=0 / USE_MPFR=0 graceful-degrade policy.
#   macOS (Homebrew): brew install gmp-ecm
#   Ubuntu/Debian:    sudo apt install libecm-dev
ifeq ($(USE_ECM), 1)
  ECM_PROBE := $(shell printf '\#include <ecm.h>\nint main(void){ecm_params p;ecm_init(p);return 0;}\n' > /tmp/mathilda_ecmprobe.c 2>/dev/null && \
    $(CC) /tmp/mathilda_ecmprobe.c -o /tmp/mathilda_ecmprobe -I/usr/local/include -I/opt/homebrew/include \
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
    CFLAGS  += -DUSE_LAPACK -DMATHILDA_USE_ACCELERATE
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
    LDFLAGS += $(shell $(PKG_CONFIG) --libs raylib)
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

SRC_DIR = src
SRC = $(wildcard $(SRC_DIR)/*.c) $(wildcard $(SRC_DIR)/list/*.c) $(wildcard $(SRC_DIR)/linalg/*.c) $(wildcard $(SRC_DIR)/numbertheory/*.c) $(wildcard $(SRC_DIR)/poly/*.c) $(wildcard $(SRC_DIR)/simp/*.c) $(wildcard $(SRC_DIR)/calculus/*.c) $(wildcard $(SRC_DIR)/sum/*.c) $(wildcard $(SRC_DIR)/product/*.c) $(wildcard $(SRC_DIR)/special_functions/*.c) $(wildcard $(SRC_DIR)/numerical_calculus/*.c) $(wildcard $(SRC_DIR)/numerical_roots/*.c) $(wildcard $(SRC_DIR)/graphics/*.c) $(wildcard $(SRC_DIR)/strings/*.c) $(wildcard $(SRC_DIR)/strings/regex/*.c)
ifneq ($(USE_GRAPHICS), 1)
SRC := $(filter-out $(SRC_DIR)/graphics/render.c $(SRC_DIR)/graphics/render3d.c $(SRC_DIR)/graphics/hershey_font.c, $(SRC))
endif
OBJ = $(SRC:.c=.o)
TARGET = Mathilda

TEST_BINARIES = eval_tests expr_tests parse_tests test_ld test_ops test_pattern list_tests stats_tests expand_tests
TEST_DIR = tests
CMAKE_TEST_BINARIES = comparisons_tests eval_tests expr_tests match_tests match_extensive_tests parse_tests regression_tests symtab_tests list_tests trig_tests hyperbolic_tests logexp_tests piecewise_tests purefunc_tests stats_tests expand_tests

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(SRC_DIR)/boolean.o: $(SRC_DIR)/boolean.c
	$(CC) $(CFLAGS) -c $< -o $@

# hershey_font.c includes the generated glyph table; the pattern rule above only
# tracks the .c, so name the .inc as an extra prerequisite to force a rebuild
# when the font data is regenerated (tools/gen_hershey.py).
$(SRC_DIR)/graphics/hershey_font.o: $(SRC_DIR)/graphics/hershey_glyphs.inc

clean:
	rm -f $(OBJ) $(TARGET)
	rm -rf *.dSYM
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
