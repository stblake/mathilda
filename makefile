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

# Homebrew prefix (issue #90). Mathilda's optional dependencies (readline,
# gmp-ecm, flint, raylib, pcre2, fftw, ...) are typically installed via Homebrew
# on macOS, whose keg layout sits outside the default compiler / pkg-config
# search paths; the probes below add the prefix's include/lib/pkgconfig dirs so a
# plain `make` finds them. Rather than hardcode /opt/homebrew (the Apple-Silicon
# default, wrong on Intel and on a custom prefix), resolve HOMEBREW once:
#
#   1. An explicit HOMEBREW=... (command line or environment) always wins. This
#      is also the documented opt-out: HOMEBREW= (empty) or HOMEBREW=/nonexisting
#      makes every $(wildcard $(HOMEBREW)/...) guard below miss, disabling all
#      Homebrew discovery.
#   2. Otherwise `brew --prefix` if a brew is on PATH — this handles the
#      Apple-Silicon /opt/homebrew, the Intel /usr/local, a Linuxbrew prefix, and
#      any custom install uniformly.
#   3. Otherwise, on macOS only, the first conventional prefix that exists.
#
# Only the explicit-override case ($(origin) test) skips auto-detection, so the
# two opt-out spellings in (1) are honoured verbatim.
ifeq ($(origin HOMEBREW),undefined)
  HOMEBREW := $(shell brew --prefix 2>/dev/null)
  ifeq ($(strip $(HOMEBREW)),)
    ifeq ($(BUILD_PLATFORM),Darwin)
      HOMEBREW := $(firstword $(wildcard /opt/homebrew /usr/local))
    endif
  endif
endif

# Homebrew-derived search paths, each added only if it actually exists under the
# resolved prefix. All empty when HOMEBREW is empty/absent (opt-out, Linux
# without Linuxbrew, Windows), so the consumers below collapse to the
# no-Homebrew behaviour.
HOMEBREW_CFLAGS  :=
HOMEBREW_LDFLAGS :=
HOMEBREW_PC_PATH :=
ifneq ($(strip $(HOMEBREW)),)
  ifneq ($(wildcard $(HOMEBREW)/include),)
    HOMEBREW_CFLAGS := -I$(HOMEBREW)/include
  endif
  ifneq ($(wildcard $(HOMEBREW)/lib),)
    HOMEBREW_LDFLAGS := -L$(HOMEBREW)/lib
  endif
  ifneq ($(wildcard $(HOMEBREW)/lib/pkgconfig),)
    HOMEBREW_PC_PATH := $(HOMEBREW)/lib/pkgconfig
  endif
endif

# Conventional Unix site prefix (/usr/local), searched for headers and libraries
# by the compile/link lines and the ECM/LAPACK probes below. Unix hosts only: on
# MSYS2 `/usr` is the MSYS (Cygwin) runtime root, so `-I/usr/include` or
# `-I/usr/local/include` there makes the native UCRT64/MINGW64 toolchain read
# Cygwin headers ahead of its own (passagemath PR #2966). The default is empty on
# Windows. An explicit SYSTEM_PREFIX=... always wins; SYSTEM_PREFIX= (empty) is the
# packager opt-out, the twin of HOMEBREW= above, so a build against a private
# prefix searches nothing it was not given. /usr/include itself is never added:
# on Linux/macOS it is a system directory GCC already searches (and ignores as an
# -I), and on MSYS2 it is the Cygwin tree.
ifeq ($(origin SYSTEM_PREFIX),undefined)
  ifeq ($(BUILD_PLATFORM),Windows)
    SYSTEM_PREFIX :=
  else
    SYSTEM_PREFIX := /usr/local
  endif
endif
SYSPREFIX_CFLAGS  :=
SYSPREFIX_LDFLAGS :=
ifneq ($(strip $(SYSTEM_PREFIX)),)
  ifneq ($(wildcard $(SYSTEM_PREFIX)/include),)
    SYSPREFIX_CFLAGS := -I$(SYSTEM_PREFIX)/include
  endif
  ifneq ($(wildcard $(SYSTEM_PREFIX)/lib),)
    SYSPREFIX_LDFLAGS := -L$(SYSTEM_PREFIX)/lib
  endif
endif

# Compiler. Mathilda must be built with a REAL GCC, never Apple's clang shim: on
# macOS the plain `gcc` is a symlink to Apple clang, and LLVM both miscompiles some
# of our -std=c99 paths and hides the glibc portability class the Linux CI gates on.
# "GCC only, never clang" is a HARD rule for this project — enforced below as an
# error, not a warning.
#
# Auto-detect a genuine GCC by trying the newest Homebrew `gcc-NN` first and falling
# back to a plain `gcc` (which IS real GCC on most Linux distros). We only pick a
# name when CC was not set explicitly, so `make CC=gcc-14` (or any other override on
# the command line / in the environment) is still honoured. A versioned name is fine
# here because it is *tried*, not hardcoded — hosts without it fall through.
# Platform-specific advice for the two compiler errors below. CLANG64/CLANGARM64
# on MSYS2 are deliberately unsupported: there is no GCC in those environments.
ifeq ($(BUILD_PLATFORM),Darwin)
  GCC_HINT := Install Homebrew GCC (brew install gcc) and use the default CC or CC=gcc-16
else ifeq ($(BUILD_PLATFORM),Linux)
  GCC_HINT := Install your distribution's gcc package and use the default CC or CC=gcc
else
  GCC_HINT := On Windows build in the MSYS2 UCRT64 (or MINGW64) environment with mingw-w64-ucrt-x86_64-gcc; the CLANG64/CLANGARM64 environments are not supported
endif
ifeq ($(origin CC),default)
  CC := $(shell for c in gcc-16 gcc-15 gcc-14 gcc-13 gcc; do \
                  command -v $$c >/dev/null 2>&1 && { echo $$c; break; }; \
                done)
  ifeq ($(strip $(CC)),)
    $(error Mathilda: no C compiler found. Mathilda builds with real GCC only. $(GCC_HINT))
  endif
endif
# Enforce the rule for BOTH the autodetected name and an explicit CC=...: if the
# selected compiler reports itself as clang, STOP rather than silently produce an
# unsupported binary. The fix is always `brew install gcc`, then build with the
# default CC (or CC=gcc-NN).
ifneq ($(shell $(CC) --version 2>/dev/null | grep -ci clang),0)
  $(error Mathilda: '$(CC)' is clang, not GCC — Mathilda must be built with real GCC. $(GCC_HINT))
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
         -Werror=implicit-int -Werror=unused-function -g -I./src -I./src/list -I./src/ml -I./src/linalg -I./src/numbertheory -I./src/numbertheory/factorinteger -I./src/poly -I./src/simp -I./src/stats -I./src/calculus -I./src/sum -I./src/product -I./src/special_functions -I./src/numerical_calculus -I./src/numerical_roots -I./src/graphics -I./src/graph -I./src/strings -I./src/strings/regex -I./src/solve -I./src/ffi -I./src/bitwise -I./src/sparse $(SYSPREFIX_CFLAGS)

# UpValues subsystem (trial). Every change is gated behind `#if UP_VALUES`, so
# UP_VALUES=0 compiles the tree back to the exact pre-UpValues behaviour. ON by
# default; build without it via `make UP_VALUES=0`.
UP_VALUES ?= 1
ifeq ($(UP_VALUES),1)
  CFLAGS += -DUP_VALUES=1
endif

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
  # macOS exposes a libedit *shim* as <readline/readline.h> (in the SDK) and
  # /usr/lib/libreadline -> libedit. That shim lacks GNU Readline entry points
  # (rl_bind_keyseq, ...) and, more importantly, cannot do the in-buffer
  # multi-line editing the REPL's completeness-driven Return relies on:
  # rebinding RET there is silently ignored. Prefer a real Homebrew GNU Readline
  # when one is installed. It is keg-only, so its include/lib dirs are NOT on the
  # default search paths and must be added explicitly; when none is found the
  # build falls back to the shim (single-line editing only). Linux resolves GNU
  # Readline through -lreadline already, so this probe is Darwin-only.
  ifeq ($(BUILD_PLATFORM),Darwin)
    # Keg-only: Homebrew symlinks readline at $(HOMEBREW)/opt/readline, so the
    # resolved prefix locates it without a second `brew` call. An empty /
    # non-existent HOMEBREW (the opt-out) finds nothing and falls back to the
    # libedit shim (single-line editing only).
    READLINE_PREFIX := $(patsubst %/include/readline/readline.h,%,\
                         $(wildcard $(HOMEBREW)/opt/readline/include/readline/readline.h))
    ifneq ($(READLINE_PREFIX),)
      CFLAGS        += -I$(READLINE_PREFIX)/include
      READLINE_LIBS := -L$(READLINE_PREFIX)/lib -lreadline
    endif
  endif
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

LDFLAGS = $(READLINE_LIBS) $(THREAD_LIBS) $(SYSPREFIX_LDFLAGS) -lgmp -lm

# Embeddable codegen (issue #88). Every object is compiled position-independent
# and with hidden default visibility, so ONE object set serves the Mathilda
# executable, libmathilda.a, and the shared libmathilda:
#   -fPIC               a Cython/Python extension is a shared object, and on
#                       Linux x86-64 non-PIC (or PIE) objects from libmathilda.a
#                       cannot be linked into one ("recompile with -fPIC").
#   -fvisibility=hidden only the MATHILDA_FFI_API entry points (mathilda_ffi.h)
#                       are dynamically exported; the ~10k internal names
#                       (evaluate, symtab_init, ...) stay out of a host process.
#                       Hidden symbols still resolve in any STATIC link (the
#                       executable, xmathilda, the mobile staticlibs), and hidden
#                       + PIC binds internal calls directly (no GOT/PLT), so the
#                       executable costs the same as the PIE build it replaces.
# Windows (PE) has neither notion — exports are driven by dllexport there.
ifneq ($(BUILD_PLATFORM),Windows)
  CFLAGS += -fPIC -fvisibility=hidden
endif
# dladdr (src/loadmodule.c: library-relative module lookup) lives in libdl on
# glibc < 2.34; on newer glibc libdl is an empty stub, so linking it is harmless.
ifeq ($(BUILD_PLATFORM),Linux)
  LDFLAGS += -ldl
endif

# Native Windows (MinGW): MemoryInUse[] / MaxMemoryUsed[] use GetProcessMemoryInfo
# from psapi.  kernel32 (QueryPerformanceCounter / GetSystemTimePreciseAsFileTime
# / Sleep, used by src/mtime.c) is linked by the MinGW driver automatically;
# psapi is not, so add it explicitly.
ifeq ($(BUILD_PLATFORM),Windows)
  LDFLAGS += -lpsapi
endif

# macOS: reserve maximum Mach-O header padding so the baked-in absolute dylib
# install names (/usr/local/opt/.../libgmp.10.dylib, …) can later be rewritten
# to LONGER relocatable paths (@loader_path/../../.dylibs/…) by install_name_tool
# WITHOUT a relink. Downstream packagers that repackage the binary into a
# relocatable bundle — Homebrew, conda-forge, and the passagemath pip wheel
# (delocate) — all do exactly that rewrite, and the longer load commands do not
# fit in the default padding the linker happens to leave, failing with
# "larger updated load commands do not fit ... use -headerpad_max_install_names".
# CMake adds this automatically on Apple platforms; a hand-written makefile must
# opt in. No-op at build time; Darwin-only (the flag is an Apple ld feature).
ifeq ($(BUILD_PLATFORM),Darwin)
LDFLAGS += -Wl,-headerpad_max_install_names
endif

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
# The '#' of each #include is written as printf's octal escape \043, NOT as
# make's \#: GNU Make 3.81 strips the backslash of \# inside a function call
# but Make >= 4.3 keeps it, so the probe source began `\#include` and failed to
# compile on every current Linux distro (issue #93). Make leaves \043 alone and
# printf turns it into '#' under any make. The $$$$ (the shell's PID) keeps
# concurrent builds from clobbering each other's probe files.
ifeq ($(USE_ECM), 1)
  ECM_PROBE := $(shell f=/tmp/mathilda_ecmprobe_$$$$; \
    printf '\043include <stdio.h>\n\043include <ecm.h>\nint main(void){ecm_params p;ecm_init(p);return 0;}\n' > $$f.c 2>/dev/null && \
    $(CC) $$f.c -o $$f $(SYSPREFIX_CFLAGS) $(HOMEBREW_CFLAGS) \
      $(SYSPREFIX_LDFLAGS) $(HOMEBREW_LDFLAGS) -lecm -lgmp 2>/dev/null && echo y; \
    rm -f $$f.c $$f)
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
# unsupported operations.  MPFR is default-ON and has no autodetect/degrade
# (unlike FLINT/LAPACK/...), so on native Windows install the package first
# (`pacman -S mingw-w64-ucrt-x86_64-mpfr`) or build with USE_MPFR=0 -- otherwise
# the link fails on -lmpfr.  See docs/building.md (§Windows).
ifeq ($(USE_MPFR), 1)
CFLAGS  += -DUSE_MPFR
LDFLAGS += -lmpfr
endif

# Homebrew include/lib dirs (resolved above; empty when HOMEBREW is opted out or
# absent). SYSTEM_PREFIX (/usr/local on Unix) is already on CFLAGS/LDFLAGS, so
# on an Intel-prefix install these may duplicate it harmlessly.
CFLAGS  += $(HOMEBREW_CFLAGS)
LDFLAGS += $(HOMEBREW_LDFLAGS)

# pkg-config wrapper. On macOS, Homebrew installs its .pc files outside the
# search path used by a MacPorts/system pkg-config, so `pkg-config --exists
# raylib` (below) fails on a stock shell even when raylib is installed. Prepend
# the resolved Homebrew prefix's pkgconfig dir (HOMEBREW_PC_PATH, computed above)
# so plain `make` detects raylib — and graphics is built by default — without a
# manual PKG_CONFIG_PATH override. Empty when Homebrew is opted out or absent
# (Linux, Windows), so this collapses to plain pkg-config there; respects any
# PKG_CONFIG_PATH the user already exported.
ifneq ($(HOMEBREW_PC_PATH),)
  PKG_CONFIG = PKG_CONFIG_PATH="$(HOMEBREW_PC_PATH):$$PKG_CONFIG_PATH" pkg-config
else
  PKG_CONFIG = pkg-config
endif

# BLAS/LAPACK for fast machine-precision linear-algebra kernels
# (machine-precision QRDecomposition; later: Inverse, LinearSolve, Det,
# Eigenvalues, LeastSquares, SVD).  Four-tier autodetection:
#   1. Darwin              -> Apple Accelerate framework (zero install).
#   2. pkg-config lapacke  -> use pkg-config flags (OpenBLAS / MKL / etc).
#   3. /usr/include/lapacke.h or $(SYSTEM_PREFIX)/include/lapacke.h
#                          -> link -llapacke -llapack -lblas.
#   4. nothing found       -> warn and override USE_LAPACK := 0.
# The BLAS provider is found by a LINK PROBE, not assumed: Mathilda calls
# cblas_* directly, and which library exports them varies by distro -- the
# reference libblas on Debian/Ubuntu, libopenblas, a separate libcblas, or
# FlexiBLAS on Fedora/RHEL, whose netlib libblas.so has no CBLAS at all, so the
# old blind `-lblas` failed there with "undefined reference to cblas_ddot ...
# DSO missing from command line" (issues #94, #97). lapack_blas_probe links a
# '#'-free program calling cblas_ddot + dgetrf_ against the LAPACK libs plus
# each candidate in turn and prints the first that links ("ok:" when the LAPACK
# libs already carry CBLAS); empty means nothing linked, and the old -lblas
# default is kept so the final link reports the real error.
# When USE_LAPACK is off, machine-precision linalg falls back to the
# MPFR / symbolic path with a one-time runtime warning.  This matches
# the existing USE_MPFR=0 / USE_ECM=0 graceful-degrade policy so that
# `git clone && make` always succeeds, no matter the host environment.
# Tier-3 header locations: /usr/include on Unix only (on MSYS2 it is the Cygwin
# tree), plus SYSTEM_PREFIX's include dir when it exists.
LAPACKE_HDR_DIRS := $(if $(filter Windows,$(BUILD_PLATFORM)),,/usr/include) $(SYSPREFIX_CFLAGS:-I%=%)
LAPACK_BLAS_CANDIDATES := -lflexiblas -lopenblas -lcblas -lblas '-lcblas -lblas'
lapack_blas_probe = $(shell f=/tmp/mathilda_blasprobe_$$$$; \
    printf '%s\n' 'double cblas_ddot(int, const double*, int, const double*, int);' \
      'void dgetrf_(int*, int*, double*, int*, int*, int*);' \
      'int main(void){double a = 2; int n = 1, p, i; dgetrf_(&n, &n, &a, &n, &p, &i); return cblas_ddot(1, &a, 1, &a, 1) < 0;}' \
      > $$f.c 2>/dev/null && \
    for b in '' $(LAPACK_BLAS_CANDIDATES); do \
      $(CC) $$f.c -o $$f $(1) $$b 2>/dev/null && { echo "ok:$$b"; break; }; \
    done; rm -f $$f.c $$f)
# $(call lapack_with_blas,<lapack libs>) -> those libs plus the BLAS provider.
lapack_with_blas = $(1) $(if $(LAPACK_BLAS_FOUND),$(patsubst ok:%,%,$(LAPACK_BLAS_FOUND)),-lblas)

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
    # Debian/Ubuntu's lapacke.pc exposes only `-llapacke` and hides the Fortran
    # LAPACK + reference CBLAS it depends on behind Requires.private, which a
    # non-static `pkg-config --libs` never expands. Mathilda calls the Fortran
    # LAPACK routines (dgeev_, zgesv_, dgesdd_, ...) and cblas_* (cblas_dgemm,
    # cblas_ddot) DIRECTLY — it uses no LAPACKE C wrappers — so a pkg-config-only
    # link fails on Ubuntu with hundreds of undefined references even though the
    # headers resolve fine. Append the transitive -llapack + the probed BLAS
    # provider (once each), after -llapacke so ld's left-to-right resolution
    # keeps each provider behind its consumer. Homebrew's .pc and the Accelerate
    # path never hit this. Same fixup idiom as the FLINT block.
    LAPACKE_LIBS := $(shell pkg-config --libs lapacke)
    ifeq ($(filter -llapack,$(LAPACKE_LIBS)),)
      LAPACKE_LIBS += -llapack
    endif
    LAPACK_BLAS_FOUND := $(call lapack_blas_probe,$(LAPACKE_LIBS))
    LDFLAGS += $(call lapack_with_blas,$(LAPACKE_LIBS))
  else ifneq ($(wildcard $(addsuffix /lapacke.h,$(LAPACKE_HDR_DIRS))),)
    CFLAGS  += -DUSE_LAPACK
    LAPACK_BLAS_FOUND := $(call lapack_blas_probe,$(SYSPREFIX_LDFLAGS) -llapacke -llapack)
    LDFLAGS += $(call lapack_with_blas,-llapacke -llapack)
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
    # graphics_export_raster screens for a usable GUI session before InitWindow
    # (which segfaults inside GLFW otherwise) via CGSessionCopyCurrentDictionary;
    # that symbol lives in CoreGraphics, and CFRelease in CoreFoundation.
    ifeq ($(BUILD_PLATFORM),Darwin)
      LDFLAGS += -framework CoreGraphics -framework CoreFoundation
    endif
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
SRC = $(wildcard $(SRC_DIR)/*.c) $(wildcard $(SRC_DIR)/solve/*.c) $(wildcard $(SRC_DIR)/list/*.c) $(wildcard $(SRC_DIR)/ml/*.c) $(wildcard $(SRC_DIR)/linalg/*.c) $(wildcard $(SRC_DIR)/numbertheory/*.c) $(wildcard $(SRC_DIR)/numbertheory/factorinteger/*.c) $(wildcard $(SRC_DIR)/poly/*.c) $(wildcard $(SRC_DIR)/simp/*.c) $(wildcard $(SRC_DIR)/stats/*.c) $(wildcard $(SRC_DIR)/calculus/*.c) $(wildcard $(SRC_DIR)/sum/*.c) $(wildcard $(SRC_DIR)/product/*.c) $(wildcard $(SRC_DIR)/special_functions/*.c) $(wildcard $(SRC_DIR)/compile/*.c) $(wildcard $(SRC_DIR)/numerical_calculus/*.c) $(wildcard $(SRC_DIR)/numerical_roots/*.c) $(wildcard $(SRC_DIR)/graphics/*.c) $(wildcard $(SRC_DIR)/graph/*.c) $(wildcard $(SRC_DIR)/strings/*.c) $(wildcard $(SRC_DIR)/strings/regex/*.c) $(wildcard $(SRC_DIR)/io/*.c) $(wildcard $(SRC_DIR)/bitwise/*.c) $(wildcard $(SRC_DIR)/sparse/*.c)
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
DEP = $(OBJ:.o=.d) $(FFI_OBJ:.o=.d)
# MinGW GCC appends .exe to an extension-less -o target, so the produced file is
# Mathilda.exe; naming TARGET to match keeps the rule's timestamp honest (no
# perpetual relink) and gives scripts a deterministic path.
ifeq ($(BUILD_PLATFORM),Windows)
TARGET = Mathilda.exe
else
TARGET = Mathilda
endif

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

# ---------------------------------------------------------------------------
# Shared library (issue #88): the same LIB_OBJ, linked as a shared object whose
# ONLY dynamic exports are the MATHILDA_FFI_API entry points of
# src/ffi/mathilda_ffi.h (everything else is -fvisibility=hidden, see CFLAGS;
# on Windows the dllexport in that header switches off MinGW's export-all).
# The Mathilda executable does NOT link it — it keeps its objects — so `make`
# (= `all`) is unchanged; `make lib` builds both libraries.
#
# SHLIB_ABI is the soname / compatibility key. It tracks the C ABI of
# mathilda_ffi.h, NOT $VersionNumber (which moves on every commit): bump it only
# for an incompatible change to that header.
#
# The library is linked against the same autodetected dependencies as the
# executable, so a consumer needs only -lmathilda; a static consumer gets the
# transitive list from the installed mathilda.pc (`pkg-config --static`).
SHLIB_ABI = 0
MATHILDA_VERSION_STRING := $(shell sed -n 's/^\#define[[:space:]]*MATHILDA_VERSION_STRING[[:space:]]*"\(.*\)".*/\1/p' $(SRC_DIR)/version.h)
ifeq ($(BUILD_PLATFORM),Darwin)
  SHLIB          = libmathilda.$(SHLIB_ABI).dylib
  SHLIB_DEVLINK  = libmathilda.dylib
  # @rpath in the build tree; `make install` rewrites the id to the absolute
  # installed path (the Homebrew convention), so plain -L consumers just work.
  # (No -current_version: Mach-O packs it as X.Y.Z with Y,Z <= 255, which a
  # 0.408-style release number does not fit; the ABI lives in the file name.)
  SHLIB_LDFLAGS  = -dynamiclib -install_name @rpath/$(SHLIB)
else ifeq ($(BUILD_PLATFORM),Windows)
  SHLIB          = libmathilda.dll
  SHLIB_IMPLIB   = libmathilda.dll.a
  SHLIB_LDFLAGS  = -shared -Wl,--out-implib,$(SHLIB_IMPLIB)
else
  SHLIB          = libmathilda.so.$(SHLIB_ABI)
  SHLIB_DEVLINK  = libmathilda.so
  # --no-undefined: an unresolved dependency fails HERE, not at a consumer's
  # dlopen time.
  SHLIB_LDFLAGS  = -shared -Wl,-soname,$(SHLIB) -Wl,--no-undefined
endif

$(SHLIB): $(LIB_OBJ)
	$(CC) $(CFLAGS) $(SHLIB_LDFLAGS) -o $@ $(LIB_OBJ) $(LDFLAGS_DEDUP) $(EXTRA_LIBS)
ifdef SHLIB_DEVLINK
	ln -sf $(SHLIB) $(SHLIB_DEVLINK)
endif

lib: libmathilda.a $(SHLIB)

# -MMD writes a .d file next to each .o listing the (non-system) headers it
# includes; -MP adds a phony target for each header so deleting a header does
# not break the build with a "no rule to make target" error. The .d files are
# consumed by the `-include $(DEP)` line below.
$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -f $(OBJ) $(DEP) $(TARGET)
	rm -f libmathilda.a $(SHLIB) $(SHLIB_DEVLINK) $(SHLIB_IMPLIB)
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

# ---------------------------------------------------------------------------
# Installation (issue #77).
#
# The Mathilda executable is NOT self-contained: its kernel is part C, part its
# own language, and the `src/internal/` tree of `.m` bootstrap files is loaded
# at RUNTIME (init.m/bessel.m/distributions.m at startup; simp/FullSimplify.m +
# simp/transforms/*.m, CRCMathTablesIntegrals.m, and mixed/ParallelMixed.m
# lazily). So an install must ship BOTH the binary and that module tree.
#
# We install into the FHS layout the loader already searches EXE-RELATIVE
# (src/loadmodule.c: <exe>/../share/mathilda/internal/), so the result is
# relocatable and needs no $MATHILDA_HOME: wherever `bin/Mathilda` ends up, it
# finds `../share/mathilda/internal/`. That also makes staged/packaged installs
# (Debian, SageMath, …) work via DESTDIR with no baked-in absolute path.
#
# Since issue #88 the install also ships the EMBEDDING API, so a downstream
# binding (passagemath's Cython module, a Jupyter kernel, any C host) can build
# against an installed Mathilda without the source tree:
#   $(includedir)/mathilda_ffi.h      the C API (self-contained header)
#   $(libdir)/libmathilda.a           static kernel (PIC: links into a .so)
#   $(libdir)/libmathilda.so.0 (+.so) shared kernel, exports only mathilda_ffi_*
#     (macOS libmathilda.0.dylib; Windows bin/libmathilda.dll + lib/*.dll.a)
#   $(pkgconfigdir)/mathilda.pc       -I/-L/-lmathilda, the build's transitive
#                                     libraries as Libs.private (for
#                                     `pkg-config --static`), and `internaldir`
# The shared library locates the module tree LIBRARY-relative (dladdr:
# <libdir>/../share/mathilda/internal/), so it needs no $MATHILDA_HOME in a host
# process either; a static consumer reads `pkg-config --variable=internaldir
# mathilda` and calls mathilda_ffi_set_home().
#
# Variables (GNU-standard). Note the LOWER-case `prefix`: giving the UPPER-case
# PREFIX a default here would trip the opt-in `ifdef PREFIX` above and bake
# -DMATHILDA_PREFIX into every ordinary build. `$(or ...)` lets both
# `make install PREFIX=/opt` and `make install prefix=/opt` work, while a bare
# `make` still defines nothing. Shared libraries (GMP, MPFR, Readline, and any
# optional FLINT/ECM/PCRE2/FFTW/Raylib/LAPACK detected at build time) are a
# separate runtime dependency on the target — see the README install section.
prefix          ?= $(or $(PREFIX),/usr/local)
exec_prefix     ?= $(prefix)
bindir          ?= $(exec_prefix)/bin
datarootdir     ?= $(prefix)/share
datadir         ?= $(datarootdir)
pkgdatadir      ?= $(datadir)/mathilda
libdir          ?= $(exec_prefix)/lib
includedir      ?= $(prefix)/include
pkgconfigdir    ?= $(libdir)/pkgconfig
INSTALL         ?= install
INSTALL_PROGRAM ?= $(INSTALL)
INSTALL_DATA    ?= $(INSTALL) -m 644

# Libs.private for mathilda.pc: exactly what the kernel was linked against,
# minus linker-only directives (-Wl,..., e.g. macOS header padding).
PC_LIBS_PRIVATE = $(strip $(filter-out -Wl%,$(LDFLAGS_DEDUP)) $(EXTRA_LIBS))

install: $(TARGET) libmathilda.a $(SHLIB)
	$(INSTALL) -d "$(DESTDIR)$(bindir)"
	$(INSTALL_PROGRAM) $(TARGET) "$(DESTDIR)$(bindir)/$(TARGET)"
	@# Recreate the runtime module tree, preserving simp/, simp/transforms/,
	@# mixed/. Only *.m files are copied (editor backups etc. are skipped).
	cd $(SRC_DIR)/internal && find . -type d \
	    -exec $(INSTALL) -d "$(DESTDIR)$(pkgdatadir)/internal/{}" \;
	cd $(SRC_DIR)/internal && find . -type f -name '*.m' \
	    -exec $(INSTALL_DATA) '{}' "$(DESTDIR)$(pkgdatadir)/internal/{}" \;
	@# Embedding API (issue #88): header, static + shared library, pkg-config.
	$(INSTALL) -d "$(DESTDIR)$(includedir)" "$(DESTDIR)$(libdir)" "$(DESTDIR)$(pkgconfigdir)"
	$(INSTALL_DATA) $(SRC_DIR)/ffi/mathilda_ffi.h "$(DESTDIR)$(includedir)/mathilda_ffi.h"
	$(INSTALL_DATA) libmathilda.a "$(DESTDIR)$(libdir)/libmathilda.a"
ifeq ($(BUILD_PLATFORM),Windows)
	$(INSTALL_PROGRAM) $(SHLIB) "$(DESTDIR)$(bindir)/$(SHLIB)"
	$(INSTALL_DATA) $(SHLIB_IMPLIB) "$(DESTDIR)$(libdir)/$(SHLIB_IMPLIB)"
else
	$(INSTALL_PROGRAM) $(SHLIB) "$(DESTDIR)$(libdir)/$(SHLIB)"
	ln -sf $(SHLIB) "$(DESTDIR)$(libdir)/$(SHLIB_DEVLINK)"
endif
ifeq ($(BUILD_PLATFORM),Darwin)
	@# Absolute install name for the installed copy (the build tree keeps
	@# @rpath), then re-sign: rewriting the id invalidates the ad-hoc signature
	@# and arm64 dyld refuses to load an unsigned dylib.
	install_name_tool -id "$(libdir)/$(SHLIB)" "$(DESTDIR)$(libdir)/$(SHLIB)"
	codesign --force --sign - "$(DESTDIR)$(libdir)/$(SHLIB)"
endif
	printf '%s\n' \
	    'prefix=$(prefix)' \
	    'libdir=$(libdir)' \
	    'includedir=$(includedir)' \
	    'internaldir=$(pkgdatadir)/internal' \
	    '' \
	    'Name: mathilda' \
	    'Description: Mathilda computer algebra kernel, embeddable via mathilda_ffi.h' \
	    'URL: https://github.com/stblake/mathilda' \
	    'Version: $(MATHILDA_VERSION_STRING)' \
	    'Libs: -L$${libdir} -lmathilda' \
	    'Libs.private: $(PC_LIBS_PRIVATE)' \
	    'Cflags: -I$${includedir}' \
	    > "$(DESTDIR)$(pkgconfigdir)/mathilda.pc"
	chmod 644 "$(DESTDIR)$(pkgconfigdir)/mathilda.pc"
	@echo
	@echo "Installed:"
	@echo "  executable  -> $(DESTDIR)$(bindir)/$(TARGET)"
	@echo "  module tree -> $(DESTDIR)$(pkgdatadir)/internal/"
	@echo "  C API       -> $(DESTDIR)$(includedir)/mathilda_ffi.h"
	@echo "  libraries   -> $(DESTDIR)$(libdir)/libmathilda.a, $(SHLIB)"
	@echo "  pkg-config  -> $(DESTDIR)$(pkgconfigdir)/mathilda.pc"
	@echo "The binary and the shared library find their modules relative to"
	@echo "their own path, so no MATHILDA_HOME is needed. Ensure the shared"
	@echo "libraries they link (see 'ldd'/'otool -L $(TARGET)') are present."

uninstall:
	rm -f "$(DESTDIR)$(bindir)/$(TARGET)"
	rm -rf "$(DESTDIR)$(pkgdatadir)"
	rm -f "$(DESTDIR)$(includedir)/mathilda_ffi.h"
	rm -f "$(DESTDIR)$(libdir)/libmathilda.a"
	rm -f "$(DESTDIR)$(pkgconfigdir)/mathilda.pc"
ifeq ($(BUILD_PLATFORM),Windows)
	rm -f "$(DESTDIR)$(bindir)/$(SHLIB)" "$(DESTDIR)$(libdir)/$(SHLIB_IMPLIB)"
else
	rm -f "$(DESTDIR)$(libdir)/$(SHLIB)" "$(DESTDIR)$(libdir)/$(SHLIB_DEVLINK)"
endif

# `make check-install` — stage an install and prove the embedding API works from
# it: a C host linked against the shared library (no MATHILDA_HOME: the library
# must find its module tree itself) and statically via `pkg-config --static`,
# the static archive linked into a shared object (the Cython case; fails without
# -fPIC), and the shared library's dynamic exports limited to mathilda_ffi_*.
check-install: $(TARGET) libmathilda.a $(SHLIB)
	CC="$(CC)" MAKE="$(MAKE)" sh tools/check_install.sh

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

# Backstop gate: re-run every published In[]/Out[] example through the current
# binary and fail on any mismatch. `make docs` verifies examples at *generation*
# time; this guards the committed tree against going stale when the binary's
# behaviour changes and nobody regenerated. Needs ./Mathilda.
check-docs-examples: $(TARGET)
	python3 site/verify_docs_examples.py

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
check-packed-aware:
	python3 tools/check_packed_aware.py

# `make check-messages` — does every user-facing message support Quiet[] and
# Check[]? Those two are a two-step convention (mth_msg_note_fired for Check,
# mth_msg_suppressed for Quiet), not a choke-point, so a raw
# `fprintf(stderr, "Head::tag: ...")` bypasses BOTH — Check silently returns the
# wrong branch and Quiet leaks. This reads the raw Head::tag stderr writes out of
# the source and diffs them against a checked-in EXEMPT (permanent) + BASELINE
# (the shrinking migration backlog); it ratchets, failing on a new bypassing site
# and on a migrated entry left stale. Source-level, so same "needs python3, not
# part of `all`" status as check-c99.
check-messages:
	python3 tools/check_message_routing.py

# `make check-feyncalc` — FeynCalc's ~14k upstream unit tests, run UNMODIFIED
# under a minimal MUnit`Test (tools/feyncalc/run_fc_tests.m), one process per
# .test file, ratcheted against the local pass baseline: a test that passed must
# keep passing. The vendored package (packages/feyncalc/) is local and
# untracked; without it the check prints a skip and succeeds.
# FC_ARGS passes options through, e.g. FC_ARGS="--dir Lorentz --triage".
check-feyncalc: Mathilda
	python3 tools/check_feyncalc.py $(FC_ARGS)

# `make check-corpus-indvar` — does every DSolve corpus record ask the question
# the book printed? A converted record can PARSE cleanly and still encode a
# DIFFERENT equation, because the converter promoted a PARAMETER to the
# independent variable (`y''-2ay'+a^2y=0` read as an ODE in `a`). The harness
# then back-substitutes into the garbled equation, so the record scores PASS —
# or UNEVAL — for a question nobody asked, and neither the residual nor the gate
# baseline can see it. Eighteen records across nine sections were in that state;
# hand-grepping found the class six times (M42/M44/M48/M49/M61/M63) and each new
# section re-opened it. Two rules read off the record itself: a `_missing_x`
# classification forbids the independent variable from occurring in the body, and
# the independent variable must be a letter that names a variable by convention.
# Assert-empty, with an EXEMPT list for a deliberate case.
check-corpus-indvar:
	python3 tools/check_corpus_indvar.py

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

# `make check-pipe-protocol` — does the NDJSON pipe protocol answer as documented in src/repl.c?
#
# Two request modes share one process: plain requests (the site generator and audit tools, whose
# behaviour must not move) and notebook cells (several statements, Print output and messages as
# "stream"/"message" lines). Needs the built ./Mathilda; SKIPs without it.
check-pipe-protocol:
	python3 tools/check_pipe_protocol.py

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

# `make check-refine-stress` — the Refine assumption-engine stress suite.
#
# A 159-case first-principles corpus (tests/refine_stress/) covering every
# documented Refine surface plus adversarial soundness/robustness probes. Each
# case is decided by structural equality (SameQ) inside Mathilda, so printer
# form never matters, and each runs in its own process under a timeout so a hang
# or crash is contained. Complements tests/test_refine.c (which asserts the
# PASSing subset in C) by ALSO holding the cases Refine cannot yet reduce.
#
# Ratchets against BASELINE_UNCHANGED in run_refine_stress.py: fails on a
# soundness/hang/crash, on a reduction that regressed, and on a baseline gap
# that started passing (so the baseline is kept honest). Needs ./Mathilda.
check-refine-stress:
	python3 tests/refine_stress/run_refine_stress.py

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
# benchmarks/results/<date>.json. Experiments with a .gp / .mac half also run
# in PARI/GP / Maxima.
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

# Report the fully expanded compile and link flags. The Windows CI leg uses this
# to assert that no /usr search path reaches a native MinGW build (the MSYS2
# Cygwin-header trap, passagemath PR #2966).
print-flags:
	@echo "CFLAGS = $(CFLAGS)"
	@echo "LDFLAGS = $(LDFLAGS)"

.PHONY: all lib clean install uninstall check-install docs docs-build docs-serve check-docs-examples check-c99 check-feyncalc check-interval check-packed-aware \
        check-messages check-corpus-indvar \
        check-array-exactness check-nd-surfaces check-compile-coverage \
        check-refine-stress \
        check-fastpath-sweep check-menu-ids check-pipe-protocol bench-gap check-diophantine-heldout print-cc print-flags

# Pull in the auto-generated header dependencies. The leading `-` silences the
# "no such file" notice on a fresh tree (no .d files exist until the first
# compile). Placed last so it never becomes the default goal.
-include $(DEP)
