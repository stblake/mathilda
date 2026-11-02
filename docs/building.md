# Building Mathilda

This document covers the Mathilda **kernel** — the `Mathilda` executable and the
embeddable `libmathilda` — on macOS, Linux and Windows: dependencies, build
flags, installation, and embedding. The optional **desktop app** (a Tauri front
end that runs the kernel as a sidecar) is covered in
[Desktop app](#desktop-app-tauri), and common problems in
[Troubleshooting](#troubleshooting).

---

## Kernel

### Dependencies

| Library | Package (Debian / Homebrew) | Status |
|---------|-----------------------------|--------|
| GCC | `build-essential` / `gcc` | **required** — the build rejects Apple/LLVM clang, which hides the glibc-portability warnings the Linux CI gates on |
| GMP | `libgmp-dev` / `gmp` | **required** — arbitrary-precision integers |
| MPFR | `libmpfr-dev` / `mpfr` | on by default — arbitrary-precision reals (`USE_MPFR=0` to omit) |
| GNU Readline | `libreadline-dev` / `readline` | auto-detected — REPL line editing, history, multi-line input |
| FLINT ≥ 3.0 | `libflint-dev` / `flint` | optional — fast, rigorous polynomial arithmetic over ℚ and number fields, rigorous `acb` numerics |
| GMP-ECM | `libecm-dev` / `gmp-ecm` | optional — Elliptic Curve Method factorization |
| LAPACK / BLAS | `liblapacke-dev libopenblas-dev` / *(Accelerate, built in)* | optional — machine-precision linear algebra |
| FFTW ≥ 3 | `libfftw3-dev` / `fftw` | optional — fast `Fourier`/`FourierDCT`/`FourierDST` (naive O(n²) otherwise) |
| PCRE2 | `libpcre2-dev` / `pcre2` | optional — `RegularExpression` and the regex string functions |
| Raylib ≥ 4.0 | `libraylib-dev` / `raylib` | optional — interactive plot windows (text placeholder otherwise) |
| CMake | `cmake` / `cmake` | only for the test suite |

Every optional backend is auto-detected and **degrades gracefully** when absent.
POSIX threads (used to parallelise packed-array kernels) need no extra package.

**macOS (Homebrew):**

```bash
brew install gcc gmp mpfr readline cmake
brew install flint gmp-ecm fftw pcre2 raylib      # optional
```

LAPACK/BLAS need not be installed on macOS — the build auto-detects Apple's
**Accelerate** framework.

**Debian / Ubuntu:**

```bash
sudo apt install build-essential libgmp-dev libmpfr-dev libreadline-dev cmake
sudo apt install libflint-dev libecm-dev liblapacke-dev libopenblas-dev \
                 libfftw3-dev libpcre2-dev libraylib-dev   # optional
```

**Fedora / RHEL:** `gmp-devel`, `mpfr-devel`, `readline-devel`, `flint-devel`
(≥ 3.0), `gmp-ecm-devel`, `lapack-devel` + `lapacke-devel` (BLAS via FlexiBLAS,
detected automatically), `fftw-devel`, `pcre2-devel`, plus `cmake`.

> **FLINT versions.** Mathilda requires **FLINT ≥ 3.0** (the release that merged
> ANTIC for number-field arithmetic). Distributions that only package FLINT 2.x —
> e.g. Ubuntu 22.04 or Debian Bullseye — are detected as too old and the build
> falls back to `USE_FLINT=0`. See [Troubleshooting](#troubleshooting) for
> `tools/install-flint.sh` and the conda-forge route.

**Windows (MSYS2):** Mathilda builds as a native Windows executable
(`Mathilda.exe`) with the MSYS2 GCC toolchain — no WSL or Cygwin needed.

1. Install [MSYS2](https://www.msys2.org/) (default location `C:\msys64`).
2. Open the **MSYS2 UCRT64** shell and install the toolchain and libraries:

   ```bash
   pacman -Syu                        # update MSYS2 first (may ask you to reopen the shell)
   pacman -S --needed git make \
     mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-gmp \
     mingw-w64-ucrt-x86_64-mpfr mingw-w64-ucrt-x86_64-pkgconf
   pacman -S --needed mingw-w64-ucrt-x86_64-cmake      # optional: test suite
   ```

   This is the configuration the Windows CI job builds and smoke-tests on every
   push. The older **MINGW64** environment works too (use the
   `mingw-w64-x86_64-` package prefix). CLANG64/CLANGARM64 are unsupported.

MPFR is on by default and, unlike the auto-detected backends, does not degrade
gracefully: install it, or build with `USE_MPFR=0`. Some features are reduced
on Windows:

* **GNU Readline** isn't available under MinGW and is compiled out
  automatically; the REPL works without line editing, history or multi-line
  continuation.
* **Threads** are off (`USE_THREADS=0`), so packed-array kernels run on one core.
* **`TimeConstrained`** is cooperative only: the evaluator loop (and the
  deadline polls in heavy kernels) enforce the budget, but a single long-running
  C builtin that never returns to the loop cannot be interrupted. Timing,
  `Pause` and memory reporting all work (via `src/mtime.c` and `src/meminfo.c`).
* **Optional backends** (FLINT, LAPACK, ECM, FFTW, PCRE2, Raylib) are switched
  off in the recommended build below, because CI doesn't test them on Windows.

### Building

The `makefile` auto-discovers `src/**/*.c` and links `./Mathilda`
(`-std=c99 -O3`):

```bash
make -j$(nproc)
./Mathilda                       # interactive REPL
./Mathilda -file script.m        # run a script, printing only what it Prints
./Mathilda --help                # all options
```

Disable optional backends for a leaner binary, e.g.
`make -j$(nproc) USE_LAPACK=0 USE_GRAPHICS=0`. On Windows, from the UCRT64 shell,
use the flags the CI job uses:

```bash
make -j4 USE_ECM=0 USE_READLINE=0 USE_GRAPHICS=0 \
         USE_LAPACK=0 USE_FLINT=0 USE_REGEX=0 USE_FFTW=0
./Mathilda.exe
```

### Build flags reference

| Flag | Default | Effect |
|------|---------|--------|
| `USE_MPFR` | `1` | Arbitrary-precision reals: `N[expr, prec]`, `Precision`/`Accuracy`, precision literals. `USE_MPFR=0` disables them (with a runtime warning). |
| `USE_FLINT` | `1` | FLINT (≥ 3.0) kernels: multivariate polynomial GCD/factoring over ℚ, univariate GCD/factoring over number fields ℚ(α) (via the `gr` layer + ANTIC), the finite-field workhorse behind parametric ℚ(t)(α) work, and rigorous `acb` numerics. Auto-detected via `pkg-config` with a ≥ 3.0 floor; falls back to the classical (slower but still rigorous) path. |
| `USE_LAPACK` | `1` | Machine-precision linear algebra. Auto-detected: Apple **Accelerate** on macOS, `lapacke`/`lapack` + a CBLAS on Linux. Falls back to the pure-C path. |
| `USE_ECM` | `1` | Elliptic Curve Method factorization via the system GMP-ECM library. Auto-detected via a compile-link probe; auto-disabled when not found. |
| `USE_FFTW` | `1` | FFTW-backed `Fourier`/`FourierDCT`/`FourierDST`. Auto-detected via `pkg-config fftw3`; falls back to a naive O(n²) transform. |
| `USE_REGEX` | `1` | PCRE2-backed `RegularExpression`, `StringMatchQ`, `StringCases`, `StringReplace`, `StringSplit`. Auto-detected via `pkg-config libpcre2-8`; those heads warn and stay unevaluated without it. |
| `USE_THREADS` | `1` | POSIX-thread parallelism for element-wise packed-array kernels on large arrays. Off on Windows. |
| `USE_GRAPHICS` | `1` | Interactive 2D/3D plot windows via Raylib, auto-detected via `pkg-config raylib`. Without it, `Plot`/`Show`/… print a text placeholder. |
| `USE_READLINE` | `1` | GNU Readline for the REPL. Auto-detected; off on Windows. On macOS the makefile prefers Homebrew's GNU Readline over Apple's libedit shim, which cannot edit multi-line input. |
| `HOMEBREW=<prefix>` | auto (`brew --prefix`) | Homebrew prefix whose `include`, `lib`, and `lib/pkgconfig` directories are added to the search paths for the optional dependencies above. Auto-detected with `brew --prefix`, so the Apple-Silicon `/opt/homebrew`, the Intel `/usr/local`, a Linuxbrew prefix, and any custom install all work without editing the makefile. Pass it explicitly to force a prefix, or **opt out of all Homebrew discovery** with `HOMEBREW=` (empty) or `HOMEBREW=/nonexisting`. Honoured identically by the makefile and the CMake test build. |
| `SYSTEM_PREFIX=<prefix>` | `/usr/local` (empty on Windows) | Conventional site prefix whose `include` and `lib` directories (when they exist) are added to the compile/link lines and to the GMP-ECM and LAPACK probes. Empty by default on native Windows, where MSYS2's `/usr` is the Cygwin runtime tree. Opt out with `SYSTEM_PREFIX=` (empty). `/usr/include` itself is never added. |
| `EXTRA_LIBS="…"` | empty | Extra link libraries, e.g. a BLAS the probe does not know. |

`HOMEBREW` and `SYSTEM_PREFIX` are read from the command line or the environment (`make HOMEBREW=/opt/homebrew`,
`export HOMEBREW=`); the CMake test build also accepts `-DHOMEBREW=...`. A packager building
against a private prefix can pass `HOMEBREW= SYSTEM_PREFIX=` so that the makefile searches only
the paths it is given (`CFLAGS`, `pkg-config`). `make print-flags` prints the fully expanded
compile and link flags.

### Installing

**The `Mathilda` executable is not self-contained.** Its kernel is written partly
in C and partly in its own language: a tree of `.m` bootstrap modules under
`src/internal/` is loaded *at runtime* — some at startup, some lazily on first
use. Copy the binary alone and it will start, but `Simplify`/`FullSimplify`, the
CRC integral table behind `Integrate`, the Bessel and distribution rules, and the
`ParallelMixedTower` integration method will all be missing (you'll see a
`LoadModule::nofile` diagnostic on stderr). An install must therefore provide the
runtime module tree (the whole `src/internal/` directory) and the shared
libraries the binary was linked against.

#### `make install` (recommended)

```bash
make                                   # build (PREFIX optional here — see note)
sudo make install PREFIX=/usr/local    # install binary + module tree under PREFIX
```

| What | Where |
|------|-------|
| the `Mathilda` executable | `$(PREFIX)/bin/Mathilda` |
| the runtime `.m` module tree | `$(PREFIX)/share/mathilda/internal/` |
| the embedding C API header | `$(PREFIX)/include/mathilda_ffi.h` |
| the static kernel library | `$(PREFIX)/lib/libmathilda.a` |
| the shared kernel library | `$(PREFIX)/lib/libmathilda.so.0` + `libmathilda.so` (macOS: `libmathilda.0.dylib` + `libmathilda.dylib`; Windows: `bin/libmathilda.dll` + `lib/libmathilda.dll.a`) |
| the pkg-config file | `$(PREFIX)/lib/pkgconfig/mathilda.pc` |

The binary locates its modules **relative to its own path**
(`<bindir>/../share/mathilda/internal/`), so the install is **relocatable** and
needs no environment variable. `PREFIX` on the **build** line only bakes in a
compile-time absolute fallback (`-DMATHILDA_PREFIX`) for environments where the
executable can't determine its own location; `PREFIX` on the **install** line is
what chooses the destination. The standard GNU variables are honoured for
packaging:

```bash
make install DESTDIR=/tmp/stage PREFIX=/usr    # staged install (Debian, SageMath, …)
make install prefix=/opt/mathilda              # or override bindir=/… datadir=/… directly
make uninstall PREFIX=/usr/local               # remove everything `make install` put there
```

Because resolution is exe-relative, `DESTDIR` staging bakes in no absolute path.

#### How the binary finds its modules

`src/loadmodule.c` searches these locations in order (first hit wins):

1. **`$MATHILDA_HOME`** — point a binary at any module tree:
   `MATHILDA_HOME=/opt/app/share/mathilda/internal ./Mathilda`.
2. **`<exe>/src/internal/`** — running straight from a build tree (the repo root).
3. **`<exe>/../share/mathilda/internal/`** — the layout `make install` creates.
   The same path is also tried relative to the *shared library* holding the
   kernel, so an installed `libmathilda` loaded into any host finds its modules.
4. **compile-time `-DMATHILDA_PREFIX`** — `<prefix>/share/mathilda/internal/`,
   baked in when you build with `make PREFIX=…`.

#### Manual install

```bash
mkdir -p /opt/app/bin /opt/app/share/mathilda
cp Mathilda /opt/app/bin/
cp -R src/internal /opt/app/share/mathilda/internal
/opt/app/bin/Mathilda          # finds ../share/mathilda/internal automatically
```

Or skip the layout entirely and set `MATHILDA_HOME` to wherever you placed
`src/internal`.

#### Shared libraries

The binary is dynamically linked against every backend detected at build time, so
the **target machine needs those same libraries** (`ldd ./Mathilda` on Linux,
`otool -L ./Mathilda` on macOS):

* **Always:** GMP, MPFR, GNU Readline, plus the C runtime.
* **When detected at build:** FLINT, GMP-ECM, PCRE2, FFTW, Raylib, LAPACK/BLAS.
  Turn optional backends off for fewer dependencies
  (e.g. `make USE_GRAPHICS=0 USE_FFTW=0 USE_ECM=0`).
* **Linux / macOS:** install the *runtime* packages on the target — the
  `-dev`/`-devel` headers are only needed to build. Accelerate ships with macOS.
* **Windows:** put the DLLs `Mathilda.exe` depends on (GMP, MPFR, the GCC
  runtime, …) beside the executable or on `PATH` (`C:\msys64\ucrt64\bin`, or
  `C:\msys64\mingw64\bin` for a MinGW64 build); `ldd ./Mathilda.exe` in the
  MSYS2 shell lists them. The desktop app bundles the sidecar and its DLLs.

### Embedding library (libmathilda)

The kernel can run **in-process** in another program (a Python/Cython extension,
a Jupyter kernel, a mobile app). The same objects, minus `repl.o`, build a
library behind the small C API in `src/ffi/mathilda_ffi.h`:
`mathilda_ffi_init`, `mathilda_ffi_eval` / `_eval_latex` / `_eval_json` /
`_eval_cell`, `mathilda_ffi_free`, and a few editor helpers.

```bash
make lib                     # libmathilda.a + the shared library, in the tree
make install prefix=/opt/m   # also installs include/mathilda_ffi.h, lib/libmathilda.*,
                             #   lib/pkgconfig/mathilda.pc
make check-install           # install to a scratch prefix and link C hosts against it
```

| Platform | Shared library | Notes |
|----------|----------------|-------|
| Linux    | `libmathilda.so.0` (soname) + `libmathilda.so` | linked `-Wl,--no-undefined` |
| macOS    | `libmathilda.0.dylib` + `libmathilda.dylib` | `@rpath` id in the tree; `make install` rewrites it to the absolute path and re-signs |
| Windows  | `libmathilda.dll` + import lib `libmathilda.dll.a` | DLL goes to `bin/` |

The trailing `0` is `SHLIB_ABI` in the makefile. It tracks the C ABI of
`mathilda_ffi.h`, not `$VersionNumber`, so bump it only for an incompatible change
to that header. Every object is compiled `-fPIC -fvisibility=hidden` (except on
Windows). That lets the static archive link into a Python extension module, and
it limits the shared library's exports to the `MATHILDA_FFI_API` functions. On
Windows the `__declspec(dllexport)` in the header does the same job.

A static consumer gets its link line from `pkg-config --static --libs mathilda`
(the dependencies detected at build time). It also has to call
`mathilda_ffi_set_home()` with `pkg-config --variable=internaldir mathilda`. The
shared library finds its module tree relative to its own location.

---

## Desktop app (Tauri)

The desktop app wraps the kernel as a sidecar process (it talks to it over the
NDJSON pipe protocol), so first install the kernel dependencies above, then the
Tauri toolchain below.

### Prerequisites (all platforms)

- **Rust toolchain** — install via [rustup.rs](https://rustup.rs/)
- **Node.js 18+** — for the Tauri frontend
- **Tauri CLI** — `cargo install tauri-cli`

---

### macOS

#### Dependencies

```bash
brew install gmp readline
# Optional: MPFR arbitrary-precision reals
brew install mpfr
# Optional: FLINT (>= 3.0) — fast, rigorous algebraic-extension arithmetic + acb numerics
brew install flint
# Optional: GMP-ECM — advanced integer factorization
brew install gmp-ecm
# Optional: hardware-accelerated linear algebra (already in Xcode CLT as Accelerate)
# Optional: 2D graphics rendering for Plot[]/Show[]
brew install raylib
```

#### Build

```bash
# 1. Build the Mathilda sidecar and install it into the Tauri bundle
cd frontend
./build-sidecar.sh

# 2. Package the Tauri app
cargo tauri build
```

The `.app` bundle is in `frontend/src-tauri/target/release/bundle/macos/`.

---

### Linux (Ubuntu 22.04+ / Debian Bookworm+)

#### Dependencies

```bash
sudo apt-get update
sudo apt-get install -y \
  build-essential \
  libgmp-dev libreadline-dev \
  libwebkit2gtk-4.1-dev \
  libgtk-3-dev \
  libayatana-appindicator3-dev \
  librsvg2-dev \
  pkg-config \
  curl wget file

# Optional: MPFR for arbitrary-precision numerics
sudo apt-get install -y libmpfr-dev

# Optional: FLINT (>= 3.0) for fast, rigorous algebraic-extension arithmetic
sudo apt-get install -y libflint-dev    # >= 3.0 only on Ubuntu 24.04+/Debian Bookworm+
# On Ubuntu 22.04 (apt ships FLINT 2.x, which the >= 3.0 gate rejects), install
# 3.x from source instead — see ./tools/install-flint.sh and the FLINT note in
# the Troubleshooting section below. Without it Mathilda still builds (USE_FLINT=0).

# Optional: GMP-ECM for advanced integer factorization
sudo apt-get install -y libecm-dev

# Optional: graphics rendering for Plot[]/Show[]
sudo apt-get install -y libraylib-dev   # or build raylib from source

# Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.cargo/env
cargo install tauri-cli
```

> Note: Ubuntu 20.04 ships `libwebkit2gtk-4.0`, which is incompatible with
> Tauri v2. Ubuntu 22.04 or newer is required.

#### Build

```bash
cd frontend
chmod +x build-sidecar.sh
./build-sidecar.sh
cargo tauri build
```

The `.deb` and `.AppImage` are in `frontend/src-tauri/target/release/bundle/`.

---

### Windows (MSYS2 + MinGW64)

Mathilda uses GNU make and GCC. The recommended Windows build environment is
**MSYS2** with the MinGW64 toolchain. The Tauri side uses the standard Windows
Rust/MSVC toolchain.

#### Dependencies

1. Install [MSYS2](https://www.msys2.org/) (default path `C:\msys64`).

2. In the **MSYS2 MinGW64** shell, install the C toolchain, GMP and MPFR:
   ```bash
   pacman -S --noconfirm \
     mingw-w64-x86_64-gcc \
     mingw-w64-x86_64-gmp \
     mingw-w64-x86_64-mpfr \
     mingw-w64-x86_64-make \
     mingw-w64-x86_64-pkg-config
   ```
   MPFR is on by default (`USE_MPFR=1`) and, unlike the pkg-config-autodetected
   libraries, has no graceful-degrade — so it must be installed, or you must
   pass `USE_MPFR=0`, otherwise the link fails on `-lmpfr`. (On the **UCRT64**
   environment, substitute the `mingw-w64-ucrt-x86_64-` package prefix; this is
   what the CI Windows job uses.)

   **Supported environments: UCRT64 and MINGW64 (GCC).** The CLANG64 and CLANGARM64
   environments are deliberately unsupported: Mathilda is built with real GCC only,
   and the makefile stops with an error when `CC` is clang. The makefile never adds
   `/usr/include` or `/usr/local` on Windows, because under MSYS2 those paths hold the
   Cygwin (`msys2-runtime-devel`) headers, and a native build that read them would fail.

3. Install [Rust for Windows](https://rustup.rs/) using the
   `x86_64-pc-windows-msvc` host (the default Windows installer selects this).

4. Install [Node.js LTS](https://nodejs.org/) for Windows.

5. `cargo install tauri-cli`

6. Install [Visual Studio Build Tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/)
   (the MSVC linker is required by the Rust/MSVC toolchain).

7. Add the MSYS2 MinGW64 binaries to your Windows `PATH`:
   ```
   C:\msys64\mingw64\bin
   C:\msys64\usr\bin
   ```
   This lets `make`, `gcc`, and `pkg-config` work from PowerShell/CMD.

#### Build the C backend (PowerShell or CMD)

```powershell
# From the repo root — GNU make from MSYS2 must be on PATH
make USE_ECM=0 USE_READLINE=0 USE_GRAPHICS=0 USE_LAPACK=0 -j4
```

This produces `Mathilda.exe`. Readline is not available on Windows: the
`NO_READLINE` flag is set automatically when `BUILD_PLATFORM` is detected as
Windows, or explicitly by passing `USE_READLINE=0`. Pipe mode (used by the Tauri
frontend) does not require readline.

**Platform degradations on Windows.** Native MinGW has no POSIX signals, so
`TimeConstrained[expr, t]` runs *cooperative-only*: it enforces the budget
through the evaluator loop (and the deadline polls in heavy kernels), exactly
like the WSL-1 path, but cannot forcibly interrupt a single long-running C
builtin that never returns to the loop. POSIX threads are off (`USE_THREADS=0`),
so the parallel NDArray kernels run serially. These match the readline and
output-capture omissions — the kernel is otherwise fully functional. Timing
(`AbsoluteTiming`, `SessionTime`), `Pause`, and memory reporting (`MemoryInUse`,
`MaxMemoryUsed`) all work, via `QueryPerformanceCounter` / `Sleep` /
`GetProcessMemoryInfo` in `src/mtime.c` and `src/meminfo.c`.

#### Install the sidecar and package the app

```powershell
cd frontend
.\build-sidecar.bat
cargo tauri build
```

The `.msi` installer is in `frontend\src-tauri\target\release\bundle\msi\`.

---

### Dev Mode (all platforms)

```bash
cd frontend
./build-sidecar.sh       # macOS/Linux (use build-sidecar.bat on Windows)
npm install
cargo tauri dev
```

---

## Troubleshooting

**"sidecar not found" at Tauri startup**
Run `build-sidecar.sh` (or `build-sidecar.bat`) first and confirm that
`frontend/src-tauri/binaries/mathilda-<triple>` (or `.exe`) exists.

**GMP not found**
- macOS: `brew install gmp`
- Ubuntu/Debian: `sudo apt-get install libgmp-dev`
- Windows: install `mingw-w64-x86_64-gmp` via pacman in MSYS2 MinGW64

**Readline link error on Windows**
Build with `make USE_READLINE=0`. GNU Readline is not available in MinGW.
The Tauri frontend communicates with Mathilda via the NDJSON pipe protocol
which does not use readline, so interactive readline is never needed in
sidecar mode.

**webkit2gtk not found on Linux**
The required version is `libwebkit2gtk-4.1`. Ubuntu 20.04 only ships `4.0`.
Upgrade to Ubuntu 22.04 (Jammy) or Debian Bookworm.

**MPFR link error**
If you pass `USE_MPFR=1` (the default), MPFR must be installed. The makefile
links `-lmpfr` automatically when `USE_MPFR=1`. Omit it by passing
`USE_MPFR=0` (disables arbitrary-precision numerics with a runtime warning).

**`undefined reference to symbol 'cblas_ddot'` / `DSO missing from command line`**
Mathilda calls CBLAS (`cblas_*`) directly, and which library exports it varies:
reference `libblas` (Debian/Ubuntu), `libopenblas`, a separate `libcblas`, or
FlexiBLAS (Fedora/RHEL, whose netlib `libblas.so` carries no CBLAS). The makefile
link-probes `-lflexiblas`, `-lopenblas`, `-lcblas`, `-lblas` in turn and uses the
first that links. If your BLAS lives somewhere else, name it explicitly:
`make EXTRA_LIBS="-lmyblas"`, or build without it: `make USE_LAPACK=0`.

**"GMP-ECM not detected" although it is installed**
Fixed for GNU Make ≥ 4.3 (issue #93; the probe source used to start with a
literal `\#include`). If it persists, check that `ecm.h` and `libecm` are on the
default compiler paths (`gcc -lecm -lgmp` on a one-line program that includes
`<stdio.h>` and `<ecm.h>`), or pass them via `CFLAGS`/`LDFLAGS`.

**FLINT not detected / builds without the accelerated paths**
Mathilda needs **FLINT ≥ 3.0** — the release that merged ANTIC, Arb, and Calcium
into FLINT and added the generic-ring `gr` layer. Everything Mathilda accelerates
through FLINT (multivariate GCD/factoring over ℚ, number-field arithmetic over
ℚ(α), `RootReduce` via `qqbar`, rigorous `acb` numerics) depends on subsystems
that **did not exist in FLINT 2.x**, so 2.x cannot be used as a partial backend.
If `pkg-config --exists 'flint >= 3.0'` fails, the makefile prints a warning and
continues with `USE_FLINT=0`, using the classical (slower but still rigorous)
code — **no feature is lost, only speed**. Install a recent FLINT:

- **macOS:** `brew install flint`
- **Ubuntu 24.04+ / Debian Bookworm+:** `sudo apt-get install libflint-dev`
- **Ubuntu 22.04 / Debian Bullseye and older** ship FLINT **2.x** through apt,
  which the `>= 3.0` gate correctly rejects. Two easy ways to get 3.x without
  waiting for a distro upgrade:

  1. **Build from source** (Ubuntu 22.04 already ships GMP 6.2.1 + MPFR 4.1.0,
     which meet FLINT 3's minimums, so this is quick and clean):

     ```bash
     ./tools/install-flint.sh            # into /usr/local (uses sudo)
     # or a user-local install with no sudo:
     PREFIX="$HOME/.local" ./tools/install-flint.sh
     ```

     The script downloads a FLINT 3.x release, `./configure && make && make
     install`s it, refreshes `ldconfig`, and verifies `pkg-config` can see it.
     Override the version with `FLINT_VERSION=3.1.3 ./tools/install-flint.sh`.
     If you installed to a non-system prefix, export the `PKG_CONFIG_PATH` the
     script prints before rebuilding Mathilda.

  2. **conda-forge** (no compilation, no sudo):

     ```bash
     conda install -c conda-forge libflint pkg-config
     ```

     With the conda env active, `pkg-config` finds the conda FLINT automatically.

Either way, rebuild Mathilda with `make` and the FLINT paths light up on their
own. Confirm the version pkg-config sees with `pkg-config --modversion flint`.

**make: command not found (Windows)**
Ensure `C:\msys64\usr\bin` is on your Windows `PATH` so the MSYS2 `make`
utility is found from PowerShell/CMD.
