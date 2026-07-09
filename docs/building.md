# Building Mathilda

This document covers building the Mathilda C backend and the Tauri desktop
frontend on macOS, Linux, and Windows.

---

## Prerequisites (all platforms)

- **Rust toolchain** — install via [rustup.rs](https://rustup.rs/)
- **Node.js 18+** — for the Tauri frontend
- **Tauri CLI** — `cargo install tauri-cli`

---

## macOS

### Dependencies

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

### Build

```bash
# 1. Build the Mathilda sidecar and install it into the Tauri bundle
cd frontend
./build-sidecar.sh

# 2. Package the Tauri app
cargo tauri build
```

The `.app` bundle is in `frontend/src-tauri/target/release/bundle/macos/`.

---

## Linux (Ubuntu 22.04+ / Debian Bookworm+)

### Dependencies

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
sudo apt-get install -y libflint-dev    # needs >= 3.0 (Ubuntu 24.04+/Debian Bookworm+)

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

### Build

```bash
cd frontend
chmod +x build-sidecar.sh
./build-sidecar.sh
cargo tauri build
```

The `.deb` and `.AppImage` are in `frontend/src-tauri/target/release/bundle/`.

---

## Windows (MSYS2 + MinGW64)

Mathilda uses GNU make and GCC. The recommended Windows build environment is
**MSYS2** with the MinGW64 toolchain. The Tauri side uses the standard Windows
Rust/MSVC toolchain.

### Dependencies

1. Install [MSYS2](https://www.msys2.org/) (default path `C:\msys64`).

2. In the **MSYS2 MinGW64** shell, install the C toolchain and GMP:
   ```bash
   pacman -S --noconfirm \
     mingw-w64-x86_64-gcc \
     mingw-w64-x86_64-gmp \
     mingw-w64-x86_64-make \
     mingw-w64-x86_64-pkg-config
   ```

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

### Build the C backend (PowerShell or CMD)

```powershell
# From the repo root — GNU make from MSYS2 must be on PATH
make USE_ECM=0 USE_READLINE=0 USE_GRAPHICS=0 USE_LAPACK=0 -j4
```

Readline is not available on Windows. The `NO_READLINE` flag is set
automatically when `BUILD_PLATFORM` is detected as Windows, or explicitly by
passing `USE_READLINE=0`. Pipe mode (used by the Tauri frontend) does not
require readline.

### Install the sidecar and package the app

```powershell
cd frontend
.\build-sidecar.bat
cargo tauri build
```

The `.msi` installer is in `frontend\src-tauri\target\release\bundle\msi\`.

---

## C backend only (no Tauri frontend)

```bash
make USE_ECM=0 -j$(nproc)   # macOS/Linux
./Mathilda                   # interactive REPL
```

On Windows (MinGW64):
```bash
make USE_ECM=0 USE_READLINE=0 -j4
./Mathilda.exe
```

---

## Dev Mode (all platforms)

```bash
cd frontend
./build-sidecar.sh       # macOS/Linux (use build-sidecar.bat on Windows)
npm install
cargo tauri dev
```

---

## Build flags reference

| Flag | Default | Effect |
|------|---------|--------|
| `USE_ECM=0` | `1` | Disable GMP-ECM factorization. When `1`, the build autodetects the system GMP-ECM library (`gmp-ecm` / `libecm-dev`) via a compile-link probe and links `-lecm`; auto-disabled when it is not found. |
| `USE_MPFR=0` | `1` | Disable MPFR arbitrary-precision reals |
| `USE_FLINT=0` | `1` | Disable FLINT (≥ 3.0) kernels (falls back to the classical algebraic-extension GCD/Factor and non-FLINT numerics). Auto-disabled when FLINT ≥ 3.0 is not found by `pkg-config`. |
| `USE_READLINE=0` | `1` | Disable GNU Readline (required on Windows; auto-detected) |
| `USE_GRAPHICS=0` | `1` | Disable Raylib 2D graphics engine |
| `USE_LAPACK=0` | `1` | Disable BLAS/LAPACK linear-algebra kernels |

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

**FLINT not detected / builds without the accelerated paths**
Mathilda needs **FLINT ≥ 3.0** (the release that merged ANTIC for number-field
arithmetic). If `pkg-config --exists 'flint >= 3.0'` fails, the makefile prints a
warning and continues with `USE_FLINT=0`, using the classical (slower but still
rigorous) algebraic-extension code. Install a recent FLINT:
- macOS: `brew install flint`
- Ubuntu 24.04+/Debian Bookworm+: `sudo apt-get install libflint-dev`
- Older distros ship FLINT 2.x, which lacks ANTIC — build FLINT ≥ 3.0 from source
  or use a backport. Confirm with `pkg-config --modversion flint`.

**make: command not found (Windows)**
Ensure `C:\msys64\usr\bin` is on your Windows `PATH` so the MSYS2 `make`
utility is found from PowerShell/CMD.
