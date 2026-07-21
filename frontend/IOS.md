# Running Mathilda Notebook on iOS

The desktop notebook runs the CAS as a **sidecar process** (spawns the
`mathilda` binary, talks NDJSON over stdio). **iOS forbids spawning child
processes**, so on iOS the kernel is compiled to a static library and linked
**in-process** through a C FFI. This document is the runbook.

## Architecture

```
Desktop:  Svelte UI ──IPC──▶ Rust (kernel.rs) ──spawn/stdio──▶ mathilda sidecar
iOS:      Svelte UI ──IPC──▶ Rust (kernel_ffi.rs) ──FFI──▶ libmathilda.a (in-process)
```

| Piece | File |
|-------|------|
| C ABI (init / eval / eval_latex / set_home / free) | `../src/ffi/mathilda_ffi.{c,h}` |
| Static archive target | `../makefile` → `libmathilda.a` |
| Safe Rust bindings | `src-tauri/src/ffi.rs` |
| In-process kernel (mobile) | `src-tauri/src/kernel_ffi.rs` |
| Backend selection (`#[cfg(mobile)]`) | `src-tauri/src/lib.rs` |
| iOS link flags | `src-tauri/build.rs` |
| GMP/MPFR + kernel cross-compile | `build-ios-lib.sh` |

The kernel is **not reentrant** (process-global symbol table); `kernel_ffi.rs`
serializes every evaluation behind a mutex and runs it on a blocking thread.

## Prerequisites (one-time, requires human action)

Command Line Tools are **not** enough — you need full Xcode:

```bash
# Install Xcode from the Mac App Store, then:
sudo xcode-select -s /Applications/Xcode.app
sudo xcodebuild -license accept
xcodebuild -downloadPlatform iOS          # a Simulator runtime

# Rust iOS targets:
rustup target add aarch64-apple-ios aarch64-apple-ios-sim

# Verify the SDK is visible (should print a path, not an error):
xcrun --sdk iphonesimulator --show-sdk-path
```

Also download GMP and MPFR sources (neither cross-compiles via Homebrew):

```bash
curl -LO https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
curl -LO https://www.mpfr.org/mpfr-current/mpfr-4.2.1.tar.xz
```

## Build & run in the Simulator

```bash
cd frontend

# 1. Cross-compile GMP + MPFR + libmathilda.a for the Simulator (arm64).
GMP_TARBALL=/path/to/gmp-6.3.0.tar.xz \
MPFR_TARBALL=/path/to/mpfr-4.2.1.tar.xz \
  ./build-ios-lib.sh sim
# → src-tauri/gen/ios-libs/aarch64-apple-ios-sim/{libmathilda,libgmp,libmpfr}.a

# 2. One-time: generate the Xcode project.
npm install
npm run ios:init

# 3. Bundle the internal module tree (init.m + tables) as an app resource so
#    the kernel finds its bootstrap on-device. Point tauri.conf.json bundle
#    resources at ../../src/internal → "internal" (see "Bundling init.m" below),
#    or copy it into src-tauri/resources/internal before building.

# 4. Run in the Simulator.
MATHILDA_IOS_LIBDIR="$(pwd)/src-tauri/gen/ios-libs/aarch64-apple-ios-sim" \
  npm run ios:dev
```

For a device / TestFlight build use `./build-ios-lib.sh device` and
`npm run ios:build` (device builds require an Apple signing identity).

## Bundling init.m

`kernel_ffi.rs` sets `MATHILDA_HOME` to `<app resources>/internal` before
initializing, and the loader resolves `init.m` from there. Ensure the
`src/internal` tree is bundled — e.g. in `tauri.conf.json`:

```jsonc
"bundle": {
  "resources": { "../../src/internal": "internal" }
}
```

(Keep this out of the desktop bundle config if you don't want the tree shipped
twice; the desktop sidecar locates `init.m` on its own.)

## Capabilities & knobs

- `MATHILDA_IOS_LIBDIR` — dir holding the cross-built `.a` files (consumed by
  `build.rs`). Defaults to `gen/ios-libs`.
- `MATHILDA_IOS_WITH_MPFR=0` — drop the MPFR link (only valid once the
  `USE_MPFR=0` minimal kernel build compiles; see below).
- `WITH_MPFR=0 ./build-ios-lib.sh` — build the kernel without MPFR (same caveat).

## Known limitations

- **Full Xcode required** — nothing iOS compiles or runs without the SDK.
- **`USE_MPFR=0` does not yet build** — several subsystems (linalg eigen,
  numerical roots, Gröbner, some special functions) reference MPFR
  unconditionally. iOS therefore ships GMP+MPFR. Fixing the guards is tracked
  as follow-up and would slim the binary.
- **No graphics / FLINT / ECM / LAPACK on iOS** — the mobile kernel is built
  with those subsystems off (text-placeholder graphics). BLAS/LAPACK via
  Accelerate could be re-enabled later (the framework is linked).
- **`interrupt` / `restart` are no-ops** in-process: a synchronous C evaluation
  can't be cancelled from another thread, and the symbol table has no reset
  entry point yet. A `mathilda_ffi_reset()` would enable true "Run All".
