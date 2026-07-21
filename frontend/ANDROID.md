# Running Mathilda Notebook on Android

Like iOS, Android forbids spawning the `mathilda` sidecar process, so the CAS
kernel is compiled to a static library and linked **in-process** via the same
C FFI (`src/ffi/mathilda_ffi.{c,h}`) and the shared `#[cfg(mobile)]` Rust
backend (`kernel_ffi.rs`). This document is the Android runbook.

## Architecture

```
Desktop:  Svelte UI ─IPC─▶ Rust (kernel.rs)     ─spawn/stdio─▶ mathilda sidecar
Mobile:   Svelte UI ─IPC─▶ Rust (kernel_ffi.rs) ────FFI──────▶ libmathilda.a (in-process)
```

Android-specific pieces:

| Piece | File |
|-------|------|
| NDK cross-compile of GMP + MPFR + kernel | `build-android-lib.sh` |
| Android link flags (static libs) | `src-tauri/build.rs` (`target_os = "android"`) |
| Drop sidecar for Android | `src-tauri/tauri.android.conf.json` |
| `init.m` bootstrap on Android | embedded in the binary + extracted at startup (see below) |

### init.m on Android

Android app resources live inside the APK and are **not** reachable via the C
loader's `fopen()`. So the `internal/` module tree (init.m + derivative/integral
tables, ~160 KB) is **embedded in the binary** with `include_dir!` and extracted
once to the app's data dir at startup (`kernel_ffi.rs`, `#[cfg(target_os =
"android")]`); `MATHILDA_HOME` then points at that real path. iOS bundles it as
a normal app resource instead.

## Prerequisites

- **Android SDK + NDK** (r27 tested). Set `ANDROID_HOME` and `NDK_HOME` (or let
  the script auto-detect `$ANDROID_HOME/ndk/<ver>`).
- **JDK 17** for Gradle. The generated project ships Gradle 8.14, which does
  **not** support JDK 25 (`Unsupported class file major version 69`). Point
  `JAVA_HOME` at a 17/21 JDK:
  ```bash
  export JAVA_HOME=/opt/homebrew/opt/openjdk@17/libexec/openjdk.jdk/Contents/Home
  ```
- **Rust targets**: `rustup target add aarch64-linux-android x86_64-linux-android`
- GMP + MPFR source tarballs (same as iOS).

## Build & run in the emulator

Verified on: NDK 27.2, Gradle 8.14 + JDK 17, Pixel 8 AVD (API 34, arm64-v8a).

```bash
cd frontend
npm install
export ANDROID_HOME=~/Library/Android/sdk
export NDK_HOME=$ANDROID_HOME/ndk/27.2.12479018
export JAVA_HOME=/opt/homebrew/opt/openjdk@17/libexec/openjdk.jdk/Contents/Home
export PATH="$JAVA_HOME/bin:$ANDROID_HOME/platform-tools:$ANDROID_HOME/emulator:$PATH"

# 1. Cross-compile GMP + MPFR + libmathilda.a for arm64-v8a (API 26; the kernel
#    uses complex.h clog/cpow/csqrt, gated __INTRODUCED_IN(26) by Bionic).
GMP_TARBALL=/path/to/gmp-6.3.0.tar.xz \
MPFR_TARBALL=/path/to/mpfr-4.2.1.tar.xz \
  ./build-android-lib.sh arm64
# → src-tauri/gen/android-libs/aarch64-linux-android/{libmathilda,libgmp,libmpfr}.a

# 2. One-time: generate the Gradle project, then set minSdk = 26.
npm run android:init
sed -i '' 's/minSdk = 24/minSdk = 26/' src-tauri/gen/android/app/build.gradle.kts

# 3. Boot an emulator (API >= 26, arm64-v8a on Apple Silicon) and run.
emulator -avd Pixel_8 &
npm run android:dev            # builds, installs, launches with live reload

# Or build a debug APK and install it manually:
npx tauri android build --debug --target aarch64 --apk
adb install -r src-tauri/gen/android/app/build/outputs/apk/universal/debug/app-universal-debug.apk
adb shell monkey -p com.mathilda.notebook -c android.intent.category.LAUNCHER 1
```

`x86_64` emulators (Intel hosts): `./build-android-lib.sh x86_64` and
`--target x86_64`.

## Knobs

- `MATHILDA_ANDROID_LIBDIR` — override the lib search dir (`build.rs` default is
  `gen/android-libs/<target-triple>`).
- `MATHILDA_ANDROID_WITH_MPFR=0` — drop the MPFR link (only once the `USE_MPFR=0`
  minimal kernel build compiles upstream).
- `API=… ABI=… ./build-android-lib.sh` — override NDK API level / ABI.

## Known limitations

- **minSdk 26** (Android 8.0) because of Bionic's complex-math gating. Covers
  ~99% of active devices.
- No graphics / FLINT / ECM / LAPACK (built with those subsystems off).
- `interrupt` / `restart` are in-process no-ops (shared with iOS): a synchronous
  C evaluation can't be cancelled from another thread, and the symbol table has
  no reset entry point yet.
- Release/Play Store builds require a signing keystore (`tauri android build`
  without `--debug`).
- Release builds targeting API 35+ also need 16 KB page alignment
  (`-Wl,-z,max-page-size=16384` on the cdylib link); debug builds only warn.
- The webview uses `user-scalable=no` (so OS pinch-zoom doesn't fight the
  canvas's own pinch-zoom), which disables system accessibility zoom inside the
  app — an intentional trade-off worth revisiting if a11y zoom is needed.
