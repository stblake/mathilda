fn main() {
    // On iOS the Mathilda kernel is linked in-process as a static library
    // (there is no sidecar). The cross-compiled archives are produced by
    // `frontend/build-ios-lib.sh`; point this build at their output directory
    // with MATHILDA_IOS_LIBDIR (defaults to the conventional location).
    //
    // The kernel is built with USE_MPFR=1 by default (build-ios-lib.sh), so
    // libmpfr.a is linked unless MATHILDA_IOS_WITH_MPFR=0 is set explicitly.
    let target_os = std::env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    if target_os == "ios" {
        // Default to the conventional location build-ios-lib.sh writes to:
        //   <src-tauri>/gen/ios-libs/<target-triple>/lib*.a
        // Use an ABSOLUTE path (via CARGO_MANIFEST_DIR) because `tauri ios
        // build` invokes cargo from an Xcode script phase whose cwd and
        // environment differ — a relative path or an exported env var is not
        // reliably visible there. MATHILDA_IOS_LIBDIR still overrides.
        let libdir = std::env::var("MATHILDA_IOS_LIBDIR").unwrap_or_else(|_| {
            let manifest = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_default();
            let target = std::env::var("TARGET").unwrap_or_default();
            format!("{manifest}/gen/ios-libs/{target}")
        });
        println!("cargo:rerun-if-env-changed=MATHILDA_IOS_LIBDIR");
        println!("cargo:rerun-if-env-changed=MATHILDA_IOS_WITH_MPFR");
        println!("cargo:rustc-link-search=native={libdir}");
        println!("cargo:rustc-link-lib=static=mathilda");
        println!("cargo:rustc-link-lib=static=gmp");
        if std::env::var("MATHILDA_IOS_WITH_MPFR").as_deref() != Ok("0") {
            println!("cargo:rustc-link-lib=static=mpfr");
        }
        // Accelerate provides BLAS/LAPACK on Apple platforms (harmless if the
        // kernel was built without USE_LAPACK — no symbols are referenced).
        println!("cargo:rustc-link-lib=framework=Accelerate");
    } else if target_os == "android" {
        // Android also links the kernel in-process (no sidecar). Archives are
        // produced per-ABI by `frontend/build-android-lib.sh` into
        //   <src-tauri>/gen/android-libs/<target-triple>/lib*.a
        // Absolute path via CARGO_MANIFEST_DIR: the Gradle/cargo integration
        // runs from its own cwd, so a relative path is unreliable.
        let libdir = std::env::var("MATHILDA_ANDROID_LIBDIR").unwrap_or_else(|_| {
            let manifest = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_default();
            let target = std::env::var("TARGET").unwrap_or_default();
            format!("{manifest}/gen/android-libs/{target}")
        });
        println!("cargo:rerun-if-env-changed=MATHILDA_ANDROID_LIBDIR");
        println!("cargo:rerun-if-env-changed=MATHILDA_ANDROID_WITH_MPFR");
        println!("cargo:rustc-link-search=native={libdir}");
        println!("cargo:rustc-link-lib=static=mathilda");
        println!("cargo:rustc-link-lib=static=gmp");
        if std::env::var("MATHILDA_ANDROID_WITH_MPFR").as_deref() != Ok("0") {
            println!("cargo:rustc-link-lib=static=mpfr");
        }
        // libm/libc come from the NDK sysroot automatically; no extra libs.
    }

    tauri_build::build()
}
