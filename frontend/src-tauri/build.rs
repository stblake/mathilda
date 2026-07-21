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
        let libdir = std::env::var("MATHILDA_IOS_LIBDIR")
            .unwrap_or_else(|_| "gen/ios-libs".to_string());
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
    }

    tauri_build::build()
}
