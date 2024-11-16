use std::path::PathBuf;

fn main() {
    let Ok(mopac_path) = std::env::var("MOPAC_PATH") else {
        eprintln!("MOPAC_PATH must be set to dir containing build/libmopac.so");
        std::process::exit(1);
    };

    println!("cargo:rerun-if-env-changed=MOPAC_PATH");

    println!("cargo:rustc-link-search={mopac_path}/build");
    println!("cargo:rustc-link-lib=mopac");

    println!("cargo:rustc-link-arg=-Wl,-rpath,{mopac_path}/build");

    let bindings = bindgen::Builder::default()
        .header("wrapper.h")
        .clang_arg(format!("-I{mopac_path}/include"))
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .unwrap();

    let out_path = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .unwrap();
}
