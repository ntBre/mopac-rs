fn main() {
    let Ok(mopac_path) = std::env::var("MOPAC_PATH") else {
        eprintln!("MOPAC_PATH must be set to dir containing build/libmopac.so");
        std::process::exit(1);
    };

    println!("cargo:rerun-if-env-changed=MOPAC_PATH");

    println!("cargo:rustc-link-arg=-Wl,-rpath,{mopac_path}/build");
}
