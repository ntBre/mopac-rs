use std::ffi::CString;

use approx::assert_abs_diff_eq;
use psqs::program::{
    mopac::{Mopac, KCALHT},
    Program,
};
use tempfile::NamedTempFile;

use mopac_rs::run_mopac_from_input;

use common::C3H2_ENERGY_KCAL_MOL;

mod common;

#[test]
fn main() {
    let f = NamedTempFile::with_suffix(".mop").unwrap();
    std::fs::copy("testfiles/try.mop", &f).unwrap();

    let path_name = f.path().to_str().unwrap();
    let name = CString::new(path_name).unwrap();
    let ptr = name.into_raw();
    unsafe {
        run_mopac_from_input(ptr);
        drop(CString::from_raw(ptr));
    }

    let got = Mopac::read_output(f.path().with_extension("").to_str().unwrap())
        .unwrap();

    assert!(got.cart_geom.is_some());
    assert_eq!(got.cart_geom.as_ref().unwrap().len(), 5);
    assert_abs_diff_eq!(
        got.energy * KCALHT,
        C3H2_ENERGY_KCAL_MOL,
        epsilon = 1e-8
    );
}
