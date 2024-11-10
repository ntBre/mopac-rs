include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

#[cfg(test)]
mod tests {
    use std::ffi::CString;

    use psqs::program::{mopac::Mopac, Program};

    use super::*;

    #[test]
    fn run_mopac_from_input_works() {
        let name = CString::new("testfiles/try.mop").unwrap();
        let ptr = name.into_raw();
        unsafe {
            run_mopac_from_input(ptr);
            drop(CString::from_raw(ptr));
        }

        let got = Mopac::read_output("testfiles/try").unwrap();

        assert!(got.cart_geom.is_some());
        assert_eq!(got.cart_geom.as_ref().unwrap().len(), 5);
        assert_eq!(got.energy, 0.204457543878102);
    }
}
