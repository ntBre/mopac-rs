include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

#[cfg(test)]
mod tests {
    use std::ffi::CString;

    use super::*;

    #[test]
    fn run_mopac_from_input_works() {
        let name = CString::new("testfiles/try.mop").unwrap();
        let ptr = name.into_raw();
        unsafe {
            run_mopac_from_input(ptr);
            drop(CString::from_raw(ptr));
        }
    }
}
