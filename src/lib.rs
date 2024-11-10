include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

#[cfg(test)]
mod tests {
    use std::{
        ffi::{CStr, CString},
        mem::MaybeUninit,
        ptr::null_mut,
    };

    use psqs::program::{mopac::Mopac, Program};

    use super::*;

    fn cleanup_mopac_output(basename: &str) -> std::io::Result<()> {
        let b = std::path::Path::new(basename);
        for ext in ["arc", "aux", "out"] {
            let p = b.with_extension(ext);
            if p.exists() {
                std::fs::remove_file(p)?;
            }
        }
        Ok(())
    }

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

        cleanup_mopac_output("testfiles/try").unwrap();
    }

    #[test]
    fn mopac_scf_works() {
        let mut atom = [6, 6, 6, 1, 1];
        #[rustfmt::skip]
        let mut coord = [
                0.000000000000, 0.003768239200, -1.686245109400,
                0.000000000000, 1.243805099800, 0.688097726900,
                0.000000000000, -1.242134139500, 0.700109341300,
                0.000000000000, 2.998368900600, 1.719180578800,
                0.000000000000, -3.003808100200, 1.718996398400,
        ];
        // this needs to be fully initialized
        let mut system = mopac_system {
            natom: 5,
            natom_move: 5,
            charge: 0,
            spin: 0,
            model: 3,     // PM6
            epsilon: 1.0, // no solvent
            atom: atom.as_mut_ptr(),
            coord: coord.as_mut_ptr(),
            nlattice: 0,
            nlattice_move: 0,
            pressure: 0.0,
            lattice: null_mut(),
            tolerance: 1e-8,
            max_time: 1,
        };

        let mut state = unsafe {
            let mut state = MaybeUninit::uninit();
            create_mopac_state(state.as_mut_ptr());
            state.assume_init()
        };

        // this can be uninitialized and should be allocated by mopac
        let mut props = MaybeUninit::uninit();
        unsafe {
            mopac_scf(&mut system, &mut state, props.as_mut_ptr());
            let props = props.assume_init();
            dbg!(&props);
        }
    }
}
