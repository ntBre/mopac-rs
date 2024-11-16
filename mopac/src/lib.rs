use std::{ffi::CStr, mem::MaybeUninit};

use mopac_sys::{
    create_mopac_state, mopac_properties, mopac_scf, mopac_system,
};

pub use symm::Molecule;

pub struct System {
    /// These slices must outlive `system` because `system` contains raw
    /// pointers to them
    _atomic_numbers: Box<[i32]>,
    _coordinates: Box<[f64]>,
    system: mopac_system,
}

impl System {
    pub fn new(mol: Molecule, charge: i32, spin: i32) -> Self {
        let natoms = mol.atoms.len();
        let mut atom = {
            let v: Vec<i32> =
                mol.atomic_numbers().into_iter().map(|n| n as i32).collect();
            v.into_boxed_slice()
        };
        let mut coord = {
            let v: Vec<_> = mol.atoms.iter().flat_map(|a| a.coord()).collect();
            v.into_boxed_slice()
        };
        Self {
            system: mopac_system {
                natom: natoms as _,
                natom_move: natoms as _,
                charge,
                spin,
                model: 3,     // PM6
                epsilon: 1.0, // no solvent
                atom: atom.as_mut_ptr(),
                coord: coord.as_mut_ptr(),
                nlattice: 0,
                nlattice_move: 0,
                pressure: 0.0,
                lattice: std::ptr::null_mut(),
                tolerance: 1e-8,
                max_time: 172_800, // default time limit from the docs
            },
            _atomic_numbers: atom,
            _coordinates: coord,
        }
    }

    pub fn scf(&mut self) -> Result<Properties, Vec<String>> {
        let mut state = unsafe {
            let mut state = MaybeUninit::uninit();
            create_mopac_state(state.as_mut_ptr());
            let mut state = state.assume_init();
            state.mpack = 0;
            state
        };

        let props = unsafe {
            let mut props = MaybeUninit::uninit();
            mopac_scf(&mut self.system, &mut state, props.as_mut_ptr());
            props.assume_init()
        };

        if props.nerror > 0 {
            let mut ret = Vec::new();
            for i in 0..props.nerror {
                let s = unsafe {
                    CStr::from_ptr(*props.error_msg.offset(i as isize))
                };
                ret.push(s.to_string_lossy().to_string());
            }
            return Err(ret);
        }

        Ok(Properties(props))
    }
}

pub struct Properties(mopac_properties);

impl Properties {
    /// Return the final heat of formation in kcal/mol
    pub fn final_energy(&self) -> f64 {
        self.0.heat
    }
}

impl Drop for Properties {
    fn drop(&mut self) {
        unsafe { mopac_sys::destroy_mopac_properties(&mut self.0) }
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;
    use symm::molecule;

    use super::*;

    #[test]
    fn c3h2_scf() {
        let mut mol = molecule! {
            C      0.000000000000      0.003768239200     -1.686245109400
            C      0.000000000000      1.243805099800      0.688097726900
            C      0.000000000000     -1.242134139500      0.700109341300
            H      0.000000000000      2.998368900600      1.719180578800
            H      0.000000000000     -3.003808100200      1.718996398400
        };
        mol.to_angstrom();
        let mut system = System::new(mol, 0, 0);
        let props = system.scf().unwrap();

        assert_abs_diff_eq!(props.final_energy(), 128.29901260, epsilon = 1e-8);
    }
}
