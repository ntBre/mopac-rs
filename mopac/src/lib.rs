use std::{ffi::CStr, mem::MaybeUninit};

use mopac_sys::{
    create_mopac_state, mopac_properties, mopac_relax, mopac_scf, mopac_state,
    mopac_system,
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
        let mut state = State::default();

        let props = unsafe {
            let mut props = MaybeUninit::uninit();
            mopac_scf(&mut self.system, &mut state.0, props.as_mut_ptr());
            Properties {
                inner: props.assume_init(),
                natoms: self.system.natom as usize,
            }
        };

        props.check_errors()?;

        Ok(props)
    }

    pub fn optimize(&mut self) -> Result<Properties, Vec<String>> {
        let mut state = State::default();

        let props = unsafe {
            let mut props = MaybeUninit::uninit();
            mopac_relax(&mut self.system, &mut state.0, props.as_mut_ptr());
            Properties {
                inner: props.assume_init(),
                natoms: self.system.natom as usize,
            }
        };

        props.check_errors()?;

        Ok(props)
    }
}

pub struct Properties {
    natoms: usize,
    inner: mopac_properties,
}

impl Properties {
    /// Return the final heat of formation in kcal/mol
    pub fn final_energy(&self) -> f64 {
        self.inner.heat
    }

    pub fn coordinates(&self) -> &[f64] {
        // Safety: `self` can only be constructed by a successful return of one
        // of the [State::scf] or [State::optimize] calls, which both pass
        // [mopac_system::natom] as [mopac_system::natom_move]. According to the
        // API docs, `coord_update` should thus be initialized with the proper
        // length.
        unsafe {
            &*std::ptr::slice_from_raw_parts(
                self.inner.coord_update,
                3 * self.natoms,
            )
        }
    }

    fn check_errors(&self) -> Result<(), Vec<String>> {
        if self.inner.nerror > 0 {
            let mut ret = Vec::new();
            for i in 0..self.inner.nerror {
                let s = unsafe {
                    CStr::from_ptr(*self.inner.error_msg.offset(i as isize))
                };
                ret.push(s.to_string_lossy().to_string());
            }
            return Err(ret);
        }

        Ok(())
    }
}

impl Drop for Properties {
    fn drop(&mut self) {
        unsafe { mopac_sys::destroy_mopac_properties(&mut self.inner) }
    }
}

pub struct State(mopac_state);

impl Default for State {
    fn default() -> Self {
        let mut state = unsafe {
            let mut state = MaybeUninit::uninit();
            create_mopac_state(state.as_mut_ptr());
            state.assume_init()
        };
        state.mpack = 0;
        Self(state)
    }
}

impl Drop for State {
    fn drop(&mut self) {
        unsafe { mopac_sys::destroy_mopac_state(&mut self.0) }
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;
    use symm::molecule;

    use super::*;

    fn c3h2() -> Molecule {
        let mut mol = molecule! {
            C      0.000000000000      0.003768239200     -1.686245109400
            C      0.000000000000      1.243805099800      0.688097726900
            C      0.000000000000     -1.242134139500      0.700109341300
            H      0.000000000000      2.998368900600      1.719180578800
            H      0.000000000000     -3.003808100200      1.718996398400
        };
        mol.to_angstrom();
        mol
    }

    #[test]
    fn c3h2_scf() {
        let mut system = System::new(c3h2(), 0, 0);
        let props = system.scf().unwrap();

        assert_abs_diff_eq!(props.final_energy(), 128.29901260, epsilon = 1e-8);
    }

    #[test]
    fn c3h2_opt() {
        let mut system = System::new(c3h2(), 0, 0);

        let props = system.optimize().unwrap();
        let coords_init = &system._coordinates[..];
        let coords_final = props.coordinates();

        // Make sure some optimization occurred
        assert_ne!(coords_init, coords_final);

        // but not too much
        assert_abs_diff_eq!(coords_init, props.coordinates(), epsilon = 1e-1);

        // Check the final energy, note the difference from the unoptimized
        // geometry above
        assert_abs_diff_eq!(props.final_energy(), 126.60240431, epsilon = 1e-8);
    }
}
