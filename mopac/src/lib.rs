//! A safe Rust interface to the [MOPAC](https://github.com/openmopac/mopac)
//! chemistry package.
//!
//! This package is based on the C API first released in MOPAC
//! [23.0.0](https://github.com/openmopac/mopac/releases/tag/v23.0.0), which
//! contains a simple diskless API for running MOPAC calculations. As the
//! "diskless" part implies, instead of using the normal MOPAC input file,
//! callers must set up a [mopac_system] in code and invoke one of the
//! [mopac_scf], [mopac_relax], or [mopac_vibe] functions to obtain the
//! corresponding [mopac_properties]. Calling these functions directly is
//! unsafe, so this package attempts to wrap the unsafety in a safe, idiomatic
//! Rust interface.
//!
//! # Examples
//!
//! Each of the examples below starts with construction of an input [Molecule].
//! A [Molecule] is just a sequence of [Atom]s. One of the easiest ways to
//! construct one is with the [molecule] macro:
//!
//! ```
//! use mopac::molecule;
//!
//! let mol = molecule! {
//!     C      0.000000000000      0.003768239200     -1.686245109400
//!     C      0.000000000000      1.243805099800      0.688097726900
//!     C      0.000000000000     -1.242134139500      0.700109341300
//!     H      0.000000000000      2.998368900600      1.719180578800
//!     H      0.000000000000     -3.003808100200      1.718996398400
//! };
//! ```
//!
//! With a [Molecule] in hand, you can construct a [System]. You also need to
//! specify the `charge` and `spin`, which are often 0:
//!
//! ```
//! # use mopac::{molecule, System};
//! # let mol = molecule! {
//! #    C      0.000000000000      0.003768239200     -1.686245109400
//! #    C      0.000000000000      1.243805099800      0.688097726900
//! #    C      0.000000000000     -1.242134139500      0.700109341300
//! #    H      0.000000000000      2.998368900600      1.719180578800
//! #    H      0.000000000000     -3.003808100200      1.718996398400
//! # };
//! let mut system = System::new(mol, 0, 0);
//! ```
//!
//! Note that `mut` is necessary for the method calls coming in the next couple
//! of sections.
//!
//! ## Single-point energies
//!
//! From the system, you can compute a single-point energy with the
//! [System::scf] method and return the final heat of formation in kcal/mol
//! with [Properties::final_energy]:
//!
//! ```
//! # use mopac::{molecule, System};
//! # let mol = molecule! {
//! #    C      0.000000000000      0.001993398537     -0.892023662873
//! #    C      0.000000000000      0.657972897794      0.364003697530
//! #    C      0.000000000000     -0.657088959796      0.370357841548
//! #    H      0.000000000000      1.586137148417      0.909446526185
//! #    H      0.000000000000     -1.589014485006      0.909349094754
//! # };
//! # let mut system = System::new(mol, 0, 0);
//! let props = system.scf().unwrap();
//! assert!((props.final_energy() - 128.33088672).abs() < 1e-8);
//! ```
//!
//! ## Geometry optimization
//!
//! Similarly, you can run a geometry optimization with [System::optimize] and
//! retrieve the results with [Properties::coordinates]:
//!
//! ```
//! use mopac::{molecule, System};
//!
//! let mol = molecule! {
//!    C      0.000000000000      0.001993398537     -0.892023662873
//!    C      0.000000000000      0.657972897794      0.364003697530
//!    C      0.000000000000     -0.657088959796      0.370357841548
//!    H      0.000000000000      1.586137148417      0.909446526185
//!    H      0.000000000000     -1.589014485006      0.909349094754
//! };
//! let mut system = System::new(mol, 0, 0);
//! let props = system.optimize().unwrap();
//! assert_ne!(props.coordinates(), system.coordinates());
//! ```
//!
//! ## Harmonic frequencies
//!
//! Finally, you can compute the harmonic vibrational frequencies with
//! [System::frequencies] and retrieve the results with
//! [Properties::frequencies]:
//!
//! ```
//! use mopac::{molecule, System};
//!
//! let mol = molecule! {
//!    C      0.000000000000      0.001993398537     -0.892023662873
//!    C      0.000000000000      0.657972897794      0.364003697530
//!    C      0.000000000000     -0.657088959796      0.370357841548
//!    H      0.000000000000      1.586137148417      0.909446526185
//!    H      0.000000000000     -1.589014485006      0.909349094754
//! };
//! // SCF calculations often fail at the default 1e-8 tolerance for some reason
//! let mut system = System::new(mol, 0, 0).tolerance(1e-4);
//! let props = system.frequencies().unwrap();
//! let freqs = props.frequencies().to_vec();
//! let mut freqs: Vec<_> = freqs
//!     .into_iter()
//!     .filter_map(|f| if f > 10.0 { Some(f.round()) } else { None })
//!     .collect();
//! freqs.reverse();
//! assert_eq!(
//!     freqs,
//!     [2774.0, 2739.0, 1954.0, 1291.0, 1090.0, 1011.0, 994.0, 989.0, 970.0],
//! );
//! ```

use std::{ffi::CStr, mem::MaybeUninit};

use mopac_sys::{
    create_mopac_state, mopac_properties, mopac_relax, mopac_scf, mopac_state,
    mopac_system, mopac_vibe,
};

pub use symm::{molecule, Atom, Molecule};

pub struct System {
    /// These slices must outlive `system` because `system` contains raw
    /// pointers to them
    _atomic_numbers: Box<[i32]>,
    coordinates: Box<[f64]>,
    system: mopac_system,
}

impl System {
    /// Construct a MOPAC system from `mol`. The coordinates in `mol` must be in
    /// Ångstroms.
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
            coordinates: coord,
        }
    }

    pub fn tolerance(mut self, tolerance: f64) -> Self {
        self.system.tolerance = tolerance;
        self
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

    /// Compute the harmonic frequencies for `self`.
    pub fn frequencies(&mut self) -> Result<Properties, Vec<String>> {
        let mut state = State::default();

        let props = unsafe {
            let mut props = MaybeUninit::uninit();
            mopac_vibe(&mut self.system, &mut state.0, props.as_mut_ptr());
            Properties {
                inner: props.assume_init(),
                natoms: self.system.natom as usize,
            }
        };

        props.check_errors()?;

        Ok(props)
    }

    pub fn coordinates(&self) -> &[f64] {
        &self.coordinates
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
        assert!(!self.inner.coord_update.is_null());
        unsafe {
            &*std::ptr::slice_from_raw_parts(
                self.inner.coord_update,
                3 * self.natoms,
            )
        }
    }

    pub fn frequencies(&self) -> &[f64] {
        assert!(!self.inner.freq.is_null());
        unsafe {
            &*std::ptr::slice_from_raw_parts(self.inner.freq, 3 * self.natoms)
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

#[cfg_attr(test, derive(Debug))]
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
        let coords_init = &system.coordinates[..];
        let coords_final = props.coordinates();

        // Make sure some optimization occurred
        assert_ne!(coords_init, coords_final);

        // but not too much
        assert_abs_diff_eq!(coords_init, props.coordinates(), epsilon = 1e-1);

        // Check the final energy, note the difference from the unoptimized
        // geometry above
        assert_abs_diff_eq!(props.final_energy(), 126.60240432, epsilon = 1e-8);
    }

    #[test]
    fn c3h2_freqs() {
        let mut system = System::new(c3h2(), 0, 0);
        // for some reason scf won't converge with the default 1e-8
        system.system.tolerance = 1e-4;
        let props = system.frequencies().unwrap();

        let want = [
            -233.0134837743131,
            -202.55035772854717,
            -166.25090736634493,
            0.029256741968343082,
            0.03073549436336908,
            0.03290405665074747,
            970.8366647056642,
            988.6960714808775,
            995.2425137409377,
            1011.892959945535,
            1089.2667229480483,
            1289.7297315059488,
            1951.9313829222758,
            2737.278109703735,
            2772.1784613501504,
        ];

        assert_abs_diff_eq!(&want[..], props.frequencies(), epsilon = 1e-1);
    }
}
