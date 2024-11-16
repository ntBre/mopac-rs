use std::mem::MaybeUninit;

use approx::assert_abs_diff_eq;
use mopac_rs::{create_mopac_state, mopac_scf, mopac_system};

use common::C3H2_ENERGY_KCAL_MOL;

mod common;

#[test]
fn main() {
    let mut atom = [6, 6, 6, 1, 1];
    #[rustfmt::skip]
    let coord = [
        0.000000000000, 0.003768239200, -1.686245109400,
        0.000000000000, 1.243805099800, 0.688097726900,
        0.000000000000, -1.242134139500, 0.700109341300,
        0.000000000000, 2.998368900600, 1.719180578800,
        0.000000000000, -3.003808100200, 1.718996398400,
    ];
    let mut coord = coord.map(|x| x * 0.529_177_210_9);
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
        lattice: std::ptr::null_mut(),
        tolerance: 1e-8,
        max_time: 1,
    };

    let mut state = unsafe {
        let mut state = MaybeUninit::uninit();
        create_mopac_state(state.as_mut_ptr());
        state.assume_init()
    };

    // this can be uninitialized and should be allocated by mopac
    let props = unsafe {
        let mut props = MaybeUninit::uninit();
        mopac_scf(&mut system, &mut state, props.as_mut_ptr());
        props.assume_init()
    };

    assert_eq!(props.nerror, 0);
    assert_abs_diff_eq!(props.heat, C3H2_ENERGY_KCAL_MOL, epsilon = 1e-8);
}
