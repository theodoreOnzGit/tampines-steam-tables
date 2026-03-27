use uom::ConstZero;
use uom::si::f64::*;
use uom::si::force::newton;
use uom::si::pressure::{bar, millibar, pascal};
use uom::si::ratio::ratio;
use uom::si::volume::cubic_meter;

use crate::prelude::functional_programming::hs_flash_eqm::v_hs_eqm;
use crate::prelude::functional_programming::ph_flash_eqm::{s_ph_eqm, w_ph_eqm};
use crate::prelude::functional_programming::ps_flash_eqm::{h_ps_eqm, v_ps_eqm};
use crate::prelude::{TampinesSteamTableCV};


// for isentropic/adiabatic diffusers/nozzles, we follow the mach number
// formulae 
//
pub fn get_dp_isentropic_nozzle_diffuser(
    a1: Area,
    a2: Area,
    p1: Pressure,
    h1: AvailableEnergy,
    v1: Velocity) -> Pressure {

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_1 = TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);
    let rho1: MassDensity = state_1.get_rho();
    let mach_number_at_inlet: Ratio = state_1.get_mach_number(v1);
    let da = a2 - a1;
    let ratio_one = Ratio::new::<ratio>(1.0);

    let dp = da/a1 * (rho1 * v1 * v1)/(ratio_one - mach_number_at_inlet *mach_number_at_inlet);

    dp

}


// for isentropic/adiabatic diffusers/nozzles, we follow the mach number
// formulae 
//
pub fn get_dv_isentropic_nozzle_diffuser(
    a1: Area,
    a2: Area,
    p1: Pressure,
    h1: AvailableEnergy,
    v1: Velocity) -> Velocity {

    let ref_vol = Volume::new::<cubic_meter>(1.0);
    let state_1 = TampinesSteamTableCV::new_from_ph(p1, h1, ref_vol);
    let mach_number_at_inlet: Ratio = state_1.get_mach_number(v1);
    let da = a2 - a1;
    let ratio_one = Ratio::new::<ratio>(1.0);

    let dv = -da/a1 * v1/(ratio_one - mach_number_at_inlet *mach_number_at_inlet);

    dv

}
