use crate::constants::{p_crit_water, t_crit_water};

use super::*;

pub fn region_fwd_eqn_two_phase(
    t: ThermodynamicTemperature,
    p: Pressure,
    steam_quality: f64) -> FwdEqnRegion {

    // for the two phase version, I need to determine if it's in the 
    // saturation region 
    // otherwise, revert back to the normal boundary 
    //
    // if in the saturation region, then x 
    // should be between 0 and 1

    let mut x = steam_quality;

    if x < 0.0 {
        x = 0.0;
    };
    if x > 1.0 {
        x = 1.0;
    };

    // now that that's settled, let's check the 
    // points
    //
    // if beyond critical point, don't even need to consider two phase
    //
    // this will cover region 5 automatically
    if t > t_crit_water() {
        return region_fwd_eqn_single_phase(t, p);
    };

    if p > p_crit_water() {
        return region_fwd_eqn_single_phase(t, p);
    };

    // if below critical temperature and pressure
    // check if it lies within the sat temperature and pressure 

    let p_sat_reg4 = sat_pressure_4(t);

    // if pressure is at exactly the saturation pressure,
    // and steam quality is between 0 and 1
    // then return region 4
    
    let multiphase_steam: bool = x > 0.0 && x < 1.0;

    if p == p_sat_reg4 && multiphase_steam {
        return FwdEqnRegion::Region4;
    };

    // if this condition is not satisfied, return single phase
    
    return region_fwd_eqn_single_phase(t, p);

}
