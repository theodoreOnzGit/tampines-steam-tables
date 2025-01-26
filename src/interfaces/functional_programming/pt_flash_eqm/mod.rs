use uom::si::{f64::*, pressure::pascal, thermodynamic_temperature::kelvin};

use crate::{region_1_subcooled_liquid::h_tp_1, region_2_vapour::h_tp_2, region_3_single_phase_plus_supercritical_steam::p_boundary_2_3, region_4_vap_liq_equilibrium::sat_pressure_4, region_5_superheated_steam::h_tp_5};

/// an enum to help represent the appropriate 
/// regions in the forward equations
pub enum FwdEqnRegion {
    /// this is from T = 273.15 K to T=623.15K 
    Region1,
    Region2,
    Region3,
    Region4,
    Region5,
}

/// Determines which region of the pT chart
/// a point belongs to.
///
/// Returns an error if the point is outside the
/// bounds of the IAPWS-IF97 correlations.
///
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
///
pub fn region_fwd_eqn(t: ThermodynamicTemperature, p: Pressure) -> FwdEqnRegion {
    let p_sat_reg4 = sat_pressure_4(t);

    let p_boundary_23 = p_boundary_2_3(t);

    let t_kelvin = t.get::<kelvin>();
    let p_pascal = p.get::<pascal>();
    let p_boundary_23_pascal = p_boundary_23.get::<pascal>();
    let p_sat_reg4_pascal = p_sat_reg4.get::<pascal>();

    match (t_kelvin, p_pascal) {
        (temp, pres)
            if (1073.15..=2273.15).contains(&temp) && (0.0..=50.0e6).contains(&pres) =>
            {
                FwdEqnRegion::Region5
            }
        (temp, pres) if (273.15..647.096).contains(&temp) && pres == p_sat_reg4_pascal => {
            FwdEqnRegion::Region4
        }
        (temp, pres)
            if (623.15..=863.15).contains(&temp) && (p_boundary_23_pascal..100e6).contains(&pres) =>
            {
                FwdEqnRegion::Region3
            }
        (temp, pres)
            if ((273.15..=623.15).contains(&temp) && (0.0..=p_sat_reg4_pascal).contains(&pres))
                || ((623.15..=863.15).contains(&temp)
                    && (0.0..=p_boundary_23_pascal).contains(&pres))
                    || ((863.15..=1073.15).contains(&temp) && (0.0..100e6).contains(&pres)) =>
            {
                FwdEqnRegion::Region2
            }
        (temp, pres)
            if (273.15..=623.15).contains(&temp) && (p_sat_reg4_pascal..=100e6).contains(&pres) =>
            {
                FwdEqnRegion::Region1
            }
        _ => panic!("out of bounds!"),
    }
}

/// returns the enthalpy given temperature and pressure
pub fn h_tp_eqm(t: ThermodynamicTemperature, p: Pressure) -> AvailableEnergy {
    let region = region_fwd_eqn(t, p);

    match region {
        FwdEqnRegion::Region1 => h_tp_1(t, p),
        FwdEqnRegion::Region2 => h_tp_2(t, p),
        FwdEqnRegion::Region3 => todo!("(t,p) flash for region 3 not yet implemented"),
        FwdEqnRegion::Region4 => todo!("cannot find enthalpy of mixture without steam quality"),
        FwdEqnRegion::Region5 => h_tp_5(t, p),
    }
}
