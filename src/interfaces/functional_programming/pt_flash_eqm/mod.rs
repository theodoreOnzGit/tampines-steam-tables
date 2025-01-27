use uom::si::{f64::*, pressure::pascal, thermodynamic_temperature::kelvin};

use crate::{region_1_subcooled_liquid::{alpha_v_tp_1, cp_tp_1, cv_tp_1, h_tp_1, kappa_t_tp_1, kappa_tp_1, s_tp_1, u_tp_1, v_tp_1, w_tp_1, InversePressure}, region_2_vapour::{alpha_v_tp_2, cp_tp_2, cv_tp_2, h_tp_2, kappa_t_tp_2, kappa_tp_2, s_tp_2, u_tp_2, v_tp_2, w_tp_2}, region_3_single_phase_plus_supercritical_steam::{alpha_v_tp_3, cp_tp_3, cv_tp_3, h_tp_3, kappa_t_tp_3, kappa_tp_3, p_boundary_2_3, s_tp_3, u_tp_3, v_tp_3, w_tp_3}, region_4_vap_liq_equilibrium::sat_pressure_4, region_5_superheated_steam::{alpha_v_tp_5, cp_tp_5, cv_tp_5, h_tp_5, kappa_t_tp_5, kappa_tp_5, s_tp_5, u_tp_5, v_tp_5, w_tp_5}};

#[derive(Debug,PartialEq, Eq, PartialOrd, Ord)]
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
        FwdEqnRegion::Region3 => h_tp_3(t, p),
        FwdEqnRegion::Region4 => todo!("cannot find enthalpy of mixture without steam quality"),
        FwdEqnRegion::Region5 => h_tp_5(t, p),
    }
}

/// returns the internal energy given temperature and pressure
pub fn u_tp_eqm(t: ThermodynamicTemperature, p: Pressure) -> AvailableEnergy {
    let region = region_fwd_eqn(t, p);

    match region {
        FwdEqnRegion::Region1 => u_tp_1(t, p),
        FwdEqnRegion::Region2 => u_tp_2(t, p),
        FwdEqnRegion::Region3 => u_tp_3(t, p),
        FwdEqnRegion::Region4 => todo!("cannot find enthalpy of mixture without steam quality"),
        FwdEqnRegion::Region5 => u_tp_5(t, p),
    }
}


/// returns the specific entropy given temperature and pressure
pub fn s_tp_eqm(t: ThermodynamicTemperature, p: Pressure) -> SpecificHeatCapacity {
    let region = region_fwd_eqn(t, p);

    match region {
        FwdEqnRegion::Region1 => s_tp_1(t, p),
        FwdEqnRegion::Region2 => s_tp_2(t, p),
        FwdEqnRegion::Region3 => s_tp_3(t, p),
        FwdEqnRegion::Region4 => todo!("cannot find entropy of mixture without steam quality"),
        FwdEqnRegion::Region5 => s_tp_5(t, p),
    }
}

/// returns the isobaric (const pressure) heat capacitygiven temperature and pressure
pub fn cp_tp_eqm(t: ThermodynamicTemperature, p: Pressure) -> SpecificHeatCapacity {
    let region = region_fwd_eqn(t, p);

    match region {
        FwdEqnRegion::Region1 => cp_tp_1(t, p),
        FwdEqnRegion::Region2 => cp_tp_2(t, p),
        FwdEqnRegion::Region3 => cp_tp_3(t, p),
        FwdEqnRegion::Region4 => todo!("cannot find cp of mixture without steam quality"),
        FwdEqnRegion::Region5 => cp_tp_5(t, p),
    }
}


/// returns the isochoric (const vol) heat capacity given temperature and pressure
pub fn cv_tp_eqm(t: ThermodynamicTemperature, p: Pressure) -> SpecificHeatCapacity {
    let region = region_fwd_eqn(t, p);

    match region {
        FwdEqnRegion::Region1 => cv_tp_1(t, p),
        FwdEqnRegion::Region2 => cv_tp_2(t, p),
        FwdEqnRegion::Region3 => cv_tp_3(t, p),
        FwdEqnRegion::Region4 => todo!("cannot find cv of mixture without steam quality"),
        FwdEqnRegion::Region5 => cv_tp_5(t, p),
    }
}




/// returns the specific volume given temperature and pressure
pub fn v_tp_eqm(t: ThermodynamicTemperature, p: Pressure) -> SpecificVolume {
    let region = region_fwd_eqn(t, p);

    match region {
        FwdEqnRegion::Region1 => v_tp_1(t, p),
        FwdEqnRegion::Region2 => v_tp_2(t, p),
        // note that for region 3, the backward eqn is used
        FwdEqnRegion::Region3 => v_tp_3(t, p),
        FwdEqnRegion::Region4 => todo!("cannot find specific vol of mixture without steam quality"),
        FwdEqnRegion::Region5 => v_tp_5(t, p),
    }
}


/// returns the speed of sound given temperature and pressure
pub fn w_tp_eqm(t: ThermodynamicTemperature, p: Pressure) -> Velocity {
    let region = region_fwd_eqn(t, p);

    match region {
        FwdEqnRegion::Region1 => w_tp_1(t, p),
        FwdEqnRegion::Region2 => w_tp_2(t, p),
        FwdEqnRegion::Region3 => w_tp_3(t, p),
        FwdEqnRegion::Region4 => todo!("cannot find speed of sound of mixture without steam quality"),
        FwdEqnRegion::Region5 => w_tp_5(t, p),
    }
}


/// returns the isentropic exponent 
pub fn kappa_tp_eqm(t: ThermodynamicTemperature, p: Pressure) -> Ratio {
    let region = region_fwd_eqn(t, p);

    match region {
        FwdEqnRegion::Region1 => kappa_tp_1(t, p),
        FwdEqnRegion::Region2 => kappa_tp_2(t, p),
        FwdEqnRegion::Region3 => kappa_tp_3(t, p),
        FwdEqnRegion::Region4 => todo!("cannot find isentropic exponent of mixture without steam quality"),
        FwdEqnRegion::Region5 => kappa_tp_5(t, p),
    }
}

/// returns the isobaric cubic expansion coefficient
pub fn alpha_v_tp_eqm(t: ThermodynamicTemperature, p: Pressure) -> TemperatureCoefficient {
    let region = region_fwd_eqn(t, p);

    match region {
        FwdEqnRegion::Region1 => alpha_v_tp_1(t, p),
        FwdEqnRegion::Region2 => alpha_v_tp_2(t, p),
        FwdEqnRegion::Region3 => alpha_v_tp_3(t, p),
        FwdEqnRegion::Region4 => todo!("cannot find isobaric cubic exp coeff of mixture without steam quality"),
        FwdEqnRegion::Region5 => alpha_v_tp_5(t, p),
    }
}


/// returns the isothermal compressibility
pub fn kappa_t_tp_eqm(t: ThermodynamicTemperature, p: Pressure) -> InversePressure {
    let region = region_fwd_eqn(t, p);

    match region {
        FwdEqnRegion::Region1 => kappa_t_tp_1(t, p),
        FwdEqnRegion::Region2 => kappa_t_tp_2(t, p),
        FwdEqnRegion::Region3 => kappa_t_tp_3(t, p),
        FwdEqnRegion::Region4 => todo!("cannot find isothermal compressibility of mixture without steam quality"),
        FwdEqnRegion::Region5 => kappa_t_tp_5(t, p),
    }
}

/// re-exports the relative pressure coeff function for 
/// region 3 relative pressure coeff (other regions don't have it)
pub use crate::region_3_single_phase_plus_supercritical_steam::alpha_p_rho_t_3;

/// re-exports the isothermal stress coeff function for 
/// region 3 isothermal stress coeff (other regions don't have it)
pub use crate::region_3_single_phase_plus_supercritical_steam::beta_p_rho_t_3;

