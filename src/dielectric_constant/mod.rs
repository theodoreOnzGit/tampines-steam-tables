use uom::si::{f64::*, molar_mass::kilogram_per_mole, ratio::ratio, thermodynamic_temperature::kelvin};

use crate::constants::{avogadro_number_na, boltzmann_constant_k, molecular_dipole_moment_mu, permittivity_of_vacuum_eps_0, rho_crit_water, t_crit_water, water_mean_molecular_polarisability_alpha};

pub fn water_dielectric_const_rho_t(rho: MassDensity, t: ThermodynamicTemperature) -> f64 {
    let capital_a = captial_a(rho, t);
    let capital_b = captial_b(rho);

    let exponent = 9.0 
        + 2.0 * capital_a 
        + 18.0 * capital_b
        + capital_a.powi(2) 
        + 10.0 * capital_a * capital_b 
        + 9.0 * capital_b.powi(2);
    let den = 4.0 * (1.0 - capital_b);

    let num = 1.0 
        + capital_a 
        + 5.0 * capital_b 
        + exponent.sqrt();
    
    return num/den;
}

/// molar mass only for dielectric constant calculations,
/// slightly different from molar mass in constants page
#[inline]
fn molar_mass_water_for_dielectric_const() -> MolarMass {
    MolarMass::new::<kilogram_per_mole>(0.018_015_268)
}

/// based on table 3.11
const G_BAR_COEFFS_DIELECTRIC_CONST: [[f64; 3]; 11] = [
    [1.0, 0.25, 0.978_224_486_826],
    [1.0, 1.0, -0.957_771_379_375],
    [1.0, 2.5, 0.237_511_794_148],
    [2.0, 1.5, 0.714_692_244_396],
    [3.0, 1.5, -0.298_217_036_956],
    [3.0, 2.5, -0.108_863_472_196],
    [4.0, 2.0, 0.949_327_488_264e-1],
    [5.0, 2.0, -0.980_469_816_509e-2],
    [6.0, 5.0, 0.165_167_634_970e-4],
    [7.0, 0.5, 0.937_359_795_772e-4],
    [10.0, 10.0, -0.123_179_218_720e-9],
];

#[inline]
fn captial_a(rho: MassDensity, t: ThermodynamicTemperature) -> f64 {
    // constants
    let na = avogadro_number_na();
    let mu = molecular_dipole_moment_mu();
    let mol_wt_water = molar_mass_water_for_dielectric_const();
    let epsilon_0 = permittivity_of_vacuum_eps_0();
    let k = boltzmann_constant_k();
    let rho_c = rho_crit_water();
    let t_c = t_crit_water();
    let n12: f64 = 0.196_096_504_426e-2;

    // intermediate quantities
    let tau: f64 = (t_c/t).get::<ratio>();
    let delta: f64 = (rho/rho_c).get::<ratio>();
    let t_c_by_228_k: f64 = t_c.get::<kelvin>()/(228.0);

    // harris alder g-bar factor
    let mut g_bar: f64 = 1.0;
    // seems the term 2 was OUTSIDE the summation
    //
    // thank you Jesus and EJT
    let term_2 = n12 * delta * (t_c_by_228_k * tau.recip() - 1.0).powf(-1.2);
    g_bar += term_2;
    for coeffs in G_BAR_COEFFS_DIELECTRIC_CONST {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        // wow, even just term 1 by itself gets u pretty close!
        let term_1 = ni * delta.powf(ii) * tau.powf(ji);


        g_bar += term_1; 


    }



    let num = na * mu * mu * rho * g_bar;
    let den = mol_wt_water * epsilon_0 * k * t;

    return (num/den).get::<ratio>();

}


#[inline]
fn captial_b(rho: MassDensity) -> f64 {

    let na = avogadro_number_na();
    let alpha = water_mean_molecular_polarisability_alpha();
    let mol_wt_water = molar_mass_water_for_dielectric_const();
    let epsilon_0 = permittivity_of_vacuum_eps_0();

    let num = na * alpha * rho;
    let den = 3.0 * mol_wt_water * epsilon_0;

    return (num/den).get::<ratio>();
}

#[cfg(test)]
mod tests;
