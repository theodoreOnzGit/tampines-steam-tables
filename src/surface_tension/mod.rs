use uom::si::{f64::*, radiant_exposure::joule_per_square_meter, ratio::ratio};

use crate::constants::t_crit_water;

/// function for surface tension 
/// units are newtons per meter
///
/// newtons = kg * m s^(-2)
///
/// so newtons per meter is 
/// newtons/m = kg * s^(-2)
///
/// this is the same unit as RadiantExposure
/// Joule per m^2
///
/// However, the UOM crate doesn't have surface tension per se
/// So I'll use RadiantExposure as the return type
pub fn water_surf_tension(t: ThermodynamicTemperature) -> RadiantExposure {
    let sigma_star = RadiantExposure::new::<joule_per_square_meter>(1.0e-3);
    let t_crit = t_crit_water();

    let theta: f64 = (t/t_crit).get::<ratio>();

    let one_minus_theta = 1.0 - theta;

    let dimensionless_surf_tension = 
        235.8 * one_minus_theta.powf(1.256) * (1.0 - 0.625 * one_minus_theta);

    return dimensionless_surf_tension * sigma_star;

}

#[cfg(test)]
mod tests;
