use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::ratio::ratio;
use uom::si::f64::*;
use uom::si::available_energy::kilojoule_per_kilogram;

/// this is for eq 2.44 on page 84
/// based on table 2.73
const HB13_PRIME_S_BOUNDARY_EQN_COEFFS: [[f64; 3]; 6] = [
    [0.0, 0.0, 0.913_965_547_600_543],
    [1.0, -2.0, -0.430_944_856_041_991e-4],
    [1.0, 2.0, 0.603_235_694_765_419e2],
    [3.0, -12.0, 0.117_518_273_082_168e-17],
    [5.0, -4.0, 0.220_000_904_781_292],
    [6.0, -3.0, -0.690_815_545_851_641e2],
];


/// this function represents the saturated liquid line
/// for hs flashing between region 1 and region 4
pub fn hb13_s_boundary_enthalpy(
    s: SpecificHeatCapacity) -> AvailableEnergy {

    let s_ref = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(3.8);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(1700.0);
    let sigma: f64 = (s/s_ref).get::<ratio>();

    let mut eta: f64 = 0.0;

    for coeffs in HB13_PRIME_S_BOUNDARY_EQN_COEFFS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        eta += ni * (sigma - 0.884).powf(ii) * (sigma - 0.864).powf(ji);
    }

    return h_ref * eta;

}



