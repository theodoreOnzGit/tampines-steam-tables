use float_equations::t_ps_1_kelvin;
use uom::si::{f64::*, pressure::pascal, specific_heat_capacity::kilojoule_per_kilogram_kelvin, thermodynamic_temperature::kelvin};

/// float equations from rust-steam
/// legacy and ported over
pub mod float_equations;


fn theta_1(p: Pressure, s: SpecificHeatCapacity) -> f64 {
    let s_kj_per_kg_k = s.get::<kilojoule_per_kilogram_kelvin>();
    let p_pascals = p.get::<pascal>();
    let t_float_kelvin = t_ps_1_kelvin(p_pascals, s_kj_per_kg_k);

    return t_float_kelvin * 1.0;
}

pub fn t_ps_1(p: Pressure, s: SpecificHeatCapacity) -> ThermodynamicTemperature {
    let theta_1 = theta_1(p, s);

    return theta_1 * ThermodynamicTemperature::new::<kelvin>(1.0);
}

