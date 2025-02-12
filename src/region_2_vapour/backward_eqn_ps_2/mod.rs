use subregion_2a::t_ps_2a;
use subregion_2b::t_ps_2b;
use subregion_2c::t_ps_2c;
use uom::si::{f64::*, pressure::megapascal, specific_heat_capacity::kilojoule_per_kilogram_kelvin};

pub mod subregion_2a;
pub mod subregion_2b;
pub mod subregion_2c;

#[inline]
pub fn t_ps_2(p: Pressure, s: SpecificHeatCapacity) -> ThermodynamicTemperature {

    let p_boundary_2a2b = Pressure::new::<megapascal>(4.0);
    let s_boundary_2b2c = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(5.85);

    if p <= p_boundary_2a2b {
        return t_ps_2a(p, s);
    };
    
    if s >= s_boundary_2b2c {
        return t_ps_2b(p, s);
    } else {
        return t_ps_2c(p, s);
    };
}
