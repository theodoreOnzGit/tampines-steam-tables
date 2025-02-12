use t_ps_flash::{t_ps_3a, t_ps_3b};
use uom::si::{f64::*, specific_heat_capacity::kilojoule_per_kilogram_kelvin};
use v_ps_flash::{v_ps_3a, v_ps_3b};
#[inline]
pub fn s_3a3b_backwards_ps_boundary() -> SpecificHeatCapacity {
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
        4.412_021_482_234_76
    );
    s

}

pub mod v_ps_flash;
pub mod t_ps_flash;

#[inline]
pub fn v_ps_3(p: Pressure, s: SpecificHeatCapacity) -> SpecificVolume {

    let s_boundary = s_3a3b_backwards_ps_boundary();

    if s > s_boundary {
        return v_ps_3b(p, s);
    } else {
        return v_ps_3a(p,s);
    };
}


#[inline]
pub fn t_ps_3(p: Pressure, s: SpecificHeatCapacity) -> ThermodynamicTemperature {

    let s_boundary = s_3a3b_backwards_ps_boundary();

    if s > s_boundary {
        return t_ps_3b(p, s);
    } else {
        return t_ps_3a(p,s);
    };
}
