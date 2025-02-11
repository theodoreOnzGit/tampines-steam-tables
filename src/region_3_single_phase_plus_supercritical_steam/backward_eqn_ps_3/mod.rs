use uom::si::{f64::*, specific_heat_capacity::kilojoule_per_kilogram_kelvin};
#[inline]
pub fn s_3a3b_backwards_ps_boundary() -> SpecificHeatCapacity {
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
        4.412_021_482_234_76
    );
    s

}

pub mod v_ps_flash;
pub mod t_ps_flash;
