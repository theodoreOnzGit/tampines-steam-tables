
pub(crate) enum BackwardPTRegion3 {
    SubregionA,
    SubregionB,
    SubregionC,
    SubregionD,
    SubregionE,
    SubregionF,
    SubregionG,
    SubregionH,
    SubregionI,
    SubregionJ,
    SubregionK,
    SubregionL,
    SubregionM,
    SubregionN,
    SubregionO,
    SubregionP,
    SubregionQ,
    SubregionR,
    SubregionS,
    SubregionT,
    SubregionU,
    SubregionV,
    SubregionW,
    SubregionX,
    SubregionY,
    SubregionZ,
}

/// I imported these from rusteam,
/// as these seem to be tested already 
/// what I'm only adding is the use of 
/// units of measure (uom) to wrap these eqns
pub(crate) mod floating_point_eqns_for_specific_vol;
use floating_point_eqns_for_specific_vol::*;
use uom::si::{f64::*, pressure::pascal, specific_volume::cubic_meter_per_kilogram, thermodynamic_temperature::kelvin};


/// obtains volume for region 3 based on pt flash 
///
/// then using vt flash, you can get everything else
#[inline]
pub fn v_tp_3(t: ThermodynamicTemperature, p: Pressure) -> SpecificVolume {
    let t_kelvin = t.get::<kelvin>();
    let p_pascal = p.get::<pascal>();

    let v_m3_per_kg = v_tp_3_float(t_kelvin, p_pascal);

    return SpecificVolume::new::<cubic_meter_per_kilogram>(
        v_m3_per_kg);
}

/// these are needed for ph flashing at saturation line 
/// liq phase specific vol
/// from 623.15 K to 634.659 K
#[inline]
pub fn v_tp_3c(t: ThermodynamicTemperature, p: Pressure) -> SpecificVolume {
    let t_kelvin = t.get::<kelvin>();
    let p_pascal = p.get::<pascal>();

    let v_m3_per_kg = subregion_c(t_kelvin, p_pascal);

    return SpecificVolume::new::<cubic_meter_per_kilogram>(
        v_m3_per_kg);
}

///  these are needed for ph flashing at saturation line 
/// liq phase specific vol
/// from 634.659 K to 643.15 K
#[inline]
pub fn v_tp_3s(t: ThermodynamicTemperature, p: Pressure) -> SpecificVolume {
    let t_kelvin = t.get::<kelvin>();
    let p_pascal = p.get::<pascal>();

    let v_m3_per_kg = subregion_s(t_kelvin, p_pascal);

    return SpecificVolume::new::<cubic_meter_per_kilogram>(
        v_m3_per_kg);
}
/// these are needed for ph flashing at saturation line 
/// vap phase specific vol
/// from 623.15 K to 640.961 K
#[inline]
pub fn v_tp_3t(t: ThermodynamicTemperature, p: Pressure) -> SpecificVolume {
    let t_kelvin = t.get::<kelvin>();
    let p_pascal = p.get::<pascal>();

    let v_m3_per_kg = subregion_t(t_kelvin, p_pascal);

    return SpecificVolume::new::<cubic_meter_per_kilogram>(
        v_m3_per_kg);
}
/// these are needed for ph flashing at saturation line 
/// vap phase specific vol
/// from 640.931 K to 643.15
#[inline]
pub fn v_tp_3r(t: ThermodynamicTemperature, p: Pressure) -> SpecificVolume {
    let t_kelvin = t.get::<kelvin>();
    let p_pascal = p.get::<pascal>();

    let v_m3_per_kg = subregion_r(t_kelvin, p_pascal);

    return SpecificVolume::new::<cubic_meter_per_kilogram>(
        v_m3_per_kg);
}


pub mod intensive_properties;
pub use intensive_properties::*;
