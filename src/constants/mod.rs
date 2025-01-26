use uom::si::{f64::*, specific_heat_capacity::kilojoule_per_kilogram_kelvin};

/// gas constant for water 
pub const R_KJ_PER_KG_KELVIN: f64 = 0.461526;
/// critical temp for water
pub const T_C_KELVIN: f64 = 647.096;
/// critical pressure for water
pub const P_C_MPA: f64 = 22.064;
/// critical vol for water
pub const RHO_C_KG_PER_M3: f64 = 322.0;


/// triple pt temp for water
pub const T_TRIPLE_PT_KELVIN: f64 = 273.16;
/// triple pt pressure for water
pub const P_TRIPLE_PT_PASCAL: f64 = 611.657;

/// boiling pt temp for water at 1 atm (normal condition)
pub const T_NORMAL_BP_KELVIN: f64 = 373.1243;
/// 1 atmosphere, this is the pressure for normal boiling pt
pub const P_NORMAL_BP_MPA: f64 = 0.101325;

/// molecular weight of water in g/gmol
pub const MOLAR_MASS_WATER_G_PER_GMOL: f64 = 18.015257;
/// molar gas constant (R_M) in joules/(mol kelvin)
pub const R_M_J_PER_MOL_KELVIN: f64 = 8.31451;

// dimensioned properties 

/// returns the specific gas constant of water 
/// in proper dimensioned units using uom
#[inline]
pub fn specific_gas_constant_of_water() -> SpecificHeatCapacity {
    let r = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
        R_KJ_PER_KG_KELVIN
    );

    r
}
