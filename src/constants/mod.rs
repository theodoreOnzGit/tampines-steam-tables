use uom::si::amount_of_substance::mole;
use uom::si::area::square_meter;
use uom::si::electric_charge::coulomb;
use uom::si::electric_dipole_moment::coulomb_meter;
use uom::si::electric_permittivity::farad_per_meter;
use uom::si::energy::joule;
use uom::si::ratio::ratio;
use uom::si::{f64::*, heat_capacity::joule_per_kelvin, mass_density::kilogram_per_cubic_meter, pressure::megapascal, specific_heat_capacity::kilojoule_per_kilogram_kelvin, thermodynamic_temperature::kelvin};


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

/// returns the dimensioned critical temperature of water
#[inline]
pub fn t_crit_water() -> ThermodynamicTemperature {
    ThermodynamicTemperature::new::<kelvin>(T_C_KELVIN)
}

/// returns the dimensioned critical pressure of water 
#[inline]
pub fn p_crit_water() -> Pressure {
    Pressure::new::<megapascal>(P_C_MPA)
}


/// returns the dimensioned critical density of water
#[inline]
pub fn rho_crit_water() -> MassDensity {
    MassDensity::new::<kilogram_per_cubic_meter>(RHO_C_KG_PER_M3)
}


/// returns dimensioned boltzmann_constant_k
#[inline]
pub fn boltzmann_constant_k() -> HeatCapacity {

    HeatCapacity::new::<joule_per_kelvin>(1.380_658e-23)
}

// to make the inverse pressure type 
// it is m s^2 / kg 
use uom::si::{ISQ, SI, Quantity};
use uom::typenum::{Z0, N1, N3, P4, P2};

// quantity is defined
// ## Generic Parameters
// * `L`: Length dimension.
// * `M`: Mass dimension.
// * `T`: Time dimension.
// * `I`: Electric current dimension.
// * `Th`: Thermodynamic temperature dimension.
// * `N`: Amount of substance dimension.
// * `J`: Luminous intensity dimension.
// * `K`: Kind.
/// type alias for per mole for avogadro_number_na 
/// mol^(-1)
pub(crate) type PerMole = Quantity<ISQ<Z0, Z0, Z0, Z0, Z0, N1, Z0>, SI<f64>, f64>;
/// returns dimensioned avogadro number 
#[inline]
pub fn avogadro_number_na() -> PerMole {
    let one_mole = AmountOfSubstance::new::<mole>(1.0);
    return Ratio::new::<ratio>(6.022_136_7e23)/one_mole;
}

/// returns dimensioned electric permittivity of vacuum
#[inline]
pub fn permittivity_of_vacuum_eps_0() -> ElectricPermittivity {
    ElectricPermittivity::new::<farad_per_meter>(8.854_187_817e-12)
}


/// returns molectular dipole moment
#[inline]
pub fn molecular_dipole_moment_mu() -> ElectricDipoleMoment {
    ElectricDipoleMoment::new::<coulomb_meter>(6.138e-30)
}
/// units for MolecularPolarisability,
/// here's the working:
/// this is C^2 J^(-1) m^(2)
///
/// joule is 
/// kg m^2 s^(-2)
/// per joule is 
/// kg^(-1) m^(-2) s^2
///
/// together is C^2 kg^(-1) m^(-2) s^2 m^(2)
/// C^2 kg^(-1) s^2 
///
/// coulomb is Ampere second
/// A^2 kg^(-1) s^4
/// ## Generic Parameters
/// * `L`: Length dimension.
/// * `M`: Mass dimension.
/// * `T`: Time dimension.
/// * `I`: Electric current dimension.
/// * `Th`: Thermodynamic temperature dimension.
/// * `N`: Amount of substance dimension.
/// * `J`: Luminous intensity dimension.
/// * `K`: Kind.
/// type alias for per mole for molecular polarisability
pub(crate) type MolecularPolarisability = Quantity<ISQ<Z0, N1, P4, P2, Z0, Z0, Z0>, SI<f64>, f64>;

/// returns dimensioned mean molecular polarisability alpha 
#[inline]
pub fn water_mean_molecular_polarisability_alpha() -> MolecularPolarisability {
    let one_coulomb = ElectricCharge::new::<coulomb>(1.0);
    let one_joule = Energy::new::<joule>(1.0);
    let one_m2 = Area::new::<square_meter>(1.0);

    let dimensioned_correct_unit = one_coulomb * one_coulomb * one_m2 / one_joule;

    return 1.636e-40_f64 * dimensioned_correct_unit;
}
