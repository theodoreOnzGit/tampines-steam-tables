use std::ops::Index;

use uom::si::pressure::megapascal;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::thermal_conductivity::watt_per_meter_kelvin;
use uom::si::{f64::*, ratio::ratio};

use crate::constants::p_crit_water;
use crate::constants::rho_crit_water;
use crate::constants::t_crit_water;
use crate::dynamic_viscosity::psi_1_viscosity;
use crate::dynamic_viscosity::psi_0_viscosity;
use crate::interfaces::functional_programming::pt_flash_eqm::{cp_tp_eqm_single_phase, v_tp_eqm_two_phase};
use crate::interfaces::functional_programming::pt_flash_eqm::cv_tp_eqm_single_phase;
use crate::interfaces::functional_programming::pt_flash_eqm::kappa_t_tp_eqm;
use crate::interfaces::functional_programming::pt_flash_eqm::v_tp_eqm_single_phase;
use crate::region_5_steam_at_800_plus_degc::InversePressure;

const LAMBDA_0_COEFFS: [[f64; 2]; 5] = [
    [1.0,  0.244_322_1e-2],
    [2.0,  0.132_309_5e-1],
    [3.0,  0.677_035_7e-2],
    [4.0,  -0.345_458_6e-2],
    [5.0,  0.409_626_6e-3],
];

pub fn lambda_tp_eqm_single_phase(t: ThermodynamicTemperature,
    p: Pressure) -> ThermalConductivity {

    let rho = v_tp_eqm_single_phase(t, p).recip();
    let lambda_0 = lambda_0(t);
    let lambda_1 = lambda_1(rho, t);
    let lambda_2 = lambda_2_crit_enhancement_term_tp_single_phase(t, p);
    let lambda_star = ThermalConductivity::new::<watt_per_meter_kelvin>(1.0e-3);

    let dimensionless_lambda = lambda_0 * lambda_1 + lambda_2;

    return lambda_star * dimensionless_lambda;

}
pub fn lambda_tp_eqm_two_phase(t: ThermodynamicTemperature,
    p: Pressure,
    x: f64) -> ThermalConductivity {

    let rho = v_tp_eqm_two_phase(t, p, x).recip();
    let lambda_0 = lambda_0(t);
    let lambda_1 = lambda_1(rho, t);
    let lambda_2 = lambda_2_crit_enhancement_term_tp_two_phase_estimate(t, p, x);
    let lambda_star = ThermalConductivity::new::<watt_per_meter_kelvin>(1.0e-3);

    let dimensionless_lambda = lambda_0 * lambda_1 + lambda_2;

    return lambda_star * dimensionless_lambda;

}

pub(crate) fn lambda_0(t: ThermodynamicTemperature) -> f64 {
    let t_c = t_crit_water();
    let theta_f64: f64 = (t/t_c).get::<ratio>();

    let num = theta_f64.sqrt();

    let mut den = 0.0;

    for coeffs in LAMBDA_0_COEFFS {
        let i = coeffs[0];
        let ni = coeffs[1];

        den += ni * theta_f64.powf(1.0 - i);
    };

    return num/den;
}

const LAMBDA_1_COEFFS_NI1: [[f64; 2]; 5] = [
    [1.0,  0.160_397_357e1],
    [2.0,  0.233_771_842e1],
    [3.0,  0.219_650_529e1],
    [4.0,  -0.121_051_378e1],
    [5.0,  -0.272_033_700e1],
];

const LAMBDA_1_COEFFS_NI2: [[f64; 2]; 5] = [
    [1.0,  -0.646_013_523],
    [2.0,  -0.278_843_778e1],
    [3.0,  -0.454_580_785e1],
    [4.0,  0.160_812_989e1],
    [5.0,  0.457_586_331e1],
];

const LAMBDA_1_COEFFS_NI3: [[f64; 2]; 5] = [
    [1.0,  0.111_443_906],
    [2.0,  0.153_616_167e1],
    [3.0,  0.355_777_244e1],
    [4.0,  -0.621_178_141],
    [5.0,  -0.318_369_245e1],
];
const LAMBDA_1_COEFFS_NI4: [[f64; 2]; 5] = [
    [1.0,  0.102_997_357],
    [2.0,  -0.463_045_512],
    [3.0,  -0.140_944_978e1],
    [4.0,  0.716_373_224e-1],
    [5.0,  0.111_683_480e1],
];
const LAMBDA_1_COEFFS_NI5: [[f64; 2]; 5] = [
    [1.0,  -0.504_123_634e-1],
    [2.0,  0.832_827_019e-1],
    [3.0,  0.275_418_278],
    [4.0,  0.0],
    [5.0,  -0.192_683_050],
];
const LAMBDA_1_COEFFS_NI6: [[f64; 2]; 5] = [
    [1.0,  0.609_859_258e-2],
    [2.0,  -0.719_201_245e-2],
    [3.0,  -0.205_938_816e-1],
    [4.0,  0.0],
    [5.0,  0.129_138_420e-1],
];
pub(crate) fn lambda_1(rho: MassDensity,
    t: ThermodynamicTemperature,) -> f64 {
    let t_c = t_crit_water();
    let theta_f64: f64 = (t/t_c).get::<ratio>();
    let rho_c = rho_crit_water();
    let delta_f64: f64 = (rho/rho_c).get::<ratio>();

    fn inner_sum_over_all_j(i: usize, delta: f64) -> f64{

        let ni1 = LAMBDA_1_COEFFS_NI1.index(i-1)[1];
        let ni2 = LAMBDA_1_COEFFS_NI2.index(i-1)[1];
        let ni3 = LAMBDA_1_COEFFS_NI3.index(i-1)[1];
        let ni4 = LAMBDA_1_COEFFS_NI4.index(i-1)[1];
        let ni5 = LAMBDA_1_COEFFS_NI5.index(i-1)[1];
        let ni6 = LAMBDA_1_COEFFS_NI6.index(i-1)[1];

        let inner_sum: f64 = 
            ni1 * (delta - 1.0).powf(1.0 -  1.0)
            + ni2 * (delta - 1.0).powf(2.0 - 1.0)
            + ni3 * (delta - 1.0).powf(3.0 - 1.0)
            + ni4 * (delta - 1.0).powf(4.0 - 1.0)
            + ni5 * (delta - 1.0).powf(5.0 - 1.0)
            + ni6 * (delta - 1.0).powf(6.0 - 1.0);

        return inner_sum;

    }

    
    let mut exponent: f64 = 0.0;

    // I know I should be using a for loop, but this makes coding 
    // easier to visualise given the nested sum.. lol
    let i = 1;
    exponent += delta_f64 * (theta_f64.recip() - 1.0)
        .powi(i - 1) * 
        inner_sum_over_all_j(i as usize, delta_f64);
    let i = 2;
    exponent += delta_f64 * (theta_f64.recip() - 1.0)
        .powi(i - 1) * 
        inner_sum_over_all_j(i as usize, delta_f64);
    let i = 3;
    exponent += delta_f64 * (theta_f64.recip() - 1.0)
        .powi(i - 1) * 
        inner_sum_over_all_j(i as usize, delta_f64);
    let i = 4;
    exponent += delta_f64 * (theta_f64.recip() - 1.0)
        .powi(i - 1) * 
        inner_sum_over_all_j(i as usize, delta_f64);
    let i = 5;
    exponent += delta_f64 * (theta_f64.recip() - 1.0)
        .powi(i - 1) * 
        inner_sum_over_all_j(i as usize, delta_f64);


    exponent.exp()

}
pub(crate) fn lambda_2_crit_enhancement_term_tp_two_phase_estimate(
    t: ThermodynamicTemperature,
    p: Pressure,
    x: f64) -> f64 {

    let rho = v_tp_eqm_two_phase(t, p, x).recip();
    let t_c = t_crit_water();
    let theta_f64: f64 = (t/t_c).get::<ratio>();
    let rho_c = rho_crit_water();
    let delta_f64: f64 = (rho/rho_c).get::<ratio>();

    // this is dimensionless viscosity
    let psi = psi_0_viscosity(t) * psi_1_viscosity(t, rho);

    // these terms are independent of density
    let n1 = 0.177_851_4e3;
    let n2 = 0.636_619_772_367_581;
    let n3 = 0.135_882_142_589_674e1;
    let n4 = 0.508_474_576_271;
    let n5 = 1.5;

    // now, looks like we need to calculate certain properties 
    // such as cp, cv and kappa_t
    //
    // those require p,h or tp flashing in order to work outside region 3 
    //
    // we only have rho and t now
    // which doesn't exactly work outside region 3
    // so without iterations, this will be a problem.
    //
    // However, perhaps one can assume pressure is given since 
    // we often use tp or ph flashing, that makes things a lot easier
    //
    //
    // so a pt flash would be good.
    //
    // However, it won't work in region 4 as there are two phases 
    // to deal with
    //
    // it may be more reasonable to work with p,h flash from the get go
    //

    let mut cp = cp_tp_eqm_single_phase(t, p);

    if cp.get::<kilojoule_per_kilogram_kelvin>() < 0.0 {
        cp = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(1.0e13);
    } else if cp.get::<kilojoule_per_kilogram_kelvin>() > 1.0e13 {
        cp = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(1.0e13);
    };
    let cv = cv_tp_eqm_single_phase(t, p);
    let kappa_t = kappa_t_tp_eqm(t, p);

    let b: f64 = (cp/cv).get::<ratio>();
    let captial_a: f64 = captial_a(n2, n3, delta_f64, theta_f64, 
        kappa_t, n4, n5, b);


    let gas_constant_r = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
            0.461_518_05
        );
    
    let lambda_2 = n1 * delta_f64 * theta_f64 / psi * 
        cp/gas_constant_r * captial_a;


    return lambda_2.get::<ratio>();
}

pub(crate) fn lambda_2_crit_enhancement_term_tp_single_phase(
    t: ThermodynamicTemperature,
    p: Pressure) -> f64 {

    let rho = v_tp_eqm_single_phase(t, p).recip();
    let t_c = t_crit_water();
    let theta_f64: f64 = (t/t_c).get::<ratio>();
    let rho_c = rho_crit_water();
    let delta_f64: f64 = (rho/rho_c).get::<ratio>();

    // this is dimensionless viscosity
    let psi = psi_0_viscosity(t) * psi_1_viscosity(t, rho);

    // these terms are independent of density
    let n1 = 0.177_851_4e3;
    let n2 = 0.636_619_772_367_581;
    let n3 = 0.135_882_142_589_674e1;
    let n4 = 0.508_474_576_271;
    let n5 = 1.5;

    // now, looks like we need to calculate certain properties 
    // such as cp, cv and kappa_t
    //
    // those require p,h or tp flashing in order to work outside region 3 
    //
    // we only have rho and t now
    // which doesn't exactly work outside region 3
    // so without iterations, this will be a problem.
    //
    // However, perhaps one can assume pressure is given since 
    // we often use tp or ph flashing, that makes things a lot easier
    //
    //
    // so a pt flash would be good.
    //
    // However, it won't work in region 4 as there are two phases 
    // to deal with
    //
    // it may be more reasonable to work with p,h flash from the get go
    //

    let mut cp = cp_tp_eqm_single_phase(t, p);

    if cp.get::<kilojoule_per_kilogram_kelvin>() < 0.0 {
        cp = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(1.0e13);
    } else if cp.get::<kilojoule_per_kilogram_kelvin>() > 1.0e13 {
        cp = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(1.0e13);
    };
    let cv = cv_tp_eqm_single_phase(t, p);
    let kappa_t = kappa_t_tp_eqm(t, p);

    let b: f64 = (cp/cv).get::<ratio>();
    let captial_a: f64 = captial_a(n2, n3, delta_f64, theta_f64, 
        kappa_t, n4, n5, b);


    let gas_constant_r = 
        SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
            0.461_518_05
        );
    
    let lambda_2 = n1 * delta_f64 * theta_f64 / psi * 
        cp/gas_constant_r * captial_a;


    return lambda_2.get::<ratio>();
}

fn captial_b(delta: f64, theta: f64, kappa_t: InversePressure,
    n5: f64) -> f64 {
    let captial_c = captial_c(delta);
    let p_c = p_crit_water();
    let corrected_kappa_t: InversePressure;
    if kappa_t.value < 0.0 {
        corrected_kappa_t = 1.0e13 * Pressure::new::<megapascal>(1.0).recip();
    } else if kappa_t.recip().get::<megapascal>() < 1.0e13_f64.recip() {
        corrected_kappa_t = 1.0e13 * Pressure::new::<megapascal>(1.0).recip();
    } else {
        corrected_kappa_t = kappa_t;
    };


    let mut captial_b = (p_c * delta * corrected_kappa_t).get::<ratio>() 
        - n5 * theta.recip() * captial_c;

    if captial_b < 0.0 {
        captial_b = 0.0;
    }

    return captial_b;

    
}

fn captial_c(delta: f64) -> f64{

    let n1: f64;
    let n2: f64;
    let n3: f64;
    let n4: f64;
    let n5: f64;
    let n6: f64;

    // this takes from table 3.7
    //
    if delta <= 0.310_559_006 {
        // tested
        n1 = 0.653_786_807_199_516e1;
        n2 = -0.561_149_954_923_348e1;
        n3 = 0.339_624_167_361_325e1;
        n4 = -0.227_492_629_730_878e1;
        n5 = 0.102_631_854_662_709e2;
        n6 = 0.197_815_050_331_519e1;
    } else if delta <= 0.776_397_516 {
        n1 = 0.652_717_759_281_799e1;
        n2 = -0.630_816_983_387_575e1;
        n3 = 0.808_379_285_492_595e1;
        n4 = -0.982_240_510_197_603e1;
        n5 = 0.121_358_413_791_395e2;
        n6 = -0.554_349_664_571_295e1;
    } else if delta <= 1.242_236_025 {
        n1 = 0.535_500_529_896_124e1;
        n2 = -0.396_415_689_925_446e1;
        n3 = 0.891_990_208_918_795e1;
        n4 = -0.120_338_729_505_790e2;
        n5 = 0.919_494_865_194_302e1;
        n6 = -0.216_866_274_479_712e1;
    } else if delta <= 1.863_354_037 {
        // tested
        n1 = 0.155_225_959_906_681e1;
        n2 = 0.464_621_290_821_181;
        n3 = 0.893_237_374_861_479e1;
        n4 = -0.110_321_960_061_126e2;
        n5 = 0.616_780_999_933_360e1;
        n6 = -0.965_458_722_086_812;
    } else {
        // tested
        n1 = 0.111_999_926_419_994e1;
        n2 = 0.595_748_562_571_649;
        n3 = 0.988_952_565_078_920e1;
        n4 = -0.103_255_051_147_040e2;
        n5 = 0.466_861_294_457_414e1;
        n6 = -0.503_243_546_373_828;
    }

    let mut den = 0.0;

    den += n1 * delta.powf(1.0 - 1.0);
    den += n2 * delta.powf(2.0 - 1.0);
    den += n3 * delta.powf(3.0 - 1.0);
    den += n4 * delta.powf(4.0 - 1.0);
    den += n5 * delta.powf(5.0 - 1.0);
    den += n6 * delta.powf(6.0 - 1.0);

    return den.recip();
        

}

fn small_a(n3: f64, delta: f64, theta: f64, kappa_t: InversePressure, 
    n4: f64,
    n5: f64) -> f64 {
    let captial_b = captial_b(delta, theta, kappa_t, n5);

    let exponent = (delta * captial_b).powf(n4);

    return n3 * exponent;
}

fn captial_a(n2: f64, n3: f64, delta: f64,
    theta: f64, kappa_t: InversePressure, n4: f64,
    n5: f64, b: f64) -> f64 {
    let a = small_a(n3, delta, theta, kappa_t, n4, n5);

    // first guard clause
    if a < 1.2e-7 {
        return 0.0;
    };

    let n2_by_a = n2/a;

    let mut bracket_term: f64 = 0.0;
    bracket_term += (1.0 - b.recip()) * a.atan();
    bracket_term += a/b - 1.0;
    bracket_term += (
        -(a.recip() + 3.0_f64.recip() * a.powi(2) * delta.powi(-2))
    ).recip().exp();

    return n2_by_a * bracket_term;





}


#[cfg(test)]
mod tests;
