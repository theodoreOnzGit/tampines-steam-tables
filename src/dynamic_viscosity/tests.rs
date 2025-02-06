use uom::si::{dynamic_viscosity::pascal_second, f64::*, mass_density::kilogram_per_cubic_meter, pressure::megapascal, thermodynamic_temperature::kelvin};

use crate::dynamic_viscosity::{psi_1_viscosity, water_viscosity_p_t, water_viscosity_rho_t};

use super::psi_0_viscosity;
#[test]
pub fn psi_0_viscosity_test(){
    let t1 = ThermodynamicTemperature::new::<kelvin>(298.15);
    let t2 = ThermodynamicTemperature::new::<kelvin>(873.15);
    let t3 = ThermodynamicTemperature::new::<kelvin>(673.15);
    let t4 = ThermodynamicTemperature::new::<kelvin>(1173.15);

    let psi_0_1 = 0.970_904_552e1;
    let psi_0_2 = 0.326_046_811e2;
    let psi_0_3 = 0.244_558_002e2;
    let psi_0_4 = 0.441_936_611e2;

    let psi_test = psi_0_viscosity(t1);
    approx::assert_relative_eq!(
        psi_0_1,
        psi_test,
        max_relative=1e-7
        );
    let psi_test = psi_0_viscosity(t2);
    approx::assert_relative_eq!(
        psi_0_2,
        psi_test,
        max_relative=1e-8
        );
    let psi_test = psi_0_viscosity(t3);
    approx::assert_relative_eq!(
        psi_0_3,
        psi_test,
        max_relative=1e-8
        );
    let psi_test = psi_0_viscosity(t4);
    approx::assert_relative_eq!(
        psi_0_4,
        psi_test,
        max_relative=1e-8
        );

}


#[test]
pub fn psi_1_viscosity_test(){
    let t1 = ThermodynamicTemperature::new::<kelvin>(298.15);
    let t2 = ThermodynamicTemperature::new::<kelvin>(873.15);
    let t3 = ThermodynamicTemperature::new::<kelvin>(673.15);
    let t4 = ThermodynamicTemperature::new::<kelvin>(1173.15);

    let rho1 = MassDensity::new::<kilogram_per_cubic_meter>(0.997_047_435e3);
    let rho2 = MassDensity::new::<kilogram_per_cubic_meter>(0.549_921_814e2);
    let rho3 = MassDensity::new::<kilogram_per_cubic_meter>(0.612_391_201e3);
    let rho4 = MassDensity::new::<kilogram_per_cubic_meter>(0.377_584_847e2);

    let psi_1_1 = 0.916_694_207e2;
    let psi_1_2 = 0.104_200_938e1;
    let psi_1_3 = 0.296_900_349e1;
    let psi_1_4 = 0.102_422_324e1;

    let psi_test = psi_1_viscosity(t1, rho1);
    approx::assert_relative_eq!(
        psi_1_1,
        psi_test,
        max_relative=1e-8
        );
    let psi_test = psi_1_viscosity(t2,rho2);
    approx::assert_relative_eq!(
        psi_1_2,
        psi_test,
        max_relative=1e-8
        );
    let psi_test = psi_1_viscosity(t3,rho3);
    approx::assert_relative_eq!(
        psi_1_3,
        psi_test,
        max_relative=1e-8
        );
    let psi_test = psi_1_viscosity(t4,rho4);
    approx::assert_relative_eq!(
        psi_1_4,
        psi_test,
        max_relative=1e-8
        );

}

#[test]
pub fn psi_1_viscosity_rho_t_test(){
    let t1 = ThermodynamicTemperature::new::<kelvin>(298.15);
    let t2 = ThermodynamicTemperature::new::<kelvin>(873.15);
    let t3 = ThermodynamicTemperature::new::<kelvin>(673.15);
    let t4 = ThermodynamicTemperature::new::<kelvin>(1173.15);

    let rho1 = MassDensity::new::<kilogram_per_cubic_meter>(0.997_047_435e3);
    let rho2 = MassDensity::new::<kilogram_per_cubic_meter>(0.549_921_814e2);
    let rho3 = MassDensity::new::<kilogram_per_cubic_meter>(0.612_391_201e3);
    let rho4 = MassDensity::new::<kilogram_per_cubic_meter>(0.377_584_847e2);

    let eta_pa_s_1 = 0.890_022_551e-3;
    let eta_pa_s_2 = 0.339_743_835e-4;
    let eta_pa_s_3 = 0.726_093_560e-4;
    let eta_pa_s_4 = 0.452_641_750e-4;

    let eta_pa_s = water_viscosity_rho_t(t1, rho1)
        .get::<pascal_second>();
    approx::assert_relative_eq!(
        eta_pa_s_1,
        eta_pa_s,
        max_relative=1e-8
        );
    let eta_pa_s = water_viscosity_rho_t(t2,rho2)
        .get::<pascal_second>();
    approx::assert_relative_eq!(
        eta_pa_s_2,
        eta_pa_s,
        max_relative=1e-8
        );
    let eta_pa_s = water_viscosity_rho_t(t3,rho3)
        .get::<pascal_second>();
    approx::assert_relative_eq!(
        eta_pa_s_3,
        eta_pa_s,
        max_relative=1e-8
        );
    let eta_pa_s = water_viscosity_rho_t(t4,rho4)
        .get::<pascal_second>();
    approx::assert_relative_eq!(
        eta_pa_s_4,
        eta_pa_s,
        max_relative=1e-8
        );

}


#[test]
pub fn psi_1_viscosity_p_t_test(){
    let t1 = ThermodynamicTemperature::new::<kelvin>(298.15);
    let t2 = ThermodynamicTemperature::new::<kelvin>(873.15);
    let t3 = ThermodynamicTemperature::new::<kelvin>(673.15);
    let t4 = ThermodynamicTemperature::new::<kelvin>(1173.15);

    let _rho1 = MassDensity::new::<kilogram_per_cubic_meter>(0.997_047_435e3);
    let _rho2 = MassDensity::new::<kilogram_per_cubic_meter>(0.549_921_814e2);
    let _rho3 = MassDensity::new::<kilogram_per_cubic_meter>(0.612_391_201e3);
    let _rho4 = MassDensity::new::<kilogram_per_cubic_meter>(0.377_584_847e2);


    let p1 = Pressure::new::<megapascal>(0.1);
    let p2 = Pressure::new::<megapascal>(20.0);
    let p3 = Pressure::new::<megapascal>(60.0);
    let p4 = Pressure::new::<megapascal>(20.0);

    let eta_pa_s_1 = 0.890_022_551e-3;
    let eta_pa_s_2 = 0.339_743_835e-4;
    let eta_pa_s_3 = 0.726_093_560e-4;
    let eta_pa_s_4 = 0.452_641_750e-4;

    let eta_pa_s = water_viscosity_p_t(t1, p1)
        .get::<pascal_second>();
    approx::assert_relative_eq!(
        eta_pa_s_1,
        eta_pa_s,
        max_relative=1e-8
        );
    let eta_pa_s = water_viscosity_p_t(t2, p2)
        .get::<pascal_second>();
    approx::assert_relative_eq!(
        eta_pa_s_2,
        eta_pa_s,
        max_relative=1e-8
        );
    let eta_pa_s = water_viscosity_p_t(t3, p3)
        .get::<pascal_second>();
    approx::assert_relative_eq!(
        eta_pa_s_3,
        eta_pa_s,
        max_relative=1e-8
        );
    let eta_pa_s = water_viscosity_p_t(t4, p4)
        .get::<pascal_second>();
    approx::assert_relative_eq!(
        eta_pa_s_4,
        eta_pa_s,
        max_relative=1e-8
        );

}
