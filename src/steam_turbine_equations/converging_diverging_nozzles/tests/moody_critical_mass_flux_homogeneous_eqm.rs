// for the Marviken tests 
//
// there was a vessel at 425 m3  (fig 3.1 on page 9)
// https://www.nrc.gov/docs/ML2005/ML20052H367.pdf
// This vessel was 

// In your test module

use uom::si::area::square_foot;
use uom::si::available_energy::btu_it_per_pound;
use uom::si::f64::*;
use uom::si::mass_flux::kilogram_per_square_meter_second;
use uom::si::mass_rate::pound_per_second;
use uom::si::ratio::ratio;
use uom::si::pressure::{atmosphere, kilopascal, megapascal, pound_force_per_square_inch};
use uom::si::length::millimeter;
use uom::si::velocity::meter_per_second;

use crate::interfaces::object_oriented_programming::TampinesSteamTableCV;
use crate::steam_turbine_equations::calculate_velocity_mass_flowrate_and_state_in_cd_nozzle;


/// From Figure 1 of:
///
/// Moody, F. J. (1975). Maximum discharge rate of liquid-vapor mixtures 
/// from vessels (No.
/// NEDO--21052). General Electric Co., San Jose, CA (United States). 
/// BWR Projects Dept..0 
///
/// Downloaded at:
/// https://www.osti.gov/servlets/purl/7309475
///
/// p_ref = 0.25
///"dimensionless stagnation enthalpy","dimensionless critical mass flux"
/// 0.4902,3.8593
/// 0.8039,3.7306
/// 1.1961,3.4469
/// 1.4314,3.1135
/// 1.6275,2.7187
/// 1.7647,2.2948
/// 1.8627,1.5106
/// 1.902,1.1133
/// 1.9412,0.4991
/// 1.9804,0.3678
/// 2.0588,0.2901
/// 2.2157,0.2212
/// 2.4118,0.1867
/// 2.8431,0.1541
/// 3.1176,0.1392
/// 3.5098,0.1243
/// 3.8431,0.111
/// 4.4902,0.1014
/// 5.0784,0.0916
/// 5.5098,0.0866
/// 6.1765,0.0791
/// 6.6863,0.0756
/// 7.2549,0.0706
/// 7.8039,0.0691
/// 8.3725,0.0653
/// 8.7255,0.0631
/// 9.3333,0.061
/// 9.902,0.057
/// 10.1765,0.057
/// 10.549,0.0564
/// 11.0784,0.0551
/// 11.5098,0.0533
/// 
///
#[test]
fn isobar_pref_0_25() {


    // let's first have our datapoints 

    let isobar_pref_0_25_critical_dimensionless_critical_mass_flux = vec![
        //(0.4902,3.8593),
        //(0.8039,3.7306),
        //(1.1961,3.4469),
        //(1.4314,3.1135),
        //(1.6275,2.7187),
        //(1.7647,2.2948),
        //(1.8627,1.5106),
        //(1.902,1.1133),
        //(1.9412,0.4991),
        //(1.9804,0.3678),
        //(2.0588,0.2901),
        //(2.2157,0.2212),
        //(2.4118,0.1867),
        //(2.8431,0.1541),
        //(3.1176,0.1392),
        //(3.5098,0.1243),
        //(3.8431,0.111),
        //(4.4902,0.1014),
        //(5.0784,0.0916),
        //(5.5098,0.0866),
        //(6.1765,0.0791),
        //(6.6863,0.0756),
        //(7.2549,0.0706),
        //(7.8039,0.0691),
        //(8.3725,0.0653),
        //(8.7255,0.0631),
        //(9.3333,0.061),
        //(9.902,0.057),
        //(10.1765,0.057),
        //(10.549,0.0564),
        //(11.0784,0.0551),
        (11.5098,0.0533),
        ];


    let p_ref = Pressure::new::<pound_force_per_square_inch>(100.0);
    let dimensionless_stagnation_pressure = 0.25;
    // for clarity: it is 2.326e5 joule/kg
    // it matches closer to btu_it as opposed to 
    // btu, which is:
    // 2.324_443_707_610_621_E3 joule/kg
    let h_ref = AvailableEnergy::new::<btu_it_per_pound>(100.0);
    let g_ref: MassFlux = 
        MassRate::new::<pound_per_second>(1000.0)/
        Area::new::<square_foot>(1.0);

    // ref vol 
    let ref_vol = TampinesSteamTableCV::get_ref_vol();

    for (h_dimensionless_ptr, g_dimensionless_ptr) 
        in isobar_pref_0_25_critical_dimensionless_critical_mass_flux.iter() {

            let h0: AvailableEnergy = h_ref * (*h_dimensionless_ptr);
            let p0: Pressure = dimensionless_stagnation_pressure * p_ref;

            let g_ref_expected: MassFlux = g_ref * (*g_dimensionless_ptr);

            let state_0 = TampinesSteamTableCV::new_from_ph(
                p0, h0, ref_vol
            );

            let g_test = state_0.get_stagnation_critical_mass_flux();
            // this helps see which point we are at on the graph
            dbg!(&(*h_dimensionless_ptr,*g_dimensionless_ptr));

            approx::assert_relative_eq!(
                g_ref_expected.get::<kilogram_per_square_meter_second>(),
                g_test.get::<kilogram_per_square_meter_second>(),
                max_relative=0.02
            );


    }

    todo!()

    



}

