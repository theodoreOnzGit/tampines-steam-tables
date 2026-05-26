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
        (11.0784,0.0551),
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
            dbg!(&(*g_dimensionless_ptr,g_test/g_ref));


            approx::assert_relative_eq!(
                g_ref_expected.get::<kilogram_per_square_meter_second>(),
                g_test.get::<kilogram_per_square_meter_second>(),
                max_relative=0.02
            );


    }

    todo!()


}
// For p0/p_ref = 0.50 
// "x","y"
// 0.4902,5.4168
// 0.7647,5.2362
// 1.2353,5.0617
// 1.6471,4.6241
// 1.9412,3.9031
// 2.1765,3.1135
// 2.2549,2.269
// 2.3137,1.275
// 2.4118,0.7005
// 2.6078,0.4508
// 3,0.336
// 3.4314,0.2773
// 3.9804,0.2314
// 4.5686,0.2021
// 5.1373,0.1867
// 5.8235,0.1668
// 6.2157,0.1612
// 7.1569,0.144
// 8.1373,0.133
// 8.8627,0.1271
// 9.8431,0.1148
// 10.5686,0.111
// 11.2549,0.1073
// 11.7255,0.1037

// For p0/p_ref = 1.00
// "x","y"
// 0.451,7.6029
// 0.6667,7.6029
// 1.3137,7.1852
// 1.9412,6.5641
// 2.4314,5.0617
// 2.6471,3.9475
// 2.8235,2.7496
// 2.8627,1.9812
// 2.9216,1.205
// 3.0588,0.8393
// 3.5686,0.598
// 4.2941,0.4559
// 5.1373,0.3762
// 6.0392,0.3248
// 6.6275,0.3001
// 7.451,0.2742
// 8.2549,0.2591
// 9.2353,0.2394
// 10.1373,0.2263
// 11.1961,0.2138
// 11.8235,0.2067

// For p0/p_ref = 2.0
// "x","y"
// 0.4706,11.0394
// 0.9216,10.7927
// 1.3529,10.4329
// 1.9412,10.1997
// 2.5294,9.0074
// 2.8824,7.2669
// 3.1569,5.7968
// 3.3333,3.9923
// 3.4314,2.4009
// 3.4902,1.6723
// 3.9216,1.1781
// 4.5098,0.9397
// 5.2353,0.7843
// 6.0784,0.6695
// 7.3725,0.5716
// 8.9804,0.4879
// 10.4314,0.4358
// 11.902,0.3981

// For p0/p_ref = 4.0 
// "x","y"
// 0.7255,13.2273
// 1.1961,12.7864
// 1.6078,12.7864
// 2.0588,12.3602
// 2.5294,11.9481
// 2.8039,11.1648
// 3.1176,9.6394
// 3.3529,8.4169
// 3.5882,6.5641
// 3.7255,4.7298
// 3.7843,3.332
// 3.902,2.269
// 4.2745,1.7105
// 4.6471,1.4602
// 5.1373,1.2325
// 5.7255,1.1008
// 6.451,0.9832
// 7.451,0.8393
// 8.4118,0.7581
// 9.1569,0.7005
// 10.098,0.6771
// 10.9804,0.6327
// 11.8235,0.6256

// For p0/p_ref = 6.0 
//
// "x","y"
// 0.5882,18.5657
// 0.8235,18.5657
// 1.1765,18.5657
// 1.5098,18.3571
// 1.9608,17.9468
// 2.2941,17.7452
// 2.6275,16.9609
// 3.0196,16.3955
// 3.5882,14.1553
// 4.1373,10.6714
// 4.3922,8.0449
// 4.5098,5.7317
// 4.5686,4.321
// 4.7647,3.6062
// 5.098,2.9759
// 5.4902,2.5693
// 5.9216,2.269
// 6.6078,2.0037
// 7.2353,1.8305
// 7.6863,1.7105
// 8.4314,1.5451
// 9.0196,1.4768
// 9.7647,1.3645
// 10.549,1.275
// 11.2549,1.2187
// 11.7843,1.1517

// For p0/p_ref = 8.0
// "x","y"
// 0.6471,21.9954
// 0.8824,21.7482
// 1.451,21.2622
// 2.0196,20.787
// 2.5098,20.0941
// 3,19.2059
// 3.451,17.9468
// 3.8824,16.0291
// 4.2549,13.9963
// 4.4706,11.8139
// 4.6863,9.424
// 4.8039,7.6029
// 4.902,5.8627
// 5.0784,4.7298
// 5.3725,3.9031
// 5.7451,3.486
// 6.1961,3.0439
// 6.5686,2.7808
// 7.0196,2.5404
// 7.5294,2.3739
// 7.9412,2.2435
// 8.2745,2.1202
// 8.8039,2.0495
// 9.2941,1.9152
// 9.7059,1.8513
// 10.1765,1.7496
// 10.6471,1.73
// 11.0784,1.6165
// 11.6078,1.5804
// 11.9608,1.5804

