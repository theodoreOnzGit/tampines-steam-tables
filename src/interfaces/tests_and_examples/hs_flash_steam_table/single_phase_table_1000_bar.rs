use uom::si::{dynamic_viscosity::micropascal_second, thermal_conductivity::milliwatt_per_meter_kelvin};
use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::pressure::{bar, kilopascal, megapascal};
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::velocity::meter_per_second;

use crate::interfaces::functional_programming::hs_flash_eqm::{cp_hs_eqm, kappa_hs_eqm, lambda_hs_eqm, mu_hs_eqm, tpvx_hs_flash_eqm, w_hs_eqm};
use crate::interfaces::functional_programming::ph_flash_eqm::x_ph_flash;

/// single phase table (see page 192)
///
#[test]
#[ignore = "At 1000 bar, ph flashing goes out of bounds for some reason, yet to debug"]
pub fn single_phase_table_2_to_750_degc(){

    let steam_table: Vec<[f64; 10]> =
        vec![
        //[1000.000,0.000,0.000956687,95.386,-0.008582,3.9057,1575.5000,25.9470,1660.6000,616.6000],
        [1000.000,2.000,0.000957105,103.207,0.019946,3.915,1584.6000,26.2350,1565.3000,619.8900],
        [1000.000,4.000,0.000957548,111.046,0.048331,3.9235,1593.2000,26.5090,1478.5000,623.1300],
        [1000.000,6.000,0.000958014,118.9,0.076571,3.9312,1601.4000,26.7700,1399.1000,626.3500],
        [1000.000,8.000,0.000958503,126.77,0.10466,3.9382,1609.3000,27.0200,1326.2000,629.5300],
        [1000.000,10.000,0.000959015,134.653,0.1326,3.9446,1616.8000,27.2570,1259.3000,632.6700],
        [1000.000,12.000,0.000959548,142.548,0.16038,3.9503,1624.0000,27.4840,1197.5000,635.7800],
        [1000.000,14.000,0.000960102,150.454,0.18801,3.9556,1630.8000,27.7000,1140.5000,638.8500],
        [1000.000,16.000,0.000960677,158.37,0.21549,3.9604,1637.3000,27.9060,1087.7000,641.8900],
        [1000.000,18.000,0.000961273,166.295,0.2428,3.9649,1643.6000,28.1030,1038.7000,644.8900],
        [1000.000,20.000,0.000961888,174.229,0.26996,3.9689,1649.6000,28.2900,993.1400,647.8600],
        [1000.000,25.000,0.00096351,194.096,0.33716,3.9777,1663.4000,28.7180,892.3600,655.1000],
        [1000.000,30.000,0.000965249,214.003,0.40337,3.9849,1675.8000,29.0940,807.1300,662.1000],
        [1000.000,35.000,0.000967101,233.943,0.46861,3.9909,1686.8000,29.4190,734.3700,668.8400],
        [1000.000,40.000,0.000969061,253.911,0.53289,3.9961,1696.4000,29.6980,671.7500,675.3100],
        [1000.000,45.000,0.000971126,273.903,0.59623,4.0007,1704.9000,29.9320,617.4600,681.4900],
        [1000.000,50.000,0.000973294,293.917,0.65864,4.0048,1712.2000,30.1220,570.0700,687.3900],
        [1000.000,55.000,0.000975562,313.95,0.72016,4.0086,1718.5000,30.2720,528.4700,693.0000],
        [1000.000,60.000,0.000977929,334.003,0.78081,4.0122,1723.7000,30.3830,491.7400,698.3100],
        [1000.000,65.000,0.000980392,354.072,0.84061,4.0157,1728.0000,30.4560,459.1600,703.3100],
        [1000.000,70.000,0.00098295,374.16,0.89957,4.0191,1731.3000,30.4940,430.1200,708.0100],
        [1000.000,75.000,0.000985603,394.264,0.95774,4.0225,1733.8000,30.4980,404.1400,712.4000],
        [1000.000,80.000,0.000988348,414.385,1.0151,4.0259,1735.4000,30.4710,380.7900,716.4900],
        [1000.000,85.000,0.000991185,434.523,1.0717,4.0293,1736.2000,30.4130,359.7500,720.2700],
        [1000.000,90.000,0.000994115,454.678,1.1276,4.0327,1736.3000,30.3270,340.7100,723.7500],
        [1000.000,95.000,0.000997135,474.85,1.1828,4.0362,1735.7000,30.2140,323.4300,726.9200],
        [1000.000,100.000,0.00100025,495.04,1.2373,4.0397,1734.4000,30.0750,307.7100,729.8000],
        [1000.000,110.000,0.00100674,535.473,1.3442,4.0471,1730.0000,29.7300,280.2300,734.6900],
        [1000.000,120.000,0.00101361,575.982,1.4486,4.0548,1723.4000,29.3020,257.1100,738.4400],
        [1000.000,130.000,0.00102084,616.571,1.5505,4.063,1714.7000,28.8030,237.4800,741.1100],
        [1000.000,140.000,0.00102844,657.245,1.6502,4.0718,1704.3000,28.2420,220.6600,742.7300],
        [1000.000,150.000,0.00103643,698.01,1.7477,4.0813,1692.1000,27.6270,206.1500,743.3500],
        [1000.000,160.000,0.00104482,738.873,1.8431,4.0916,1678.5000,26.9660,193.5400,743.0300],
        [1000.000,170.000,0.00105361,779.845,1.9366,4.1028,1663.6000,26.2660,182.5000,741.8000],
        [1000.000,180.000,0.00106282,820.934,2.0283,4.1152,1647.3000,25.5320,172.7800,739.7200],
        [1000.000,190.000,0.00107246,862.152,2.1183,4.1288,1629.8000,24.7680,164.1800,736.8200],
        [1000.000,200.000,0.00108256,903.513,2.2066,4.1437,1611.2000,23.9790,156.5100,733.1400],
        [1000.000,210.000,0.00109313,945.032,2.2935,4.1603,1591.5000,23.1690,149.6400,728.7200],
        [1000.000,220.000,0.00110421,986.725,2.3789,4.1785,1570.7000,22.3420,143.4600,723.9200],
        [1000.000,230.000,0.00111581,1028.61,2.463,4.1987,1548.9000,21.5010,137.8700,718.4200],
        [1000.000,240.000,0.00112796,1070.7,2.5458,4.2207,1526.2000,20.6490,132.7800,712.2100],
        [1000.000,250.000,0.00114071,1113.03,2.6275,4.2449,1502.5000,19.7910,128.1200,705.3700],
        [1000.000,260.000,0.00115408,1155.61,2.7081,4.2714,1478.1000,18.9300,123.8500,697.9400],
        [1000.000,270.000,0.00116813,1198.47,2.7878,4.3001,1452.8000,18.0690,119.9000,689.9400],
        [1000.000,280.000,0.00118289,1241.62,2.8665,4.3311,1426.9000,17.2120,116.2300,681.4100],
        [1000.000,290.000,0.00119843,1285.1,2.9444,4.3646,1400.4000,16.3640,112.8100,672.3600],
        [1000.000,300.000,0.0012148,1328.92,3.0215,4.4003,1373.4000,15.5270,109.6100,662.8200],
        [1000.000,310.000,0.00123207,1373.11,3.098,4.4382,1346.1000,14.7070,106.5900,652.8100],
        [1000.000,320.000,0.00125031,1417.69,3.1738,4.4781,1318.5000,13.9040,103.7300,642.3600],
        [1000.000,330.000,0.00126961,1462.68,3.249,4.5197,1290.8000,13.1240,101.0100,631.4700],
        [1000.000,340.000,0.00129006,1508.09,3.3236,4.5622,1263.0000,12.3660,98.4090,620.1800],
        [1000.000,350.000,0.00131176,1553.92,3.3978,4.6048,1235.2000,11.6320,95.9150,608.5000],
        [1000.000,360.000,0.00133481,1600.21,3.4715,4.6541,1207.4000,10.9220,93.5130,596.4700],
        [1000.000,370.000,0.00135934,1647.05,3.5449,4.7144,1179.3000,10.2320,91.1910,584.1100],
        [1000.000,380.000,0.00138548,1694.5,3.6181,4.7739,1150.7000,9.5566,88.9400,571.4500],
        [1000.000,390.000,0.00141337,1742.53,3.691,4.832,1121.9000,8.9058,86.7510,558.5200],
        [1000.000,400.000,0.00144317,1791.14,3.7638,4.8917,1093.3000,8.2826,84.6190,545.3300],
        [1000.000,410.000,0.00147508,1840.37,3.8364,4.9557,1065.0000,7.6892,82.5370,531.9300],
        [1000.000,420.000,0.00150929,1890.27,3.9089,5.0253,1037.2000,7.1273,80.5000,518.3300],
        [1000.000,430.000,0.00154605,1940.9,3.9814,5.1001,1010.0000,6.5981,78.5050,504.5600],
        [1000.000,440.000,0.00158559,1992.29,4.054,5.1784,983.7000,6.1028,76.5500,490.6300],
        [1000.000,450.000,0.00162815,2044.47,4.1267,5.2581,958.4300,5.6419,74.6370,476.6000],
        [1000.000,460.000,0.00167397,2097.44,4.1994,5.3363,934.3700,5.2154,72.7660,462.5000],
        [1000.000,470.000,0.00172324,2151.18,4.2722,5.4102,911.6600,4.8230,70.9450,448.3800],
        [1000.000,480.000,0.00177614,2205.62,4.345,5.4766,890.4200,4.4639,69.1770,434.3000],
        [1000.000,490.000,0.00183278,2260.68,4.4176,5.5326,870.7500,4.1369,67.4710,420.3400],
        [1000.000,500.000,0.00189324,2316.23,4.4899,5.5757,852.7300,3.8407,65.8330,406.5600],
        [1000.000,510.000,0.0019575,2372.14,4.5618,5.6038,836.4200,3.5740,64.2720,393.0500],
        [1000.000,520.000,0.00202548,2428.25,4.633,5.6155,821.8700,3.3348,62.7940,379.8800],
        [1000.000,530.000,0.00209703,2484.39,4.7033,5.6101,809.0700,3.1215,61.4040,367.1400],
        [1000.000,540.000,0.00217192,2540.4,4.7726,5.5877,798.0300,2.9322,60.1080,354.8900],
        [1000.000,550.000,0.00224984,2596.09,4.8407,5.549,788.6900,2.7647,58.9080,343.2000],
        [1000.000,560.000,0.00233047,2651.33,4.9074,5.4954,780.9700,2.6172,57.8060,332.1200],
        [1000.000,570.000,0.00241341,2705.96,4.9726,5.4288,774.8000,2.4874,56.8000,321.6900],
        [1000.000,580.000,0.00249829,2759.87,5.0361,5.3512,770.0500,2.3735,55.8890,311.9400],
        [1000.000,590.000,0.00258472,2812.94,5.098,5.2687,769.6900,2.2920,55.0700,302.8800],
        [1000.000,600.000,0.00267226,2865.07,5.158,5.1706,766.5300,2.1988,54.3390,294.5100],
        [1000.000,650.000,0.00311448,3110.6,5.4316,4.6275,767.3100,1.8904,51.8250,262.5600],
        [1000.000,700.000,0.00354616,3330.76,5.664,4.191,784.0800,1.7336,50.7110,242.8900],
        [1000.000,750.000,0.00395319,3530.68,5.8644,3.8235,801.5000,1.6250,50.4780,234.8400],
        //[1000.000,800.000,0.00433551,3715.19,6.0405,3.5762,821.0000,1.5547,50.7810,232.2400],

        ];

        for dataset in steam_table {
            let p_bar = dataset[0];
            let t_deg_c = dataset[1];
            let v_m3_per_kg = dataset[2];
            let h_kj_per_kg = dataset[3];
            let s_kj_per_kg_k = dataset[4];
            let cp_kj_per_kg_k = dataset[5];
            let w_m_per_s = dataset[6];
            let kappa_dimensionless = dataset[7];
            let eta_micropascal_second = dataset[8];
            let lambda_milliwatt_per_meter_kelvin = dataset[9];
            assert_hs_flash(p_bar, t_deg_c, v_m3_per_kg, h_kj_per_kg, 
                s_kj_per_kg_k, cp_kj_per_kg_k, w_m_per_s, 
                kappa_dimensionless, eta_micropascal_second, 
                lambda_milliwatt_per_meter_kelvin);
        }

}

fn assert_hs_flash(
    p_bar: f64,
    t_deg_c: f64,
    v_m3_per_kg: f64,
    h_kj_per_kg: f64,
    s_kj_per_kg_k: f64,
    cp_kj_per_kg_k: f64,
    w_m_per_s: f64,
    kappa_dimensionless: f64,
    eta_micropascal_second: f64,
    lambda_milliwatt_per_meter_kelvin: f64,
){
    let p = Pressure::new::<bar>(p_bar);
    let h = AvailableEnergy::new::<kilojoule_per_kilogram>(h_kj_per_kg);
    let s = SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(
        s_kj_per_kg_k
    );

    // assert temp first to within 0.025 mk 
    // based on table 2.8
    let temp_tol_millikelvin = 25.0;

    let (t_test, p_test, v_test, x_test) = 
        tpvx_hs_flash_eqm(h, s);

    approx::assert_abs_diff_eq!(
        t_deg_c,
        t_test.get::<degree_celsius>(),
        epsilon=temp_tol_millikelvin*1e-3
        );

    // assert volume to within 0.5%  
    approx::assert_relative_eq!(
        v_m3_per_kg,
        v_test.get::<cubic_meter_per_kilogram>(),
        max_relative=5e-3
        );

    // assert pressure to within 15 kPa or 0.60%
    if p > Pressure::new::<megapascal>(2.5) {
        let pressure_error_tol = Pressure::new::<kilopascal>(15.0);
        let pressure_error_tol_bar = pressure_error_tol.get::<bar>();
        approx::assert_abs_diff_eq!(
            p_bar,
            p_test.get::<bar>(),
            epsilon=pressure_error_tol_bar
        ); 
    } else {
        dbg!(&(p_bar,t_deg_c));
        approx::assert_relative_eq!(
            p_bar,
            p_test.get::<bar>(),
            max_relative = 0.60e-2
        ); 

    }
    // the x for ph is quite reliable, i'll use that as reference 
    let x_ref = x_ph_flash(p, h);

    approx::assert_relative_eq!(
        x_ref,
        x_test.get::<ratio>(),
        max_relative = 1e-3
    );


    // cp 
    let cp_test = cp_hs_eqm(h, s);
    approx::assert_relative_eq!(
        cp_kj_per_kg_k,
        cp_test.get::<kilojoule_per_kilogram_kelvin>(),
        max_relative=1e-4
        );
    // w 
    let w_test = w_hs_eqm(h, s);
    approx::assert_relative_eq!(
        w_m_per_s,
        w_test.get::<meter_per_second>(),
        max_relative=1e-4
        );

    // kappa
    let kappa_test = kappa_hs_eqm(h, s);
    approx::assert_relative_eq!(
        kappa_dimensionless,
        kappa_test.get::<ratio>(),
        max_relative=5e-3
        );

    // dynamic_viscosity
    //
    let eta_micropascal_second_test = mu_hs_eqm(h, s)
        .get::<micropascal_second>();
    approx::assert_relative_eq!(
        eta_micropascal_second,
        eta_micropascal_second_test,
        max_relative=2e-2
        );

    // thermal thermal conductivity
    let lambda_test_milliwatt_per_meter_kelvin = 
        lambda_hs_eqm(h, s).get::<milliwatt_per_meter_kelvin>();
    approx::assert_relative_eq!(
        lambda_milliwatt_per_meter_kelvin,
        lambda_test_milliwatt_per_meter_kelvin,
        max_relative=1e-2
        );

}



