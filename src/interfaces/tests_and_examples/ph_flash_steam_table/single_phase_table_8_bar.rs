use uom::si::thermal_conductivity::milliwatt_per_meter_kelvin;
use uom::si::{available_energy::kilojoule_per_kilogram, dynamic_viscosity::micropascal_second};
use uom::si::pressure::bar;
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::velocity::meter_per_second;

use crate::dynamic_viscosity::mu_ph_eqm;
use crate::interfaces::functional_programming::ph_flash_eqm::{cp_ph_eqm, kappa_ph_eqm, lambda_ph_eqm, s_ph_eqm, t_ph_eqm, v_ph_eqm, w_ph_eqm};

/// single phase table (see page 201)
///
/// NOTE: ph flash UNABLE to do triple point liquid and vapour accurately.
/// or even close to 0 degc
#[test]
pub fn single_phase_table_2_to_800_degc(){

    let steam_table: Vec<[f64; 10]> =
        vec![
        //[8.000,0.000,0.000999801,0.77234,-0.0001013,4.216,1403.5000,2462.9000,1790.2000,556.1800],
        [8.000,2.000,0.000999706,9.19771,0.030632,4.2096,1413.3000,2497.6000,1672.2000,561.1800],
        [8.000,4.000,0.000999679,17.6115,0.0611,4.2043,1422.7000,2531.0000,1566.2000,565.9600],
        [8.000,6.000,0.000999716,26.0154,0.091314,4.1998,1431.8000,2563.1000,1470.6000,570.5600],
        [8.000,8.000,0.000999812,34.4111,0.12128,4.196,1440.4000,2593.9000,1384.0000,574.9800],
        [8.000,10.000,0.000999965,42.7999,0.15101,4.1928,1448.7000,2623.4000,1305.3000,579.2300],
        [8.000,12.000,0.00100017,51.1827,0.18052,4.1901,1456.5000,2651.5000,1233.5000,583.3400],
        [8.000,14.000,0.00100043,59.5605,0.20979,4.1878,1464.1000,2678.2000,1167.9000,587.3000],
        [8.000,16.000,0.00100073,67.9339,0.23885,4.1858,1471.2000,2703.7000,1107.7000,591.1300],
        [8.000,18.000,0.00100108,76.3037,0.2677,4.1841,1478.0000,2727.8000,1052.4000,594.8300],
        [8.000,20.000,0.00100148,84.6704,0.29634,4.1826,1484.5000,2750.7000,1001.4000,598.4200],
        [8.000,25.000,0.00100264,105.576,0.36705,4.1799,1499.3000,2802.4000,889.9300,606.9100],
        [8.000,30.000,0.00100405,126.471,0.43655,4.1781,1512.1000,2846.5000,797.2100,614.7800],
        [8.000,35.000,0.00100569,147.359,0.50489,4.1772,1523.1000,2883.4000,719.1700,622.0800],
        [8.000,40.000,0.00100753,168.243,0.57212,4.1768,1532.4000,2913.5000,652.8200,628.8600],
        [8.000,45.000,0.00100956,189.128,0.63829,4.1771,1540.2000,2937.1000,595.8900,635.1600],
        [8.000,50.000,0.00101179,210.015,0.70343,4.1779,1546.5000,2954.7000,546.6600,641.0000],
        [8.000,55.000,0.00101419,230.908,0.76759,4.1793,1551.4000,2966.5000,503.7900,646.4000],
        [8.000,60.000,0.00101676,251.809,0.8308,4.1812,1555.1000,2972.9000,466.2100,651.3800],
        [8.000,65.000,0.0010195,272.721,0.89311,4.1836,1557.5000,2974.3000,433.0900,655.9600],
        [8.000,70.000,0.0010224,293.647,0.95453,4.1866,1558.8000,2970.9000,403.7400,660.1400],
        [8.000,75.000,0.00102547,314.588,1.0151,4.19,1559.1000,2963.0000,377.6100,663.9500],
        [8.000,80.000,0.00102868,335.548,1.0749,4.194,1558.4000,2951.0000,354.2500,667.3800],
        [8.000,85.000,0.00103206,356.529,1.1339,4.1985,1556.7000,2935.0000,333.2700,670.4600],
        [8.000,90.000,0.00103559,377.533,1.1921,4.2035,1554.1000,2915.4000,314.3700,673.1800],
        [8.000,95.000,0.00103927,398.564,1.2496,4.209,1550.7000,2892.4000,297.2800,675.5700],
        [8.000,100.000,0.0010431,419.624,1.3065,4.215,1546.5000,2866.1000,281.7700,677.6100],
        [8.000,110.000,0.00105123,461.841,1.4181,4.2288,1535.9000,2804.9000,254.7900,680.7300],
        [8.000,120.000,0.00105999,504.207,1.5273,4.2448,1522.4000,2733.1000,232.1900,682.6100],
        [8.000,130.000,0.0010694,546.746,1.6341,4.2634,1506.3000,2652.1000,213.0800,683.2800],
        [8.000,140.000,0.00107948,589.484,1.7388,4.2847,1487.7000,2563.0000,196.7600,682.8100],
        [8.000,150.000,0.00109029,632.451,1.8416,4.3092,1466.9000,2466.9000,182.6900,681.2400],
        [8.000,160.000,0.00110186,675.681,1.9426,4.3373,1443.7000,2364.5000,170.4800,678.8500],
        [8.000,170.000,0.00111426,719.211,2.0419,4.3695,1418.3000,2256.6000,159.7800,675.5300],
        [8.000,170.414,0.00111479,721.018,2.046,4.3709,1417.2000,2252.0000,159.3600,675.3700],
        [8.000,170.414,0.240328,2768.3,6.6615,2.6032,498.8700,1.2944,14.6590,33.2850],
        [8.000,180.000,0.247183,2792.44,6.7154,2.4499,506.9000,1.2994,15.0790,34.0220],
        [8.000,190.000,0.254105,2816.46,6.7678,2.3621,514.3400,1.3013,15.5160,34.8130],
        [8.000,200.000,0.260868,2839.77,6.8176,2.303,521.2900,1.3021,15.9500,35.6240],
        [8.000,210.000,0.267506,2862.57,6.8653,2.2586,527.9500,1.3024,16.3820,36.4550],
        [8.000,220.000,0.274043,2884.97,6.9112,2.2236,534.3600,1.3025,16.8140,37.3060],
        [8.000,230.000,0.280496,2907.06,6.9555,2.1956,540.5800,1.3023,17.2440,38.1760],
        [8.000,240.000,0.286878,2928.9,6.9985,2.1732,546.6300,1.3020,17.6730,39.0640],
        [8.000,250.000,0.293199,2950.54,7.0403,2.1553,552.5200,1.3015,18.1010,39.9710],
        [8.000,260.000,0.299469,2972.02,7.081,2.1411,558.2800,1.3010,18.5280,40.8950],
        [8.000,270.000,0.305693,2993.37,7.1206,2.1298,563.9100,1.3003,18.9540,41.8360],
        [8.000,280.000,0.311878,3014.63,7.1594,2.1211,569.4200,1.2996,19.3800,42.7940],
        [8.000,290.000,0.31803,3035.8,7.1974,2.1145,574.8400,1.2988,19.8050,43.7680],
        [8.000,300.000,0.324151,3056.92,7.2345,2.1097,580.1500,1.2979,20.2290,44.7580],
        [8.000,310.000,0.330245,3078.0,7.271,2.1063,585.3700,1.2970,20.6530,45.7620],
        [8.000,320.000,0.336316,3099.05,7.3068,2.1042,590.5100,1.2961,21.0760,46.7820],
        [8.000,330.000,0.342366,3120.09,7.342,2.1031,595.5800,1.2951,21.4980,47.8150],
        [8.000,340.000,0.348397,3141.12,7.3765,2.1029,600.5700,1.2941,21.9200,48.8620],
        [8.000,350.000,0.354411,3162.15,7.4106,2.1035,605.4900,1.2930,22.3410,49.9220],
        [8.000,360.000,0.360409,3183.19,7.4441,2.1048,610.3400,1.2920,22.7620,50.9950],
        [8.000,370.000,0.366393,3204.25,7.4771,2.1066,615.1400,1.2909,23.1820,52.0810],
        [8.000,380.000,0.372365,3225.33,7.5096,2.109,619.8700,1.2899,23.6010,53.1780],
        [8.000,390.000,0.378324,3246.43,7.5416,2.1117,624.5500,1.2888,24.0190,54.2860],
        [8.000,400.000,0.384273,3267.56,7.5733,2.1149,629.1700,1.2877,24.4370,55.4060],
        [8.000,410.000,0.390212,3288.73,7.6045,2.1184,633.7400,1.2866,24.8550,56.5370],
        [8.000,420.000,0.396142,3309.93,7.6353,2.1223,638.2600,1.2855,25.2710,57.6780],
        [8.000,430.000,0.402063,3331.17,7.6657,2.1264,642.7400,1.2844,25.6870,58.8290],
        [8.000,440.000,0.407977,3352.46,7.6958,2.1307,647.1700,1.2832,26.1020,59.9890],
        [8.000,450.000,0.413883,3373.79,7.7255,2.1352,651.5500,1.2821,26.5170,61.1590],
        [8.000,460.000,0.419783,3395.16,7.7548,2.14,655.8900,1.2810,26.9300,62.3390],
        [8.000,470.000,0.425676,3416.59,7.7839,2.1449,660.1900,1.2799,27.3430,63.5270],
        [8.000,480.000,0.431563,3438.06,7.8126,2.15,664.4500,1.2788,27.7550,64.7240],
        [8.000,490.000,0.437444,3459.59,7.841,2.1553,668.6700,1.2777,28.1670,65.9290],
        [8.000,500.000,0.443321,3481.17,7.869,2.1606,672.8600,1.2766,28.5770,67.1420],
        [8.000,510.000,0.449192,3502.8,7.8969,2.1661,677.0000,1.2754,28.9870,68.3630],
        [8.000,520.000,0.455059,3524.49,7.9244,2.1717,681.1200,1.2743,29.3950,69.5920],
        [8.000,530.000,0.460921,3546.24,7.9516,2.1774,685.1900,1.2732,29.8030,70.8280],
        [8.000,540.000,0.46678,3568.04,7.9786,2.1833,689.2400,1.2721,30.2100,72.0710],
        [8.000,550.000,0.472635,3589.9,8.0053,2.1892,693.2500,1.2711,30.6170,73.3210],
        [8.000,560.000,0.478486,3611.82,8.0318,2.1951,697.2300,1.2700,31.0220,74.5770],
        [8.000,570.000,0.484333,3633.81,8.058,2.2012,701.1800,1.2689,31.4260,75.8410],
        [8.000,580.000,0.490178,3655.85,8.084,2.2073,705.1000,1.2678,31.8300,77.1100],
        [8.000,590.000,0.496019,3677.95,8.1098,2.2135,708.9900,1.2667,32.2330,78.3860],
        [8.000,600.000,0.501857,3700.12,8.1353,2.2197,712.8500,1.2657,32.6340,79.6680],
        [8.000,650.000,0.531012,3811.9,8.2598,2.2516,731.7500,1.2605,34.6280,86.1620],
        [8.000,700.000,0.560113,3925.29,8.3794,2.2843,750.0400,1.2555,36.5990,92.7840],
        [8.000,750.000,0.589174,4040.33,8.4947,2.3175,767.7700,1.2506,38.5440,99.5180],
        [8.000,800.000,0.618202,4157.04,8.606,2.351,785.0000,1.2460,40.4640,106.3500],

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
            assert_ph_flash(p_bar, t_deg_c, v_m3_per_kg, h_kj_per_kg, 
                s_kj_per_kg_k, cp_kj_per_kg_k, w_m_per_s, 
                kappa_dimensionless, eta_micropascal_second, 
                lambda_milliwatt_per_meter_kelvin);
        }
}

fn assert_ph_flash(
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

    // assert temp first to within 0.035 mk 
    let temp_tol_millikelvin = 35.0;

    let t_test = t_ph_eqm(p, h);

    approx::assert_abs_diff_eq!(
        t_deg_c,
        t_test.get::<degree_celsius>(),
        epsilon=temp_tol_millikelvin*1e-3
        );

    // assert volume to within 0.01%  (that's the tolerable error for 
    // backward eqn)
    let v_test = v_ph_eqm(p, h);
    approx::assert_relative_eq!(
        v_m3_per_kg,
        v_test.get::<cubic_meter_per_kilogram>(),
        max_relative=1e-4
        );

    // now entropy 
    let s_test = s_ph_eqm(p, h);
    approx::assert_relative_eq!(
        s_kj_per_kg_k,
        s_test.get::<kilojoule_per_kilogram_kelvin>(),
        max_relative=7e-3
        );

    // cp 
    let cp_test = cp_ph_eqm(p, h);
    approx::assert_relative_eq!(
        cp_kj_per_kg_k,
        cp_test.get::<kilojoule_per_kilogram_kelvin>(),
        max_relative=1e-4
        );
    // w 
    let w_test = w_ph_eqm(p, h);
    approx::assert_relative_eq!(
        w_m_per_s,
        w_test.get::<meter_per_second>(),
        max_relative=1e-4
        );

    // kappa
    let kappa_test = kappa_ph_eqm(p, h);
    approx::assert_relative_eq!(
        kappa_dimensionless,
        kappa_test.get::<ratio>(),
        max_relative=6e-3
        );

    // dynamic_viscosity
    //
    let eta_micropascal_second_test = mu_ph_eqm(p, h)
        .get::<micropascal_second>();
    approx::assert_relative_eq!(
        eta_micropascal_second,
        eta_micropascal_second_test,
        max_relative=2e-2
        );

    // thermal conductivity
    let lambda_test_milliwatt_per_meter_kelvin = 
        lambda_ph_eqm(p, h).get::<milliwatt_per_meter_kelvin>();
    approx::assert_relative_eq!(
        lambda_milliwatt_per_meter_kelvin,
        lambda_test_milliwatt_per_meter_kelvin,
        max_relative=1e-2
        );



}



