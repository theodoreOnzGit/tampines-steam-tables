use uom::si::available_energy::kilojoule_per_kilogram;
use uom::si::pressure::bar;
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::velocity::meter_per_second;

use crate::interfaces::functional_programming::ph_flash_eqm::{cp_ph_eqm, kappa_ph_eqm, s_ph_eqm, t_ph_eqm, v_ph_eqm, w_ph_eqm};

/// single phase table (see page 201)
///
/// NOTE: ph flash UNABLE to do triple point liquid and vapour accurately.
/// or even close to 0 degc
#[test]
pub fn single_phase_table_2_to_750_degc(){

    let steam_table: Vec<[f64; 10]> =
        vec![
        //[180.000,0.000,0.000991311,18.0523,0.0004737,4.1373,1431.20000,114.80000,1755.00000,568.75000],
        [180.000,2.000,0.000991335,26.3252,0.03065,4.1357,1441.00000,116.37000,1642.80000,573.32000],
        [180.000,4.000,0.000991416,34.5952,0.060598,4.1344,1450.50000,117.89000,1541.60000,577.73000],
        [180.000,6.000,0.000991551,42.8632,0.090323,4.1336,1459.50000,119.34000,1450.10000,581.98000],
        [180.000,8.000,0.000991738,51.1298,0.11983,4.133,1468.10000,120.73000,1367.00000,586.09000],
        [180.000,10.000,0.000991972,59.3954,0.14913,4.1327,1476.30000,122.06000,1291.20000,590.07000],
        [180.000,12.000,0.000992253,67.6605,0.17821,4.1325,1484.20000,123.33000,1222.00000,593.92000],
        [180.000,14.000,0.000992576,75.9254,0.2071,4.1324,1491.70000,124.54000,1158.50000,597.66000],
        [180.000,16.000,0.000992942,84.1903,0.23578,4.1325,1498.80000,125.69000,1100.20000,601.28000],
        [180.000,18.000,0.000993347,92.4553,0.26426,4.1326,1505.60000,126.78000,1046.50000,604.81000],
        [180.000,20.000,0.000993791,100.721,0.29256,4.1328,1512.10000,127.81000,996.84000,608.23000],
        [180.000,25.000,0.00099506,121.386,0.36246,4.1335,1526.80000,130.16000,888.07000,616.38000],
        [180.000,30.000,0.000996541,142.056,0.43121,4.1344,1539.70000,132.17000,797.24000,623.99000],
        [180.000,35.000,0.000998219,162.731,0.49885,4.1355,1550.90000,133.86000,720.56000,631.11000],
        [180.000,40.000,0.00100008,183.412,0.56543,4.1369,1560.40000,135.26000,655.17000,637.76000],
        [180.000,45.000,0.00100212,204.1,0.63097,4.1384,1568.50000,136.38000,598.94000,643.98000],
        [180.000,50.000,0.00100433,224.796,0.69552,4.1403,1575.10000,137.24000,550.22000,649.78000],
        [180.000,55.000,0.0010067,245.503,0.7591,4.1424,1580.40000,137.84000,507.70000,655.18000],
        [180.000,60.000,0.00100922,266.221,0.82176,4.1449,1584.50000,138.21000,470.38000,660.18000],
        [180.000,65.000,0.00101189,286.952,0.88353,4.1476,1587.50000,138.36000,437.44000,664.81000],
        [180.000,70.000,0.00101471,307.698,0.94443,4.1508,1589.30000,138.30000,408.21000,669.06000],
        [180.000,75.000,0.00101768,328.461,1.0045,4.1543,1590.20000,138.04000,382.16000,672.95000],
        [180.000,80.000,0.00102079,349.241,1.0638,4.1582,1590.10000,137.60000,358.84000,676.49000],
        [180.000,85.000,0.00102404,370.043,1.1223,4.1624,1589.00000,136.99000,337.89000,679.67000],
        [180.000,90.000,0.00102742,390.866,1.18,4.167,1587.20000,136.22000,319.00000,682.52000],
        [180.000,95.000,0.00103095,411.714,1.237,4.1721,1584.50000,135.29000,301.90000,685.04000],
        [180.000,100.000,0.00103462,432.587,1.2933,4.1775,1581.10000,134.23000,286.38000,687.23000],
        [180.000,110.000,0.00104236,474.421,1.404,4.1896,1572.10000,131.72000,259.34000,690.67000],
        [180.000,120.000,0.00105068,516.385,1.5121,4.2035,1560.40000,128.74000,236.69000,692.89000],
        [180.000,130.000,0.00105957,558.497,1.6178,4.2192,1546.30000,125.37000,217.51000,693.94000],
        [180.000,140.000,0.00106906,600.777,1.7214,4.2371,1530.00000,121.64000,201.13000,693.88000],
        [180.000,150.000,0.00107918,643.247,1.823,4.2574,1511.50000,117.61000,187.01000,692.74000],
        [180.000,160.000,0.00108996,685.933,1.9227,4.2803,1491.00000,113.31000,174.75000,690.57000],
        [180.000,170.000,0.00110145,728.864,2.0207,4.3063,1468.50000,108.77000,164.03000,687.67000],
        [180.000,180.000,0.00111371,772.071,2.1171,4.3358,1444.10000,104.03000,154.57000,683.89000],
        [180.000,190.000,0.00112679,815.593,2.2121,4.3693,1417.80000,99.11300,146.17000,679.11000],
        [180.000,200.000,0.00114077,859.473,2.3058,4.4075,1389.60000,94.03900,138.66000,673.41000],
        [180.000,210.000,0.00115575,903.761,2.3985,4.4511,1359.40000,88.83600,131.89000,666.82000],
        [180.000,220.000,0.00117185,948.516,2.4902,4.5011,1327.30000,83.52700,125.75000,659.37000],
        [180.000,230.000,0.0011892,993.807,2.5811,4.5585,1293.30000,78.13500,120.13000,651.06000],
        [180.000,240.000,0.00120796,1039.72,2.6714,4.6248,1257.10000,72.68500,114.95000,641.89000],
        [180.000,250.000,0.00122836,1086.34,2.7614,4.7018,1218.90000,67.19800,110.10000,631.87000],
        [180.000,260.000,0.00125066,1133.8,2.8513,4.7919,1178.50000,61.69400,105.60000,620.97000],
        [180.000,270.000,0.0012752,1182.23,2.9413,4.8984,1135.70000,56.19100,101.30000,609.17000],
        [180.000,280.000,0.00130243,1231.83,3.0317,5.026,1090.20000,50.69800,97.26000,596.44000],
        [180.000,290.000,0.00133297,1282.84,3.1231,5.1816,1041.70000,45.22300,93.30000,582.69000],
        [180.000,300.000,0.00136769,1335.6,3.216,5.3761,989.47000,39.76900,89.41000,567.84000],
        [180.000,310.000,0.00140788,1390.56,3.311,5.6275,932.73000,34.33000,85.54000,551.72000],
        [180.000,320.000,0.00145557,1448.44,3.4095,5.9683,870.13000,28.89800,81.59000,534.09000],
        [180.000,330.000,0.00151428,1510.43,3.5131,6.4631,800.55000,23.51200,77.44000,514.54000],
        [180.000,340.000,0.00159084,1578.71,3.6253,7.2691,722.19000,18.21400,72.89000,492.33000],
        [180.000,350.000,0.00170295,1658.65,3.7546,8.9989,619.19000,12.50800,67.45000,466.22000],
        [180.000,356.992,0.00183949,1732.02,3.8717,12.84,513.11000,7.95150,62.12000,445.07000],
        [180.000,356.992,0.00749867,2509.53,5.1055,22.966,410.33000,1.24740,24.96000,172.52000],
        [180.000,360.000,0.00810999,2566.03,5.195,15.82,428.85000,1.25990,24.72000,145.53000],
        [180.000,370.000,0.0094513,2683.67,5.3795,9.327,465.42000,1.27330,24.64000,114.66000],
        [180.000,380.000,0.0104189,2764.89,5.5048,7.1712,489.69000,1.27860,24.87000,102.90000],
        [180.000,390.000,0.0112174,2830.24,5.6041,5.9991,509.01000,1.28320,25.20000,96.49700],
        [180.000,400.000,0.0119147,2886.31,5.6881,5.2645,525.21000,1.28620,25.57000,92.55300],
        [180.000,410.000,0.0125434,2936.27,5.7618,4.7572,539.29000,1.28810,25.96000,89.96400],
        [180.000,420.000,0.013122,2981.89,5.8281,4.3826,551.82000,1.28920,26.37000,88.22100],
        [180.000,430.000,0.0136617,3024.21,5.8887,4.0939,563.20000,1.28990,26.78000,87.05700],
        [180.000,440.000,0.0141706,3063.96,5.9448,3.8649,573.67000,1.29020,27.20000,86.31300],
        [180.000,450.000,0.0146541,3101.65,5.9973,3.6792,583.40000,1.29030,27.62000,85.88600],
        [180.000,460.000,0.0151165,3137.66,6.0468,3.5262,592.52000,1.29030,28.04000,85.70900],
        [180.000,470.000,0.015561,3172.26,6.0936,3.3983,601.13000,1.29010,28.46000,85.73200],
        [180.000,480.000,0.0159901,3205.69,6.1383,3.2902,609.29000,1.28980,28.88000,85.92000],
        [180.000,490.000,0.016406,3238.12,6.1811,3.1981,617.07000,1.28940,29.29000,86.24700],
        [180.000,500.000,0.0168102,3269.69,6.2222,3.119,624.51000,1.28890,29.71000,86.69100],
        [180.000,510.000,0.0172042,3300.53,6.2618,3.0506,631.65000,1.28840,30.12000,87.23700],
        [180.000,520.000,0.017589,3330.73,6.3002,2.9911,638.52000,1.28780,30.53000,87.87100],
        [180.000,530.000,0.0179656,3360.38,6.3373,2.9391,645.16000,1.28710,30.95000,88.58300],
        [180.000,540.000,0.018335,3389.54,6.3734,2.8935,651.58000,1.28640,31.35000,89.36300],
        [180.000,550.000,0.0186977,3418.27,6.4085,2.8533,657.80000,1.28570,31.76000,90.20500],
        [180.000,560.000,0.0190544,3446.62,6.4427,2.8179,663.85000,1.28490,32.17000,91.10200],
        [180.000,570.000,0.0194055,3474.64,6.4762,2.7864,669.73000,1.28410,32.57000,92.04900],
        [180.000,580.000,0.0197517,3502.36,6.5089,2.7585,675.47000,1.28330,32.97000,93.04000],
        [180.000,590.000,0.0200933,3529.82,6.5408,2.7336,681.06000,1.28250,33.37000,94.07200],
        [180.000,600.000,0.0204306,3557.04,6.5722,2.7114,686.53000,1.28160,33.77000,95.14200],
        [180.000,650.000,0.0220638,3690.42,6.7208,2.6313,712.27000,1.27740,35.75000,100.95000],
        [180.000,700.000,0.0236279,3820.74,6.8583,2.5859,735.86000,1.27320,37.68000,107.31000],
        [180.000,750.000,0.0251418,3949.37,6.9872,2.562,757.84000,1.26910,39.59000,114.12000],
        //[180.000,800.000,0.0266179,4077.18,7.1091,2.5525,778.51000,1.26500,41.47000,121.20000],

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
        max_relative=5e-4
        );

    // now entropy 
    let s_test = s_ph_eqm(p, h);
    approx::assert_relative_eq!(
        s_kj_per_kg_k,
        s_test.get::<kilojoule_per_kilogram_kelvin>(),
        max_relative=2e-3
        );

    // cp 
    let cp_test = cp_ph_eqm(p, h);
    approx::assert_relative_eq!(
        cp_kj_per_kg_k,
        cp_test.get::<kilojoule_per_kilogram_kelvin>(),
        max_relative=1e-3
        );
    // w 
    let w_test = w_ph_eqm(p, h);
    approx::assert_relative_eq!(
        w_m_per_s,
        w_test.get::<meter_per_second>(),
        max_relative=1e-3
        );

    // kappa
    let kappa_test = kappa_ph_eqm(p, h);
    approx::assert_relative_eq!(
        kappa_dimensionless,
        kappa_test.get::<ratio>(),
        max_relative=6e-3
        );

    // eta and lambda tbd



}



