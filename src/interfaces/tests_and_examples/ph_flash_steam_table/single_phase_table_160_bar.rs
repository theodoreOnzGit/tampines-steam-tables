use uom::si::{available_energy::kilojoule_per_kilogram, dynamic_viscosity::micropascal_second};
use uom::si::pressure::bar;
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::velocity::meter_per_second;

use crate::dynamic_viscosity::mu_ph_eqm;
use crate::interfaces::functional_programming::ph_flash_eqm::{cp_ph_eqm, kappa_ph_eqm, s_ph_eqm, t_ph_eqm, v_ph_eqm, w_ph_eqm};

/// single phase table (see page 201)
///
/// NOTE: ph flash UNABLE to do triple point liquid and vapour accurately.
/// or even close to 0 degc
#[test]
pub fn single_phase_table_2_to_750_degc(){

    let steam_table: Vec<[f64; 10]> =
        vec![
        //[160.000,0.000,0.000992274,16.0651,0.0004606,4.1458,1427.90000,128.43000,1758.70000,567.34000],
        [160.000,2.000,0.000992285,24.3545,0.030697,4.1437,1437.70000,130.20000,1645.90000,571.96000],
        [160.000,4.000,0.000992354,32.6401,0.060701,4.142,1447.20000,131.90000,1544.30000,576.41000],
        [160.000,6.000,0.000992478,40.9228,0.090479,4.1408,1456.20000,133.53000,1452.30000,580.70000],
        [160.000,8.000,0.000992655,49.2033,0.12004,4.1399,1464.80000,135.10000,1368.80000,584.84000],
        [160.000,10.000,0.00099288,57.4824,0.14938,4.1392,1473.10000,136.59000,1292.70000,588.85000],
        [160.000,12.000,0.000993152,65.7603,0.17851,4.1387,1480.90000,138.01000,1223.20000,592.73000],
        [160.000,14.000,0.000993469,74.0374,0.20744,4.1385,1488.40000,139.37000,1159.50000,596.49000],
        [160.000,16.000,0.000993828,82.3141,0.23616,4.1383,1495.60000,140.66000,1101.00000,600.14000],
        [160.000,18.000,0.000994227,90.5906,0.26469,4.1382,1502.40000,141.89000,1047.10000,603.68000],
        [160.000,20.000,0.000994665,98.8671,0.29302,4.1382,1508.80000,143.05000,997.29000,607.12000],
        [160.000,25.000,0.000995923,119.559,0.36301,4.1386,1523.60000,145.68000,888.23000,615.30000],
        [160.000,30.000,0.000997395,140.253,0.43184,4.1392,1536.50000,147.94000,797.20000,622.94000],
        [160.000,35.000,0.000999068,160.952,0.49956,4.1401,1547.60000,149.84000,720.37000,630.08000],
        [160.000,40.000,0.00100093,181.655,0.56621,4.1413,1557.20000,151.41000,654.88000,636.75000],
        [160.000,45.000,0.00100297,202.365,0.63182,4.1427,1565.20000,152.66000,598.57000,642.98000],
        [160.000,50.000,0.00100518,223.083,0.69643,4.1445,1571.80000,153.61000,549.79000,648.78000],
        [160.000,55.000,0.00100755,243.81,0.76008,4.1465,1577.00000,154.28000,507.24000,654.17000],
        [160.000,60.000,0.00101008,264.549,0.8228,4.1489,1581.10000,154.68000,469.89000,659.17000],
        [160.000,65.000,0.00101276,285.3,0.88463,4.1516,1584.00000,154.84000,436.93000,663.79000],
        [160.000,70.000,0.00101559,306.066,0.94559,4.1548,1585.80000,154.76000,407.69000,668.04000],
        [160.000,75.000,0.00101856,326.848,1.0057,4.1583,1586.60000,154.46000,381.63000,671.92000],
        [160.000,80.000,0.00102168,347.649,1.065,4.1621,1586.40000,153.95000,358.31000,675.44000],
        [160.000,85.000,0.00102494,368.47,1.1236,4.1664,1585.30000,153.25000,337.36000,678.62000],
        [160.000,90.000,0.00102835,389.314,1.1814,4.1711,1583.40000,152.37000,318.46000,681.45000],
        [160.000,95.000,0.00103189,410.182,1.2384,4.1762,1580.60000,151.32000,301.37000,683.95000],
        [160.000,100.000,0.00103558,431.076,1.2948,4.1817,1577.10000,150.11000,285.85000,686.13000],
        [160.000,110.000,0.00104337,472.953,1.4056,4.1939,1567.90000,147.26000,258.82000,689.53000],
        [160.000,120.000,0.00105173,514.961,1.5138,4.208,1556.10000,143.89000,236.17000,691.71000],
        [160.000,130.000,0.00106067,557.12,1.6197,4.2241,1541.80000,140.07000,217.00000,692.72000],
        [160.000,140.000,0.00107023,599.45,1.7234,4.2423,1525.20000,135.85000,200.62000,692.61000],
        [160.000,150.000,0.00108042,641.975,1.8251,4.263,1506.50000,131.29000,186.52000,691.42000],
        [160.000,160.000,0.00109129,684.72,1.925,4.2865,1485.70000,126.42000,174.26000,689.20000],
        [160.000,170.000,0.00110288,727.715,2.0231,4.3131,1463.00000,121.29000,163.54000,686.29000],
        [160.000,180.000,0.00111524,770.994,2.1197,4.3434,1438.20000,115.92000,154.09000,682.44000],
        [160.000,190.000,0.00112845,814.596,2.2148,4.3778,1411.60000,110.36000,145.69000,677.59000],
        [160.000,200.000,0.00114258,858.566,2.3088,4.4171,1382.90000,104.62000,138.17000,671.81000],
        [160.000,210.000,0.00115774,902.956,2.4016,4.462,1352.30000,98.72800,131.40000,665.15000],
        [160.000,220.000,0.00117404,947.828,2.4935,4.5136,1319.80000,92.72200,125.25000,657.61000],
        [160.000,230.000,0.00119163,993.254,2.5847,4.573,1285.10000,86.62200,119.63000,649.21000],
        [160.000,240.000,0.00121069,1039.32,2.6754,4.6418,1248.40000,80.45600,114.44000,639.94000],
        [160.000,250.000,0.00123144,1086.13,2.7657,4.722,1209.50000,74.24800,109.61000,629.80000],
        [160.000,260.000,0.00125418,1133.81,2.856,4.8163,1168.30000,68.02000,105.07000,618.76000],
        [160.000,270.000,0.00127926,1182.51,2.9465,4.9284,1124.60000,61.78700,100.77000,606.81000],
        [160.000,280.000,0.0013072,1232.45,3.0376,5.0636,1078.00000,55.56100,96.65300,593.88000],
        [160.000,290.000,0.00133866,1283.89,3.1297,5.2301,1028.00000,49.34400,92.65200,579.91000],
        [160.000,300.000,0.00137464,1337.2,3.2236,5.4411,973.99000,43.13200,88.70900,564.77000],
        [160.000,310.000,0.00141664,1392.93,3.32,5.719,914.67000,36.91000,84.75200,548.28000],
        [160.000,320.000,0.00146711,1451.94,3.4203,6.1063,848.79000,30.69200,80.68300,530.13000],
        [160.000,330.000,0.00153049,1515.71,3.5269,6.6932,775.52000,24.56000,76.35100,509.81000],
        [160.000,340.000,0.00161627,1587.27,3.6445,7.7439,687.82000,18.29500,71.46900,486.41000],
        [160.000,347.357,0.00170954,1649.67,3.7457,9.4729,597.67000,13.05900,67.05600,466.33000],
        [160.000,347.357,0.00930813,2580.8,5.2463,15.207,429.20000,1.23690,23.36300,132.79000],
        [160.000,350.000,0.00976565,2616.99,5.3045,12.413,441.49000,1.24750,23.35500,121.52000],
        [160.000,360.000,0.0110599,2715.63,5.4616,8.1928,472.82000,1.26340,23.57000,102.25000],
        [160.000,370.000,0.0120464,2788.3,5.5755,6.5172,494.85000,1.27050,23.92200,93.94400],
        [160.000,380.000,0.0128781,2848.27,5.668,5.5538,512.76000,1.27600,24.32100,89.24300],
        [160.000,390.000,0.0136131,2900.49,5.7474,4.9289,527.98000,1.27980,24.74100,86.30100],
        [160.000,400.000,0.014281,2947.46,5.8177,4.4882,541.32000,1.28240,25.17100,84.37400],
        [160.000,410.000,0.0148991,2990.62,5.8814,4.1588,553.28000,1.28410,25.60600,83.10300],
        [160.000,420.000,0.0154783,3030.88,5.9399,3.9032,564.18000,1.28530,26.04200,82.29000],
        [160.000,430.000,0.0160263,3068.85,5.9943,3.6993,574.24000,1.28600,26.47800,81.81700],
        [160.000,440.000,0.0165486,3104.99,6.0453,3.5333,583.62000,1.28640,26.91300,81.60500],
        [160.000,450.000,0.0170494,3139.61,6.0935,3.396,592.43000,1.28660,27.34600,81.60000],
        [160.000,460.000,0.017532,3172.98,6.1393,3.2809,600.76000,1.28660,27.77700,81.76400],
        [160.000,470.000,0.0179988,3205.29,6.1831,3.1835,608.68000,1.28650,28.20600,82.07000],
        [160.000,480.000,0.018452,3236.7,6.2251,3.1004,616.23000,1.28620,28.63300,82.49400],
        [160.000,490.000,0.0188931,3267.34,6.2655,3.0289,623.47000,1.28590,29.05700,83.02100],
        [160.000,500.000,0.0193237,3297.31,6.3045,2.9671,630.42000,1.28540,29.47900,83.63700],
        [160.000,510.000,0.0197449,3326.71,6.3423,2.9134,637.13000,1.28490,29.89900,84.33000],
        [160.000,520.000,0.0201577,3355.6,6.379,2.8664,643.61000,1.28440,30.31700,85.09300],
        [160.000,530.000,0.0205628,3384.05,6.4146,2.8253,649.89000,1.28370,30.73300,85.91600],
        [160.000,540.000,0.0209611,3412.12,6.4494,2.7891,655.98000,1.28310,31.14600,86.79500],
        [160.000,550.000,0.0213532,3439.85,6.4832,2.7572,661.91000,1.28240,31.55800,87.72300],
        [160.000,560.000,0.0217396,3467.28,6.5164,2.729,667.68000,1.28160,31.96700,88.69600],
        [160.000,570.000,0.0221208,3494.44,6.5488,2.7039,673.31000,1.28090,32.37500,89.71000],
        [160.000,580.000,0.0224972,3521.37,6.5805,2.6817,678.81000,1.28010,32.78000,90.76100],
        [160.000,590.000,0.0228693,3548.08,6.6117,2.662,684.18000,1.27930,33.18400,91.84700],
        [160.000,600.000,0.0232373,3574.61,6.6422,2.6444,689.45000,1.27850,33.58700,92.96300],
        [160.000,650.000,0.0250259,3705.11,6.7876,2.5817,714.35000,1.27440,35.57300,98.93500],
        [160.000,700.000,0.0267475,3833.26,6.9228,2.5477,737.32000,1.27030,37.52300,105.39000],
        [160.000,750.000,0.02842,3960.18,7.0499,2.5316,758.80000,1.26620,39.44100,112.23000],
        //[160.000,800.000,0.0300554,4086.62,7.1706,2.5277,779.09000,1.26220,41.32800,119.31000],

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
        max_relative=2e-3
        );

    // lambda tbd
    //
    let eta_micropascal_second_test = mu_ph_eqm(p, h)
        .get::<micropascal_second>();
    approx::assert_relative_eq!(
        eta_micropascal_second,
        eta_micropascal_second_test,
        max_relative=3e-2
        );



}



