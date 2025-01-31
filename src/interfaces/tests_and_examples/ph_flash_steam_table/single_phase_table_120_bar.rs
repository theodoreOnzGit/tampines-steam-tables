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
        //[120.000,0.000,0.000994218,12.0737,0.0003931,4.1633,1421.40000,169.34000,1766.60000,564.48000],
        [120.000,2.000,0.000994202,20.397,0.030754,4.1601,1431.20000,171.70000,1652.50000,569.19000],
        [120.000,4.000,0.000994247,28.7145,0.060873,4.1575,1440.60000,173.96000,1549.70000,573.72000],
        [120.000,6.000,0.00099435,37.0275,0.09076,4.1555,1449.70000,176.12000,1456.80000,578.09000],
        [120.000,8.000,0.000994506,45.3369,0.12042,4.1539,1458.30000,178.20000,1372.50000,582.30000],
        [120.000,10.000,0.000994713,53.6433,0.14986,4.1526,1466.60000,180.19000,1295.80000,586.37000],
        [120.000,12.000,0.000994968,61.9475,0.17909,4.1516,1474.40000,182.08000,1225.80000,590.31000],
        [120.000,14.000,0.000995269,70.2499,0.2081,4.1508,1481.90000,183.88000,1161.60000,594.12000],
        [120.000,16.000,0.000995614,78.5509,0.23691,4.1502,1489.10000,185.60000,1102.60000,597.81000],
        [120.000,18.000,0.000996001,86.8509,0.26551,4.1498,1495.90000,187.23000,1048.40000,601.40000],
        [120.000,20.000,0.000996429,95.15,0.29392,4.1494,1502.40000,188.77000,998.27000,604.87000],
        [120.000,25.000,0.000997663,115.896,0.36409,4.149,1517.20000,192.26000,888.60000,613.13000],
        [120.000,30.000,0.00099912,136.641,0.43309,4.1491,1530.00000,195.26000,797.15000,620.83000],
        [120.000,35.000,0.00100078,157.387,0.50097,4.1495,1541.20000,197.78000,720.02000,628.01000],
        [120.000,40.000,0.00100264,178.136,0.56777,4.1503,1550.60000,199.84000,654.31000,634.70000],
        [120.000,45.000,0.00100468,198.891,0.63352,4.1515,1558.60000,201.49000,597.85000,640.95000],
        [120.000,50.000,0.00100689,219.652,0.69827,4.153,1565.10000,202.73000,548.95000,646.76000],
        [120.000,55.000,0.00100927,240.421,0.76205,4.1549,1570.30000,203.60000,506.32000,652.15000],
        [120.000,60.000,0.00101181,261.201,0.82489,4.1571,1574.20000,204.11000,468.92000,657.15000],
        [120.000,65.000,0.0010145,281.993,0.88684,4.1598,1577.00000,204.29000,435.91000,661.75000],
        [120.000,70.000,0.00101735,302.8,0.94792,4.1629,1578.70000,204.15000,406.65000,665.98000],
        [120.000,75.000,0.00102035,323.622,1.0082,4.1664,1579.40000,203.72000,380.57000,669.84000],
        [120.000,80.000,0.00102349,344.464,1.0676,4.1703,1579.10000,203.02000,357.24000,673.34000],
        [120.000,85.000,0.00102678,365.326,1.1263,4.1746,1577.80000,202.05000,336.28000,676.49000],
        [120.000,90.000,0.00103022,386.21,1.1842,4.1793,1575.70000,200.85000,317.39000,679.30000],
        [120.000,95.000,0.00103379,407.12,1.2414,4.1845,1572.80000,199.41000,300.30000,681.77000],
        [120.000,100.000,0.00103751,428.056,1.2978,4.1902,1569.10000,197.77000,284.78000,683.91000],
        [120.000,110.000,0.00104539,470.02,1.4088,4.2028,1559.60000,193.89000,257.76000,687.24000],
        [120.000,120.000,0.00105385,512.119,1.5173,4.2174,1547.40000,189.33000,235.13000,689.34000],
        [120.000,130.000,0.00106291,554.374,1.6234,4.234,1532.60000,184.16000,215.98000,690.26000],
        [120.000,140.000,0.0010726,596.807,1.7274,4.253,1515.60000,178.46000,199.62000,690.06000],
        [120.000,150.000,0.00108294,639.443,1.8294,4.2746,1496.40000,172.30000,185.52000,688.77000],
        [120.000,160.000,0.00109398,682.309,1.9295,4.2992,1475.00000,165.73000,173.28000,686.44000],
        [120.000,170.000,0.00110577,725.438,2.0279,4.3271,1451.60000,158.80000,162.57000,683.51000],
        [120.000,180.000,0.00111836,768.865,2.1248,4.359,1426.20000,151.56000,153.11000,679.50000],
        [120.000,190.000,0.00113184,812.632,2.2204,4.3953,1398.80000,144.05000,144.72000,674.52000],
        [120.000,200.000,0.00114628,856.788,2.3147,4.4369,1369.30000,136.31000,137.20000,668.60000],
        [120.000,210.000,0.0011618,901.391,2.408,4.4847,1337.80000,128.36000,130.42000,661.78000],
        [120.000,220.000,0.00117853,946.506,2.5004,4.5397,1304.10000,120.26000,124.26000,654.07000],
        [120.000,230.000,0.00119664,992.215,2.5921,4.6035,1268.40000,112.04000,118.61000,645.47000],
        [120.000,240.000,0.00121632,1038.61,2.6834,4.6779,1230.40000,103.72000,113.40000,635.98000],
        [120.000,250.000,0.00123783,1085.81,2.7745,4.7652,1190.10000,95.34600,108.53000,625.59000],
        [120.000,260.000,0.00126151,1133.97,2.8657,4.8688,1147.20000,86.93700,103.96000,614.27000],
        [120.000,270.000,0.00128779,1183.26,2.9573,4.9936,1101.50000,78.50800,99.60000,601.98000],
        [120.000,280.000,0.00131728,1233.94,3.0498,5.1467,1052.40000,70.06300,95.40300,588.65000],
        [120.000,290.000,0.00135084,1286.33,3.1436,5.3398,999.16000,61.58700,91.30200,574.18000],
        [120.000,300.000,0.00138976,1340.93,3.2397,5.5924,940.65000,53.05600,87.22200,558.39000],
        [120.000,310.000,0.00143614,1398.49,3.3393,5.9411,875.62000,44.48900,83.06700,541.03000],
        [120.000,320.000,0.00149369,1460.31,3.4444,6.4621,803.26000,35.99800,78.69500,521.62000],
        [120.000,324.678,0.00152633,1491.33,3.4965,6.8126,765.59000,32.00100,76.51100,511.64000],
        [120.000,324.678,0.0142689,2685.58,5.4941,8.8189,459.46000,1.23290,21.10900,91.06200],
        [120.000,330.000,0.0150236,2728.14,5.565,7.3313,474.22000,1.24740,21.35700,85.82700],
        [120.000,340.000,0.0162112,2793.47,5.6725,5.8968,494.92000,1.25910,21.83700,80.37400],
        [120.000,350.000,0.0172227,2848.01,5.7607,5.0746,511.54000,1.26610,22.32100,77.30400],
        [120.000,360.000,0.0181226,2895.87,5.8369,4.5309,525.77000,1.27110,22.80400,75.41300],
        [120.000,370.000,0.0189442,2939.15,5.9047,4.1447,538.33000,1.27480,23.28300,74.22200],
        [120.000,380.000,0.0197077,2979.09,5.9664,3.8563,549.63000,1.27740,23.75600,73.49400],
        [120.000,390.000,0.0204258,3016.49,6.0232,3.6329,559.97000,1.27930,24.22500,73.09600],
        [120.000,400.000,0.0211077,3051.9,6.0762,3.4551,569.54000,1.28060,24.68900,72.94800],
        [120.000,410.000,0.0217597,3085.7,6.1261,3.3106,578.47000,1.28150,25.14900,72.99500],
        [120.000,420.000,0.0223867,3118.19,6.1733,3.1912,586.89000,1.28210,25.60400,73.20000],
        [120.000,430.000,0.0229924,3149.59,6.2182,3.0914,594.85000,1.28250,26.05400,73.53600],
        [120.000,440.000,0.0235799,3180.07,6.2613,3.0071,602.43000,1.28260,26.50100,73.98200],
        [120.000,450.000,0.0241515,3209.77,6.3027,2.9353,609.68000,1.28260,26.94400,74.52200],
        [120.000,460.000,0.0247091,3238.81,6.3425,2.8737,616.64000,1.28240,27.38400,75.14500],
        [120.000,470.000,0.0252545,3267.28,6.3811,2.8207,623.34000,1.28210,27.82000,75.84000],
        [120.000,480.000,0.025789,3295.25,6.4185,2.7748,629.81000,1.28170,28.25300,76.59800],
        [120.000,490.000,0.0263138,3322.79,6.4548,2.7349,636.07000,1.28130,28.68300,77.41400],
        [120.000,500.000,0.0268298,3349.97,6.4902,2.7002,642.14000,1.28070,29.11100,78.28100],
        [120.000,510.000,0.0273378,3376.81,6.5247,2.6698,648.04000,1.28010,29.53600,79.19300],
        [120.000,520.000,0.0278387,3403.37,6.5584,2.6432,653.79000,1.27950,29.95800,80.14800],
        [120.000,530.000,0.0283331,3429.69,6.5914,2.6199,659.39000,1.27880,30.37800,81.14200],
        [120.000,540.000,0.0288215,3455.78,6.6237,2.5993,664.86000,1.27810,30.79600,82.17000],
        [120.000,550.000,0.0293045,3481.68,6.6553,2.5813,670.22000,1.27740,31.21100,83.23100],
        [120.000,560.000,0.0297824,3507.41,6.6864,2.5654,675.46000,1.27660,31.62500,84.32100],
        [120.000,570.000,0.0302557,3533.0,6.7169,2.5514,680.60000,1.27580,32.03600,85.43800],
        [120.000,580.000,0.0307247,3558.45,6.7469,2.5391,685.64000,1.27510,32.44600,86.58100],
        [120.000,590.000,0.0311897,3583.78,6.7764,2.5283,690.60000,1.27430,32.85300,87.74800],
        [120.000,600.000,0.0316511,3609.02,6.8055,2.5188,695.46000,1.27340,33.25900,88.93600],
        [120.000,650.000,0.0339104,3734.07,6.9448,2.4871,718.70000,1.26930,35.26400,95.15800],
        [120.000,700.000,0.0361069,3858.03,7.0756,2.4737,740.40000,1.26520,37.23100,101.74000],
        [120.000,750.000,0.0382567,3981.64,7.1994,2.4723,760.90000,1.26120,39.16500,108.60000],
        //[120.000,800.000,0.0403706,4105.4,7.3175,2.4793,780.40000,1.25710,41.06800,115.66000],

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
        max_relative=3e-4
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
        max_relative=1e-3
        );
    // w 
    let w_test = w_ph_eqm(p, h);
    approx::assert_relative_eq!(
        w_m_per_s,
        w_test.get::<meter_per_second>(),
        max_relative=3e-4
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



