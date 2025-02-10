use uom::si::{available_energy::kilojoule_per_kilogram, dynamic_viscosity::micropascal_second};
use uom::si::pressure::bar;
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::specific_heat_capacity::kilojoule_per_kilogram_kelvin;
use uom::si::specific_volume::cubic_meter_per_kilogram;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::velocity::meter_per_second;

use crate::{dynamic_viscosity::mu_ph_eqm, interfaces::functional_programming::ph_flash_eqm::{cp_ph_eqm, kappa_ph_eqm, s_ph_eqm, t_ph_eqm, v_ph_eqm, w_ph_eqm}};

/// single phase table (see page 201)
///
/// NOTE: ph flash UNABLE to do triple point liquid and vapour accurately.
/// or even close to 0 degc
#[test]
pub fn single_phase_table_2_to_750_degc(){

    let steam_table: Vec<[f64; 10]> =
        vec![
        //[200.000,0.000,0.000990355,20.0338,0.0004733,4.129,1434.50000,103.89000,1751.20000,570.14000],
        [200.000,2.000,0.000990391,28.2906,0.030591,4.1278,1444.30000,105.32000,1639.70000,574.67000],
        [200.000,4.000,0.000990484,36.5454,0.060484,4.127,1453.70000,106.68000,1539.00000,579.04000],
        [200.000,6.000,0.00099063,44.799,0.090157,4.1265,1462.80000,107.99000,1448.00000,583.26000],
        [200.000,8.000,0.000990826,53.0518,0.11962,4.1263,1471.40000,109.25000,1365.20000,587.34000],
        [200.000,10.000,0.00099107,61.3043,0.14886,4.1262,1479.60000,110.45000,1289.80000,591.28000],
        [200.000,12.000,0.000991358,69.5568,0.17791,4.1263,1487.40000,111.59000,1220.80000,595.11000],
        [200.000,14.000,0.000991689,77.8096,0.20675,4.1265,1494.90000,112.68000,1157.60000,598.82000],
        [200.000,16.000,0.000992062,86.0628,0.23539,4.1267,1502.10000,113.71000,1099.50000,602.42000],
        [200.000,18.000,0.000992473,94.3166,0.26384,4.1271,1508.90000,114.70000,1045.90000,605.93000],
        [200.000,20.000,0.000992922,102.571,0.29209,4.1274,1515.30000,115.63000,996.40000,609.33000],
        [200.000,25.000,0.000994202,123.211,0.3619,4.1285,1530.10000,117.74000,887.91000,617.44000],
        [200.000,30.000,0.00099569,143.856,0.43057,4.1297,1543.00000,119.55000,797.30000,625.03000],
        [200.000,35.000,0.000997374,164.508,0.49814,4.131,1554.10000,121.09000,720.75000,632.13000],
        [200.000,40.000,0.00099924,185.166,0.56464,4.1325,1563.70000,122.35000,655.47000,638.77000],
        [200.000,45.000,0.00100128,205.833,0.63012,4.1342,1571.80000,123.36000,599.32000,644.99000],
        [200.000,50.000,0.00100348,226.509,0.6946,4.1361,1578.40000,124.14000,550.64000,650.78000],
        [200.000,55.000,0.00100585,247.195,0.75813,4.1384,1583.80000,124.69000,508.17000,656.18000],
        [200.000,60.000,0.00100837,267.893,0.82072,4.1409,1587.90000,125.03000,470.87000,661.19000],
        [200.000,65.000,0.00101103,288.604,0.88243,4.1437,1590.90000,125.17000,437.94000,665.82000],
        [200.000,70.000,0.00101385,309.33,0.94327,4.1468,1592.90000,125.13000,408.73000,670.08000],
        [200.000,75.000,0.0010168,330.073,1.0033,4.1503,1593.80000,124.91000,382.69000,673.98000],
        [200.000,80.000,0.0010199,350.834,1.0625,4.1542,1593.70000,124.52000,359.38000,677.53000],
        [200.000,85.000,0.00102313,371.615,1.1209,4.1584,1592.80000,123.98000,338.43000,680.73000],
        [200.000,90.000,0.00102651,392.419,1.1786,4.163,1591.00000,123.29000,319.53000,683.59000],
        [200.000,95.000,0.00103002,413.246,1.2356,4.168,1588.40000,122.47000,302.43000,686.13000],
        [200.000,100.000,0.00103366,434.1,1.2918,4.1734,1585.00000,121.52000,286.91000,688.34000],
        [200.000,110.000,0.00104137,475.892,1.4024,4.1853,1576.20000,119.28000,259.87000,691.81000],
        [200.000,120.000,0.00104964,517.811,1.5104,4.1989,1564.70000,116.63000,237.20000,694.07000],
        [200.000,130.000,0.00105847,559.877,1.616,4.2144,1550.80000,113.61000,218.02000,695.16000],
        [200.000,140.000,0.0010679,602.107,1.7195,4.232,1534.70000,110.28000,201.63000,695.14000],
        [200.000,150.000,0.00107795,644.524,1.8209,4.2518,1516.50000,106.67000,187.50000,694.06000],
        [200.000,160.000,0.00108865,687.152,1.9205,4.2742,1496.20000,102.82000,175.24000,691.94000],
        [200.000,170.000,0.00110005,730.018,2.0183,4.2996,1474.00000,98.76000,164.51000,689.03000],
        [200.000,180.000,0.00111219,773.155,2.1146,4.3284,1450.00000,94.51600,155.05000,685.34000],
        [200.000,190.000,0.00112515,816.599,2.2094,4.361,1424.00000,90.11200,146.65000,680.63000],
        [200.000,200.000,0.00113899,860.391,2.303,4.3982,1396.20000,85.57100,139.14000,675.00000],
        [200.000,210.000,0.0011538,904.58,2.3954,4.4405,1366.40000,80.91400,132.37000,668.48000],
        [200.000,220.000,0.0011697,949.222,2.4868,4.4889,1334.80000,76.16200,126.23000,661.11000],
        [200.000,230.000,0.00118681,994.382,2.5775,4.5444,1301.20000,71.33600,120.62000,652.89000],
        [200.000,240.000,0.0012053,1040.14,2.6675,4.6084,1265.70000,66.45700,115.45000,643.83000],
        [200.000,250.000,0.00122536,1086.58,2.7572,4.6824,1228.10000,61.54600,110.65000,633.92000],
        [200.000,260.000,0.00124724,1133.83,2.8466,4.7687,1188.50000,56.62200,106.16000,623.15000],
        [200.000,270.000,0.00127126,1182.01,2.9362,4.8701,1146.50000,51.70000,101.91000,611.51000],
        [200.000,280.000,0.00129784,1231.29,3.0261,4.9908,1102.10000,46.79200,97.84900,598.95000],
        [200.000,290.000,0.00132752,1281.91,3.1167,5.1366,1054.80000,41.90700,93.93100,585.42000],
        [200.000,300.000,0.00136109,1334.14,3.2087,5.3168,1004.30000,37.05000,90.09800,570.83000],
        [200.000,310.000,0.00139965,1388.4,3.3025,5.5459,949.79000,32.22600,86.29300,555.06000],
        [200.000,320.000,0.00144494,1445.3,3.3993,5.8491,890.30000,27.42800,82.44300,537.89000],
        [200.000,330.000,0.00149977,1505.79,3.5004,6.274,824.41000,22.65900,78.45000,519.00000],
        [200.000,340.000,0.00156931,1571.52,3.6085,6.9239,751.10000,17.97500,74.15900,497.81000],
        [200.000,350.000,0.00166487,1645.95,3.7288,8.1062,664.96000,13.28000,69.26600,473.33000],
        [200.000,360.000,0.00182472,1740.13,3.8787,11.46,542.74000,8.07160,62.79000,443.88000],
        [200.000,365.746,0.00203865,1827.1,4.0154,23.2,422.20000,4.37190,56.19800,432.42000],
        [200.000,365.746,0.00585828,2411.39,4.9299,45.677,384.50000,1.26180,27.40000,250.80000],
        [200.000,370.000,0.00692374,2526.48,5.1095,18.66,421.11000,1.28060,26.28500,168.42000],
        [200.000,380.000,0.00825779,2659.19,5.3144,10.221,461.33000,1.28860,25.82400,126.74000],
        [200.000,390.000,0.00918976,2747.17,5.4482,7.6925,487.08000,1.29080,25.89900,111.72000],
        [200.000,400.000,0.00994958,2816.84,5.5525,6.3601,507.34000,1.29350,26.13700,103.69000],
        [200.000,410.000,0.0106082,2876.05,5.6398,5.541,524.18000,1.29510,26.44900,98.77600],
        [200.000,420.000,0.0111994,2928.51,5.716,4.982,538.73000,1.29570,26.79900,95.54200],
        [200.000,430.000,0.0117416,2976.18,5.7843,4.5718,551.66000,1.29600,27.17200,93.34100],
        [200.000,440.000,0.0122459,3020.26,5.8466,4.2568,563.37000,1.29590,27.55800,91.83500],
        [200.000,450.000,0.0127202,3061.53,5.9041,4.0074,574.13000,1.29570,27.95300,90.82700],
        [200.000,460.000,0.0131699,3100.57,5.9577,3.8056,584.12000,1.29540,28.35300,90.19200],
        [200.000,470.000,0.0135992,3137.77,6.0081,3.6395,593.47000,1.29500,28.75500,89.84700],
        [200.000,480.000,0.0140113,3173.45,6.0558,3.501,602.28000,1.29450,29.16000,89.73300],
        [200.000,490.000,0.0144085,3207.86,6.1012,3.3841,610.63000,1.29390,29.56500,89.80700],
        [200.000,500.000,0.0147929,3241.19,6.1445,3.2845,618.58000,1.29330,29.97100,90.03900],
        [200.000,510.000,0.0151662,3273.59,6.1862,3.199,626.18000,1.29270,30.37600,90.40400],
        [200.000,520.000,0.0155296,3305.21,6.2263,3.1251,633.46000,1.29200,30.78000,90.88200],
        [200.000,530.000,0.0158842,3336.13,6.265,3.0608,640.47000,1.29120,31.18400,91.45900],
        [200.000,540.000,0.0162309,3366.45,6.3026,3.0046,647.23000,1.29040,31.58600,92.12200],
        [200.000,550.000,0.0165707,3396.24,6.339,2.9552,653.76000,1.28960,31.98700,92.86100],
        [200.000,560.000,0.016904,3425.57,6.3744,2.9116,660.09000,1.28880,32.38700,93.66800],
        [200.000,570.000,0.0172316,3454.49,6.4089,2.8731,666.23000,1.28790,32.78600,94.53400],
        [200.000,580.000,0.0175539,3483.05,6.4426,2.8389,672.20000,1.28700,33.18300,95.45600],
        [200.000,590.000,0.0178714,3511.28,6.4755,2.8084,678.02000,1.28620,33.58000,96.42500],
        [200.000,600.000,0.0181844,3539.23,6.5077,2.7812,683.69000,1.28530,33.97400,97.44000],
        [200.000,650.000,0.0196942,3675.59,6.6596,2.6824,710.26000,1.28080,35.92900,103.05000],
        [200.000,700.000,0.0211327,3808.15,6.7994,2.6251,734.48000,1.27640,37.85100,109.30000],
        [200.000,750.000,0.0225198,3938.52,6.9301,2.593,756.93000,1.27210,39.74500,116.07000],
        //[200.000,800.000,0.0238685,4067.73,7.0534,2.5775,777.99000,1.26790,41.61000,123.14000],

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
        max_relative=1e-3
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

    // lambda tbd
    //
    let eta_micropascal_second_test = mu_ph_eqm(p, h)
        .get::<micropascal_second>();
    approx::assert_relative_eq!(
        eta_micropascal_second,
        eta_micropascal_second_test,
        max_relative=2e-2
        );



}



