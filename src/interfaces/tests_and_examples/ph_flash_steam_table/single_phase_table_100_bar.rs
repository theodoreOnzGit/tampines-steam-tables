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
        //[100.000,0.000,0.0009952,10.0693,0.0003384,4.1723,1418.20000,202.09000,1770.60000,563.03000],
        [100.000,2.000,0.000995171,18.41,0.030762,4.1686,1428.00000,204.90000,1655.80000,567.79000],
        [100.000,4.000,0.000995203,26.744,0.060942,4.1655,1437.40000,207.61000,1552.50000,572.37000],
        [100.000,6.000,0.000995294,35.0726,0.090884,4.1631,1446.40000,210.21000,1459.10000,576.77000],
        [100.000,8.000,0.00099544,43.3967,0.1206,4.1611,1455.10000,212.69000,1374.50000,581.02000],
        [100.000,10.000,0.000995638,51.7173,0.15009,4.1595,1463.30000,215.07000,1297.40000,585.12000],
        [100.000,12.000,0.000995884,60.0349,0.17936,4.1582,1471.20000,217.34000,1227.10000,589.09000],
        [100.000,14.000,0.000996178,68.3503,0.20842,4.1572,1478.70000,219.50000,1162.60000,592.92000],
        [100.000,16.000,0.000996516,76.6637,0.23727,4.1563,1485.90000,221.56000,1103.50000,596.64000],
        [100.000,18.000,0.000996897,84.9757,0.26592,4.1557,1492.70000,223.51000,1049.00000,600.24000],
        [100.000,20.000,0.000997318,93.2865,0.29437,4.1551,1499.20000,225.36000,998.78000,603.74000],
        [100.000,25.000,0.000998541,114.06,0.36463,4.1543,1514.00000,229.54000,888.81000,612.04000],
        [100.000,30.000,0.000999989,134.831,0.43372,4.1541,1526.80000,233.12000,797.14000,619.76000],
        [100.000,35.000,0.00100165,155.601,0.50168,4.1543,1537.90000,236.13000,719.85000,626.96000],
        [100.000,40.000,0.0010035,176.374,0.56855,4.1549,1547.40000,238.60000,654.03000,633.67000],
        [100.000,45.000,0.00100554,197.151,0.63437,4.1559,1555.30000,240.56000,597.49000,639.92000],
        [100.000,50.000,0.00100775,217.934,0.69919,4.1573,1561.80000,242.04000,548.54000,645.74000],
        [100.000,55.000,0.00101013,238.725,0.76303,4.1591,1566.90000,243.06000,505.86000,651.13000],
        [100.000,60.000,0.00101268,259.526,0.82594,4.1613,1570.80000,243.66000,468.43000,656.13000],
        [100.000,65.000,0.00101538,280.339,0.88795,4.1639,1573.60000,243.86000,435.41000,660.73000],
        [100.000,70.000,0.00101824,301.166,0.94909,4.167,1575.20000,243.67000,406.13000,664.95000],
        [100.000,75.000,0.00102125,322.009,1.0094,4.1705,1575.80000,243.14000,380.04000,668.80000],
        [100.000,80.000,0.00102441,342.871,1.0689,4.1744,1575.40000,242.27000,356.71000,672.29000],
        [100.000,85.000,0.00102771,363.754,1.1276,4.1787,1574.10000,241.09000,335.75000,675.42000],
        [100.000,90.000,0.00103116,384.659,1.1856,4.1835,1571.90000,239.62000,316.85000,678.21000],
        [100.000,95.000,0.00103475,405.59,1.2428,4.1888,1568.90000,237.88000,299.76000,680.67000],
        [100.000,100.000,0.0010385,426.548,1.2994,4.1945,1565.10000,235.89000,284.25000,682.79000],
        [100.000,110.000,0.00104641,468.555,1.4105,4.2073,1555.40000,231.20000,257.24000,686.08000],
        [100.000,120.000,0.00105493,510.701,1.519,4.2221,1543.00000,225.68000,234.61000,688.14000],
        [100.000,130.000,0.00106405,553.005,1.6253,4.2391,1528.00000,219.43000,215.46000,689.02000],
        [100.000,140.000,0.0010738,595.491,1.7294,4.2585,1510.70000,212.54000,199.11000,688.77000],
        [100.000,150.000,0.00108422,638.184,1.8315,4.2806,1491.20000,205.10000,185.02000,687.44000],
        [100.000,160.000,0.00109535,681.112,1.9318,4.3057,1469.60000,197.16000,172.79000,685.06000],
        [100.000,170.000,0.00110724,724.309,2.0304,4.3343,1445.80000,188.80000,162.07000,682.10000],
        [100.000,180.000,0.00111996,767.812,2.1274,4.367,1420.10000,180.06000,152.62000,678.03000],
        [100.000,190.000,0.00113357,811.665,2.2232,4.4044,1392.20000,170.99000,144.23000,672.97000],
        [100.000,200.000,0.00114818,855.918,2.3177,4.4472,1362.30000,161.64000,136.71000,666.98000],
        [100.000,210.000,0.00116389,900.631,2.4112,4.4965,1330.30000,152.05000,129.92000,660.08000],
        [100.000,220.000,0.00118085,945.874,2.5039,4.5535,1296.10000,142.27000,123.76000,652.28000],
        [100.000,230.000,0.00119922,991.731,2.5959,4.6196,1259.80000,132.34000,118.10000,643.58000],
        [100.000,240.000,0.00121923,1038.3,2.6876,4.697,1221.10000,122.30000,112.87000,633.98000],
        [100.000,250.000,0.00124116,1085.72,2.7791,4.7883,1180.00000,112.19000,107.99000,623.46000],
        [100.000,260.000,0.00126534,1134.13,2.8708,4.8972,1136.30000,102.03000,103.38000,611.99000],
        [100.000,270.000,0.00129228,1183.74,2.9629,5.0293,1089.40000,91.84000,98.99500,599.52000],
        [100.000,280.000,0.00132264,1234.82,3.0561,5.1931,1038.90000,81.60500,94.75500,585.97000],
        [100.000,290.000,0.00135739,1287.75,3.151,5.4023,983.78000,71.30000,90.59600,571.21000],
        [100.000,300.000,0.00139804,1343.1,3.2484,5.6816,922.76000,60.90500,86.43400,555.07000],
        [100.000,310.000,0.0014471,1401.77,3.3498,6.0782,854.92000,50.50700,82.15600,537.19000],
        [100.000,310.999,0.00145262,1407.87,3.3603,6.1275,847.74000,49.47400,81.71600,535.29000],
        [100.000,310.999,0.0180336,2725.47,5.6159,7.1472,472.44000,1.23770,20.19400,78.33800],
        [100.000,320.000,0.0192716,2782.66,5.7131,5.7468,491.71000,1.25460,20.66400,74.14800],
        [100.000,330.000,0.0204462,2835.67,5.8017,4.9228,508.20000,1.26320,21.17700,71.58000],
        [100.000,340.000,0.0214897,2882.06,5.878,4.3885,522.16000,1.26880,21.68200,70.04000],
        [100.000,350.000,0.0224422,2923.96,5.9458,4.0118,534.45000,1.27280,22.17700,69.10400],
        [100.000,360.000,0.0233274,2962.61,6.0073,3.7324,545.52000,1.27570,22.66500,68.57000],
        [100.000,370.000,0.0241605,2998.82,6.0641,3.5174,555.64000,1.27790,23.14600,68.32300],
        [100.000,380.000,0.0249522,3033.11,6.117,3.3471,565.02000,1.27940,23.62100,68.29300],
        [100.000,390.000,0.0257099,3065.87,6.1668,3.2092,573.79000,1.28060,24.08900,68.43400],
        [100.000,400.000,0.0264393,3097.38,6.2139,3.0958,582.04000,1.28130,24.55300,68.71500],
        [100.000,410.000,0.0271447,3127.85,6.2589,3.0013,589.86000,1.28180,25.01100,69.11100],
        [100.000,420.000,0.0278294,3157.45,6.3019,2.9217,597.31000,1.28200,25.46500,69.60600],
        [100.000,430.000,0.0284963,3186.32,6.3432,2.8542,604.44000,1.28210,25.91400,70.18600],
        [100.000,440.000,0.0291475,3214.57,6.3831,2.7965,611.28000,1.28200,26.36000,70.84000],
        [100.000,450.000,0.029785,3242.28,6.4217,2.747,617.87000,1.28170,26.80200,71.56000],
        [100.000,460.000,0.0304102,3269.53,6.4591,2.7043,624.24000,1.28140,27.24100,72.33800],
        [100.000,470.000,0.0310246,3296.38,6.4955,2.6674,630.41000,1.28100,27.67600,73.16700],
        [100.000,480.000,0.0316292,3322.89,6.531,2.6354,636.39000,1.28040,28.10900,74.04400],
        [100.000,490.000,0.032225,3349.11,6.5655,2.6076,642.21000,1.27990,28.53900,74.96400],
        [100.000,500.000,0.0328129,3375.06,6.5993,2.5833,647.89000,1.27920,28.96600,75.92200],
        [100.000,510.000,0.0333935,3400.78,6.6324,2.5622,653.42000,1.27860,29.39100,76.91600],
        [100.000,520.000,0.0339675,3426.31,6.6648,2.5437,658.83000,1.27790,29.81300,77.94300],
        [100.000,530.000,0.0345355,3451.67,6.6965,2.5275,664.12000,1.27710,30.23300,79.00000],
        [100.000,540.000,0.0350979,3476.87,6.7277,2.5134,669.31000,1.27640,30.65100,80.08500],
        [100.000,550.000,0.0356552,3501.94,6.7584,2.5011,674.39000,1.27560,31.06700,81.19500],
        [100.000,560.000,0.0362078,3526.9,6.7885,2.4904,679.39000,1.27480,31.48000,82.33000],
        [100.000,570.000,0.0367561,3551.75,6.8182,2.4811,684.29000,1.27400,31.89200,83.48700],
        [100.000,580.000,0.0373002,3576.52,6.8474,2.473,689.11000,1.27310,32.30200,84.66500],
        [100.000,590.000,0.0378406,3601.22,6.8761,2.466,693.86000,1.27230,32.71100,85.86200],
        [100.000,600.000,0.0383775,3625.84,6.9045,2.46,698.54000,1.27150,33.11700,87.07700],
        [100.000,650.000,0.0410163,3748.32,7.0409,2.442,720.95000,1.26720,35.12600,93.38700],
        [100.000,700.000,0.0435944,3870.27,7.1696,2.438,742.03000,1.26300,37.09800,100.00000],
        [100.000,750.000,0.0461269,3992.28,7.2918,2.4435,762.03000,1.25890,39.03800,106.87000],
        //[100.000,800.000,0.0486242,4114.73,7.4087,2.4555,781.12000,1.25480,40.94700,113.90000],

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
        max_relative=3e-4
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



