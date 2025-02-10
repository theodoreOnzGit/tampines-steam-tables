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
pub fn single_phase_table_2_to_750_degc(){

    let steam_table: Vec<[f64; 10]> =
        vec![
        //[500.000,0.000,0.00097673,49.1286,-0.001021,4.0225,1485.90000,45.21200,1704.70000,589.62000],
        [500.000,2.000,0.000976937,57.1781,0.028341,4.027,1495.60000,45.79100,1601.00000,593.57000],
        [500.000,4.000,0.000977186,65.2365,0.057522,4.0313,1504.80000,46.34500,1507.00000,597.41000],
        [500.000,6.000,0.000977475,73.3031,0.086523,4.0353,1513.60000,46.87500,1421.60000,601.15000],
        [500.000,8.000,0.000977801,81.3773,0.11534,4.039,1522.00000,47.38100,1343.60000,604.80000],
        [500.000,10.000,0.000978164,89.4588,0.14399,4.0424,1530.00000,47.86400,1272.20000,608.36000],
        [500.000,12.000,0.000978561,97.5468,0.17245,4.0456,1537.70000,48.32400,1206.70000,611.84000],
        [500.000,14.000,0.000978991,105.641,0.20074,4.0486,1545.00000,48.76300,1146.50000,615.25000],
        [500.000,16.000,0.000979453,113.741,0.22885,4.0513,1551.90000,49.18100,1090.90000,618.57000],
        [500.000,18.000,0.000979945,121.846,0.25678,4.0538,1558.60000,49.57800,1039.60000,621.83000],
        [500.000,20.000,0.000980468,129.956,0.28454,4.0562,1564.90000,49.95600,991.97000,625.01000],
        [500.000,25.000,0.000981897,150.25,0.35319,4.0614,1579.50000,50.81500,887.21000,632.68000],
        [500.000,30.000,0.000983494,170.569,0.42077,4.0659,1592.30000,51.55900,799.22000,639.94000],
        [500.000,35.000,0.000985249,190.908,0.48732,4.0698,1603.50000,52.19600,724.54000,646.82000],
        [500.000,40.000,0.000987154,211.266,0.55285,4.0733,1613.30000,52.73100,660.58000,653.33000],
        [500.000,45.000,0.000989203,231.64,0.6174,4.0765,1621.60000,53.16900,605.37000,659.47000],
        [500.000,50.000,0.000991389,252.031,0.68099,4.0797,1628.70000,53.51600,557.37000,665.26000],
        [500.000,55.000,0.000993708,272.437,0.74366,4.0828,1634.60000,53.77600,515.37000,670.70000],
        [500.000,60.000,0.000996157,292.859,0.80542,4.086,1639.30000,53.95500,478.40000,675.78000],
        [500.000,65.000,0.000998732,313.297,0.86631,4.0892,1643.00000,54.05600,445.69000,680.53000],
        [500.000,70.000,0.00100143,333.751,0.92636,4.0926,1645.60000,54.08400,416.61000,684.93000],
        [500.000,75.000,0.00100425,354.223,0.98558,4.0961,1647.30000,54.04300,390.65000,689.00000],
        [500.000,80.000,0.00100719,374.713,1.044,4.0998,1648.10000,53.93600,367.37000,692.73000],
        [500.000,85.000,0.00101025,395.221,1.1017,4.1037,1648.00000,53.76900,346.43000,696.14000],
        [500.000,90.000,0.00101343,415.75,1.1586,4.1078,1647.20000,53.54400,327.51000,699.23000],
        [500.000,95.000,0.00101672,436.3,1.2148,4.1121,1645.50000,53.26600,310.37000,702.00000],
        [500.000,100.000,0.00102013,456.872,1.2703,4.1167,1643.20000,52.93700,294.80000,704.46000],
        [500.000,110.000,0.0010273,498.087,1.3793,4.1265,1636.50000,52.14000,267.62000,708.47000],
        [500.000,120.000,0.00103494,539.405,1.4858,4.1374,1627.40000,51.17800,244.82000,711.31000],
        [500.000,130.000,0.00104306,580.838,1.5898,4.1494,1616.00000,50.07100,225.49000,713.02000],
        [500.000,140.000,0.00105167,622.397,1.6917,4.1627,1602.50000,48.84000,208.96000,713.65000],
        [500.000,150.000,0.00106078,664.096,1.7914,4.1774,1587.20000,47.49800,194.71000,713.25000],
        [500.000,160.000,0.00107043,705.951,1.8891,4.1938,1570.10000,46.06300,182.34000,711.86000],
        [500.000,170.000,0.00108062,747.979,1.9851,4.2121,1551.40000,44.54400,171.52000,709.53000],
        [500.000,180.000,0.0010914,790.199,2.0793,4.2324,1531.00000,42.95500,161.99000,706.30000],
        [500.000,190.000,0.0011028,832.635,2.1719,4.2551,1509.10000,41.30400,153.54000,702.51000],
        [500.000,200.000,0.00111486,875.311,2.2631,4.2806,1485.80000,39.60200,146.00000,697.88000],
        [500.000,210.000,0.00112763,918.256,2.3529,4.309,1460.90000,37.85600,139.23000,692.38000],
        [500.000,220.000,0.00114116,961.503,2.4415,4.3409,1434.70000,36.07400,133.12000,686.09000],
        [500.000,230.000,0.00115553,1005.09,2.529,4.3766,1407.00000,34.26600,127.57000,679.04000],
        [500.000,240.000,0.0011708,1049.05,2.6155,4.4165,1378.00000,32.43900,122.49000,671.26000],
        [500.000,250.000,0.00118708,1093.43,2.7012,4.4613,1347.70000,30.60100,117.81000,662.76000],
        [500.000,260.000,0.00120446,1138.29,2.7861,4.5115,1316.00000,28.75900,113.48000,653.57000],
        [500.000,270.000,0.00122307,1183.68,2.8704,4.5678,1283.10000,26.92300,109.45000,643.69000],
        [500.000,280.000,0.00124305,1229.67,2.9543,4.631,1249.00000,25.10000,105.66000,633.14000],
        [500.000,290.000,0.00126459,1276.33,3.0379,4.702,1213.70000,23.29800,102.08000,621.91000],
        [500.000,300.000,0.00128789,1323.74,3.1214,4.7819,1177.30000,21.52600,98.67100,610.02000],
        [500.000,310.000,0.00131322,1372.0,3.2049,4.8719,1139.90000,19.79000,95.40600,597.44000],
        [500.000,320.000,0.00134088,1421.22,3.2885,4.9736,1101.50000,18.09700,92.25400,584.19000],
        [500.000,330.000,0.00137128,1471.52,3.3726,5.0889,1062.10000,16.45300,89.18900,570.23000],
        [500.000,340.000,0.00140492,1523.05,3.4574,5.22,1021.70000,14.85900,86.18700,555.54000],
        [500.000,350.000,0.00144244,1575.98,3.543,5.37,980.05000,13.31800,83.22100,540.08000],
        [500.000,360.000,0.00148475,1630.63,3.63,5.562,936.04000,11.80200,80.26400,523.79000],
        [500.000,370.000,0.00153293,1687.35,3.7189,5.7866,891.79000,10.37600,77.29600,506.59000],
        [500.000,380.000,0.00158849,1746.51,3.8101,6.0531,846.77000,9.02780,74.29100,488.37000],
        [500.000,390.000,0.00165351,1808.6,3.9045,6.3776,801.09000,7.76230,71.22400,469.01000],
        [500.000,400.000,0.00173089,1874.31,4.0028,6.7781,755.14000,6.58890,68.06800,448.33000],
        [500.000,410.000,0.00182467,1944.47,4.1063,7.2713,709.70000,5.52070,64.80000,426.12000],
        [500.000,420.000,0.00194037,2020.07,4.2161,7.8636,666.13000,4.57370,61.40600,402.17000],
        [500.000,430.000,0.00208504,2101.99,4.3334,8.5275,626.51000,3.76510,57.90400,376.31000],
        [500.000,440.000,0.00226604,2190.53,4.4585,9.1599,593.57000,3.10960,54.37100,348.62000],
        [500.000,450.000,0.00248744,2284.44,4.5892,9.5672,569.96000,2.61200,50.97100,319.72000],
        [500.000,460.000,0.00274521,2380.52,4.7212,9.5776,556.72000,2.25800,47.92800,290.98000],
        [500.000,470.000,0.00302709,2474.69,4.8488,9.2032,552.52000,2.01700,45.40700,264.27000],
        [500.000,480.000,0.00331861,2563.86,4.968,8.6087,554.76000,1.85470,43.44800,241.00000],
        [500.000,490.000,0.00360858,2646.56,5.077,7.8969,560.89000,1.74360,41.98500,221.51000],
        [500.000,500.000,0.00388941,2722.52,5.1759,7.309,568.92000,1.66440,40.92400,205.98000],
        [500.000,510.000,0.00415994,2792.7,5.2661,6.7313,578.20000,1.60730,40.16300,193.40000],
        [500.000,520.000,0.00441748,2857.36,5.3482,6.2131,588.07000,1.56570,39.63300,183.36000],
        [500.000,530.000,0.00466236,2917.25,5.4232,5.7787,598.00000,1.53400,39.27600,175.39000],
        [500.000,540.000,0.00489567,2973.16,5.4924,5.4136,607.73000,1.50880,39.05100,169.04000],
        [500.000,550.000,0.00511848,3025.7,5.5566,5.1031,617.16000,1.48830,38.92500,163.96000],
        [500.000,560.000,0.00533183,3075.37,5.6166,4.8372,626.25000,1.47110,38.87600,159.88000],
        [500.000,570.000,0.00553667,3122.57,5.6729,4.6091,634.99000,1.45650,38.88900,156.62000],
        [500.000,580.000,0.0057339,3167.66,5.7261,4.4131,643.38000,1.44380,38.95100,154.02000],
        [500.000,590.000,0.00592433,3210.92,5.7765,4.244,651.44000,1.43270,39.05200,151.97000],
        [500.000,600.000,0.00610867,3252.61,5.8245,4.0973,659.20000,1.42270,39.18500,150.36000],
        [500.000,650.000,0.00695746,3443.48,6.0372,3.5873,694.26000,1.38560,40.16800,146.95000],
        [500.000,700.000,0.00771757,3614.76,6.218,3.2881,724.96000,1.36200,41.45200,147.95000],
        [500.000,750.000,0.00841745,3774.13,6.3777,3.1008,752.51000,1.34550,42.88300,152.05000],
        //[500.000,800.000,0.00907413,3925.96,6.5226,2.9813,777.37000,1.33190,44.39100,157.54000],

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



