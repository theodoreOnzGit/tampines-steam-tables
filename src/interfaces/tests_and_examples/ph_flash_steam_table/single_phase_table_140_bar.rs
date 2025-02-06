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
        //[140.000,0.000,0.000993243,14.0723,0.0004338,4.1545,1424.7000,145.9600,1762.6000,565.9100],
        [140.000,2.000,0.00099324,22.3784,0.030732,4.1518,1434.5000,147.9800,1649.2000,570.5800],
        [140.000,4.000,0.000993298,30.6799,0.060793,4.1497,1443.9000,149.9200,1546.9000,575.0700],
        [140.000,6.000,0.000993411,38.9776,0.090625,4.1481,1452.9000,151.7800,1454.5000,579.4000],
        [140.000,8.000,0.000993577,47.2724,0.12023,4.1468,1461.6000,153.5700,1370.6000,583.5800],
        [140.000,10.000,0.000993794,55.565,0.14962,4.1458,1469.8000,155.2700,1294.2000,587.6200],
        [140.000,12.000,0.000994058,63.8559,0.1788,4.1451,1477.7000,156.9000,1224.5000,591.5200],
        [140.000,14.000,0.000994366,72.1456,0.20777,4.1446,1485.2000,158.4500,1160.5000,595.3100],
        [140.000,16.000,0.000994718,80.4344,0.23654,4.1442,1492.3000,159.9200,1101.8000,598.9800],
        [140.000,18.000,0.000995112,88.7225,0.2651,4.1439,1499.1000,161.3200,1047.7000,602.5400],
        [140.000,20.000,0.000995544,97.0102,0.29347,4.1438,1505.6000,162.6400,997.7700,606.0000],
        [140.000,25.000,0.00099679,117.729,0.36355,4.1438,1520.4000,165.6400,888.4100,614.2200],
        [140.000,30.000,0.000998255,138.448,0.43247,4.1441,1533.3000,168.2200,797.1700,621.8900],
        [140.000,35.000,0.000999923,159.17,0.50027,4.1448,1544.4000,170.3800,720.1900,629.0500],
        [140.000,40.000,0.00100178,179.897,0.56699,4.1458,1553.9000,172.1600,654.5900,635.7300],
        [140.000,45.000,0.00100382,200.629,0.63267,4.1471,1561.9000,173.5800,598.2100,641.9600],
        [140.000,50.000,0.00100603,221.368,0.69735,4.1487,1568.4000,174.6600,549.3700,647.7700],
        [140.000,55.000,0.0010084,242.116,0.76106,4.1507,1573.7000,175.4100,506.7800,653.1700],
        [140.000,60.000,0.00101094,262.875,0.82385,4.153,1577.7000,175.8700,469.4000,658.1600],
        [140.000,65.000,0.00101363,283.647,0.88573,4.1557,1580.5000,176.0300,436.4200,662.7700],
        [140.000,70.000,0.00101647,304.433,0.94675,4.1588,1582.3000,175.9300,407.1700,667.0100],
        [140.000,75.000,0.00101945,325.235,1.0069,4.1623,1583.0000,175.5700,381.1000,670.8800],
        [140.000,80.000,0.00102258,346.056,1.0663,4.1662,1582.7000,174.9800,357.7800,674.3900],
        [140.000,85.000,0.00102586,366.898,1.1249,4.1705,1581.6000,174.1700,336.8200,677.5600],
        [140.000,90.000,0.00102928,387.762,1.1828,4.1752,1579.6000,173.1500,317.9300,680.3800],
        [140.000,95.000,0.00103284,408.65,1.2399,4.1803,1576.7000,171.9300,300.8300,682.8600],
        [140.000,100.000,0.00103654,429.566,1.2963,4.1859,1573.1000,170.5400,285.3100,685.0200],
        [140.000,110.000,0.00104437,471.485,1.4072,4.1983,1563.8000,167.2500,258.2900,688.3900],
        [140.000,120.000,0.00105278,513.539,1.5155,4.2127,1551.7000,163.3700,235.6500,690.5200],
        [140.000,130.000,0.00106179,555.745,1.6215,4.229,1537.2000,158.9700,216.4900,691.4900],
        [140.000,140.000,0.00107141,598.127,1.7254,4.2476,1520.4000,154.1100,200.1200,691.3300],
        [140.000,150.000,0.00108168,640.707,1.8272,4.2688,1501.4000,148.8700,186.0200,690.1000],
        [140.000,160.000,0.00109263,683.512,1.9272,4.2928,1480.4000,143.2700,173.7700,687.8300],
        [140.000,170.000,0.00110432,726.573,2.0255,4.32,1457.3000,137.3700,163.0500,684.9100],
        [140.000,180.000,0.00111679,769.925,2.1222,4.3511,1432.2000,131.2000,153.6000,680.9700],
        [140.000,190.000,0.00113013,813.609,2.2176,4.3864,1405.2000,124.8000,145.2000,676.0600],
        [140.000,200.000,0.00114442,857.671,2.3117,4.4269,1376.2000,118.2000,137.6900,670.2100],
        [140.000,210.000,0.00115975,902.166,2.4048,4.4732,1345.1000,111.4400,130.9100,663.4700],
        [140.000,220.000,0.00117627,947.158,2.4969,4.5265,1312.0000,104.5300,124.7600,655.8500],
        [140.000,230.000,0.00119411,992.722,2.5884,4.588,1276.8000,97.5210,119.1200,647.3500],
        [140.000,240.000,0.00121347,1038.95,2.6794,4.6595,1239.5000,90.4340,113.9200,637.9700],
        [140.000,250.000,0.0012346,1085.95,2.7701,4.7431,1199.9000,83.2980,109.0800,627.7100],
        [140.000,260.000,0.00125779,1133.86,2.8608,4.8419,1157.9000,76.1360,104.5200,616.5300],
        [140.000,270.000,0.00128346,1182.85,2.9518,4.96,1113.2000,68.9640,100.1900,604.4100],
        [140.000,280.000,0.00131214,1233.15,3.0436,5.1037,1065.4000,61.7890,96.0360,591.2900],
        [140.000,290.000,0.00134461,1285.05,3.1366,5.2826,1013.9000,54.6090,91.9870,577.0700],
        [140.000,300.000,0.00138198,1338.97,3.2315,5.5128,957.7300,47.4090,87.9800,561.6300],
        [140.000,310.000,0.00142603,1395.56,3.3293,5.8226,895.5900,40.1760,83.9310,544.7200],
        [140.000,320.000,0.00147974,1455.87,3.4319,6.2683,826.5000,32.9740,79.7230,525.9900],
        [140.000,330.000,0.00154886,1521.8,3.5421,6.983,748.3800,25.8290,75.1630,504.7700],
        [140.000,336.669,0.00160971,1570.88,3.623,7.8117,682.8400,20.6900,71.7320,488.7400],
        [140.000,336.669,0.0114889,2638.09,5.373,11.26,445.1800,1.2322,22.1350,108.1100],
        [140.000,340.000,0.0119989,2672.38,5.4291,9.4878,457.1600,1.2442,22.2400,101.2500],
        [140.000,350.000,0.0132316,2752.92,5.5595,7.0059,482.9000,1.2588,22.6320,90.6220],
        [140.000,360.000,0.0142288,2816.39,5.6605,5.7978,502.2800,1.2665,23.0680,85.3110],
        [140.000,370.000,0.0150919,2870.38,5.7452,5.0539,518.4300,1.2721,23.5200,82.1560],
        [140.000,380.000,0.0158666,2918.26,5.819,4.5496,532.4000,1.2761,23.9760,80.1540],
        [140.000,390.000,0.0165779,2961.83,5.8853,4.184,544.8000,1.2789,24.4330,78.8600],
        [140.000,400.000,0.017241,3002.23,5.9457,3.9064,556.0200,1.2808,24.8890,78.0450],
        [140.000,410.000,0.0178661,3040.16,6.0017,3.6884,566.3100,1.2822,25.3420,77.5760],
        [140.000,420.000,0.0184603,3076.14,6.0539,3.5132,575.8500,1.2831,25.7920,77.3690],
        [140.000,430.000,0.019029,3110.53,6.1032,3.3697,584.7900,1.2837,26.2390,77.3690],
        [140.000,440.000,0.0195762,3143.61,6.1499,3.2504,593.2100,1.2840,26.6830,77.5370],
        [140.000,450.000,0.0201049,3175.6,6.1945,3.15,601.2000,1.2841,27.1230,77.8440],
        [140.000,460.000,0.0206177,3206.66,6.2371,3.0648,608.8100,1.2841,27.5610,78.2680],
        [140.000,470.000,0.0211167,3236.94,6.2782,2.992,616.0900,1.2839,27.9950,78.7910],
        [140.000,480.000,0.0216034,3266.54,6.3177,2.9293,623.0800,1.2836,28.4270,79.4020],
        [140.000,490.000,0.0220794,3295.55,6.356,2.875,629.8100,1.2832,28.8550,80.0890],
        [140.000,500.000,0.0225457,3324.06,6.3931,2.8278,636.3100,1.2828,29.2810,80.8430],
        [140.000,510.000,0.0230034,3352.13,6.4292,2.7867,642.6000,1.2822,29.7050,81.6580],
        [140.000,520.000,0.0234533,3379.81,6.4643,2.7506,648.7100,1.2816,30.1260,82.5260],
        [140.000,530.000,0.0238961,3407.15,6.4986,2.719,654.6400,1.2810,30.5440,83.4420],
        [140.000,540.000,0.0243326,3434.2,6.532,2.6911,660.4200,1.2803,30.9600,84.4020],
        [140.000,550.000,0.0247632,3460.99,6.5648,2.6666,666.0500,1.2796,31.3750,85.4030],
        [140.000,560.000,0.0251885,3487.54,6.5968,2.6449,671.5500,1.2789,31.7870,86.4400],
        [140.000,570.000,0.0256089,3513.89,6.6283,2.6256,676.9300,1.2781,32.1970,87.5100],
        [140.000,580.000,0.0260247,3540.06,6.6591,2.6086,682.2000,1.2774,32.6050,88.6110],
        [140.000,590.000,0.0264364,3566.07,6.6894,2.5936,687.3700,1.2766,33.0110,89.7410],
        [140.000,600.000,0.0268442,3591.94,6.7192,2.5803,692.4300,1.2758,33.4160,90.8970],
        [140.000,650.000,0.0288338,3719.67,6.8615,2.5336,716.4900,1.2717,35.4130,97.0060],
        [140.000,700.000,0.0307586,3845.69,6.9944,2.5103,738.8300,1.2676,37.3730,103.5300],
        [140.000,750.000,0.0326354,3970.94,7.12,2.5017,759.8300,1.2636,39.2990,110.3900],
        //[140.000,800.000,0.0344758,4096.02,7.2393,2.5033,779.7200,1.2596,41.1950,117.4700],

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
        max_relative=4e-3
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
        max_relative=5e-4
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
        max_relative=2e-2
        );



}



