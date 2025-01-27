use uom::si::{available_energy::kilojoule_per_kilogram, f64::*, pressure::megapascal, thermodynamic_temperature::kelvin};

use crate::constants::P_C_MPA;
// assuming we are already in region 3
// calculate temperature given p and h
#[inline]
pub fn t_ph_3(p: Pressure, h: AvailableEnergy,) -> ThermodynamicTemperature {

    let is_region_3a = is_3a_when_in_region_3(p, h);

    if is_region_3a {
        t_ph_3a(p, h)
    } else {
        t_ph_3b(p, h)
    }

}

use super::h_3a3b_backwards_ph_boundary;
const T_PH_SUBREGION_3A_COEFFS: [[f64; 3]; 31] = [
    [-12.0,0.0,-0.133_645_667_811_215e-6],
    [-12.0,1.0,0.455_912_656_802_978e-5],
    [-12.0,2.0,-0.146_294_640_700_979e-4],
    [-12.0,6.0,0.639_341_312_970_080e-2],
    [-12.0,14.0,3.72_783_927_268_847E+02],
    [-12.0,16.0,-7.18_654_377_460_447E+03],
    [-12.0,20.0,5.73_494_752_103_400E+05],
    [-12.0,22.0,-2.67569329111439E+06],
    [-10.0,1.0,-3.34_066_283_302_614E-05],
    [-10.0,5.0,-2.45_479_214_069_597E-02],
    [-10.0,12.0,4.78_087_847_764_996E+01],
    [-8.0,0.0,7.64_664_131_818_904E-06],
    [-8.0,2.0,1.28_350_627_676_972E-03],
    [-8.0,4.0,1.71_219_081_377_331E-02],
    [-8.0,10.0,-8.51_007_304_583_213E+00],
    [-5.0,2.0,-1.36_513_461_629_781E-02],
    [-3.0,0.0,-3.84_460_997_596_657E-06],
    [-2.0,1.0,3.37_423_807_911_655E-03],
    [-2.0,3.0,-5.51_624_873_066_791E-01],
    [-2.0,4.0,7.29_202_277_107_470E-01],
    [-1.0,0.0,-9.92_522_757_376_041E-03],
    [-1.0,2.0,-1.19_308_831_407_288E-01],
    [0.0,0.0,7.93_929_190_615_421E-01],
    [0.0,1.0,4.54_270_731_799_386E-01],
    [1.0,1.0,2.09_998_591_259_910E-01],
    [3.0,0.0,-6.42_109_823_904_738E-03],
    [3.0,1.0,-0.235_155_868_604_540E-01],
    [4.0,0.0,2.52_233_108_341_612E-03],
    [4.0,3.0,-7.64885133368119E-03],
    [10.0,4.0,1.36176427574291E-02],
    [12.0,5.0,-1.33027883575669E-02],
    ];
/// eq 2.28
#[inline]
pub fn t_ph_3a(p: Pressure, h: AvailableEnergy) -> ThermodynamicTemperature {
    let p_ref = Pressure::new::<megapascal>(100.0);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(2300.0);
    let t_ref_kelvin = ThermodynamicTemperature::new::<kelvin>(760.0);

    let pi: f64 = (p/p_ref).into();
    let eta: f64 = (h/h_ref).into();

    // this is dimensionless temperature
    let mut theta = 0.0;

    for coeffs in T_PH_SUBREGION_3A_COEFFS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        theta += ni * (pi + 0.240).powi(ii as i32) * (eta - 0.615).powi(ji as i32);
    };

    return theta * t_ref_kelvin;

}

const T_PH_SUBREGION_3B_COEFFS: [[f64; 3]; 33] = [
    [-12.0,0.0,3.23_254_573_644_920E-05],
    [-12.0,1.0,-1.27_575_556_587_181E-04],
    [-10.0,0.0,-4.75851877356068E-04],
    [-10.0,1.0,1.56183014181602E-03],
    [-10.0,5.0,1.05724860113781E-01],
    [-10.0,10.0,-8.58514221132534E+01],
    [-10.0,12.0,7.24140095480911E+02],
    [-8.0,0.0,2.96475810273257E-03],
    [-8.0,1.0,-5.9272198336599E-03],
    [-8.0,2.0,-1.2630542281867E-02],
    [-8.0,4.0,-1.1571619636485E-01],
    [-8.0,10.0,8.4900096973960E+01],
    [-6.0,0.0,-1.0860226008662E-02],
    [-6.0,1.0,1.5430447532885E-02],
    [-6.0,2.0,7.5045544152447E-02],
    [-4.0,0.0,2.5252097361298E-02],
    [-4.0,1.0,-6.0250790123300E-02],
    [-3.0,5.0,-3.0762222135050E+00],
    [-2.0,0.0,-5.7401195986488E-02],
    [-2.0,4.0,5.0347136093985E+00],
    [-1.0,2.0,-9.2508188858483E-01],
    [-1.0,4.0,3.9173388291755E+00],
    [-1.0,6.0,-7.7314600713019E+01],
    [-1.0,10.0,9.4930876209859E+03],
    [-1.0,14.0,-1.4104371967941E+06],
    [-1.0,16.0,8.4916623081903E+06],
    [0.0,0.0,8.6109572944670E-01],
    [0.0,2.0,3.2334644281172E-01],
    [1.0,1.0,8.7328193602044E-01],
    [3.0,1.0,-4.3665304852668E-01],
    [5.0,1.0,2.8659671452948E-01],
    [6.0,1.0,-1.3177833127623E-01],
    [8.0,1.0,6.7668206433028E-03],
    ];

/// eq 2.29
#[inline]
pub fn t_ph_3b(p: Pressure, h: AvailableEnergy) -> ThermodynamicTemperature {
    let p_ref = Pressure::new::<megapascal>(100.0);
    let h_ref = AvailableEnergy::new::<kilojoule_per_kilogram>(2800.0);
    let t_ref_kelvin = ThermodynamicTemperature::new::<kelvin>(860.0);

    let pi: f64 = (p/p_ref).into();
    let eta: f64 = (h/h_ref).into();

    // this is dimensionless temperature
    let mut theta = 0.0;

    for coeffs in T_PH_SUBREGION_3B_COEFFS {
        let ii = coeffs[0];
        let ji = coeffs[1];
        let ni = coeffs[2];

        theta += ni * (pi + 0.298).powi(ii as i32) * (eta - 0.720).powi(ji as i32);
    };

    return theta * t_ref_kelvin;

}

/// see page 48
///
/// given we are in region 3 already, this eqn 
/// determines if we are in 3a or 3b
#[inline] 
pub(crate) fn is_3a_when_in_region_3(p: Pressure, h: AvailableEnergy) -> bool {

    // now below critical pressure (enthalpy)
    //
    // for 3a or 3b determination, we check if the 
    // enthalpy is above or below the appropriate boundary
    //
    let p_crit = Pressure::new::<megapascal>(P_C_MPA);

    // the dividing line if we are already in region 3 is 
    // in page 48 of the text, 
    // the boundary line h_3a3b_backwards_ph_boundary
    // belongs to region 3a

    let mut h_3ab_boundary = h_3a3b_backwards_ph_boundary(p);

    if p < p_crit {
        h_3ab_boundary = h_3a3b_backwards_ph_boundary(p_crit);
    };

    // so if enthalpy is greater than the boundary,
    // it is region 3b
    if h > h_3ab_boundary {
        return false;
    };

    // otherwise, it's region 3a
    return true;
}
