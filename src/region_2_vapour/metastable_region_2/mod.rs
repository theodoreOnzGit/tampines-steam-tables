/// ideal gas part
///
/// see footnotes for table 2.6
const METASTABLE_REGION_2_COEFFS_IDEAL: [[f64; 2]; 9] = [
    [0.0, -0.96937268393049e1],
    [1.0, 0.10087275970006e2],
    [-5.0, -0.56087911283020e-2],
    [-4.0, 0.71452738081455e-1],
    [-3.0, -0.40710498223928],
    [-2.0, 0.14240819171444e1],
    [-1.0, -0.43839511319450e1],
    [2.0, -0.28408632460772],
    [3.0, 0.21268463753307e-1],
];
/// residual part
///
/// I use underscores to make the demarcation clear for easy checking
const METASTABLE_REGION_2_COEFFS_RES: [[f64; 3]; 13] = [
    [1.0, 0.0, -0.733_622_601_865_06e-2],
    [1.0, 2.0, -0.882_238_319_431_46e-1],
    [1.0, 5.0, -0.723_345_552_132_45e-1],
    [1.0, 11.0, -0.408_131_785_344_55e-2],
    [2.0, 1.0, 0.200_978_033_802_07e-2],
    [2.0, 7.0, -0.530_459_218_986_42e-1],
    [2.0, 16.0, -0.761_904_090_869_70e-2],
    [3.0, 4.0, -0.634_980_376_573_13e-2],
    [3.0, 16.0, -0.860_430_930_285_88e-1],
    [4.0, 7.0, 0.753_215_815_227_70e-2],
    [4.0, 10.0, -0.792_383_754_461_39e-2],
    [5.0, 9.0, -0.228_881_607_784_47e-3],
    [5.0, 10.0, -0.264_565_014_828_10e-2],
];


/// metastable ideal gas correlations
pub mod metastable_ideal_gas_gamma;
pub use metastable_ideal_gas_gamma::*;

/// metastable residual correlations
pub mod metastable_residual_gamma;
pub use metastable_residual_gamma::*;

/// intensive properties in metastable region 2 
pub mod intensive_properties;
pub use intensive_properties::*;

