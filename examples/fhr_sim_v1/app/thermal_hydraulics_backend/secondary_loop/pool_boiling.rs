// for this dummy DNB type model
//
// I specify a maximum UA, say 700,000 W/K as in the FHR educational simulator 
//
// The minimum UA is perhaps 1,000 W/K. It should start from a minimum value,
// head to maximum at say 
//
// 30 degrees of superheat, then dip back down to 1,000 W/K
//
// We can make a dummy model of bubble accumulation and growth in the wetted 
// portion.
//
// A typical curve could look like this:
//
// https://www.sciencedirect.com/topics/engineering/pool-boiling-curve
// 
// degree of superheat (degc),heat flux (W/m^2),heat transfer coeff (W/m^2 K)
// 1,1.00E+03,1.00E+03
// 5,1.00E+04,2.00E+03
// 12,7.00E+04,5.83E+03
// 30,1.00E+06,3.33E+04
// 320,1.00E+04,3.13E+01
// 2000,1.00E+06,5.00E+02
//
// Planted on log log plots, the heat transfer coeff vs degree of superheat 
// resembles that of a quantum harmonic oscillator (one level above ground 
// state)
// https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator
//
// that is x * e^(-a x^2)
//
//
// When plotting log10(degree of superheat) vs 
// log10 (heat transfer coeff)
//
// I make a variable x
//
// x_mod = log10(degree of superheat) - 2 
//
// apparently, this curve fits pretty well
// 
//
// y = -b * x_mod * exp (-a * x_mod * x_mod) + c
//
// These parameters reproduce log10(heat transfer coeff) quite okay. 
// At least fits to within the curve
// a = 2.4
// b = 6 
// c = 3
//
//
//
//

