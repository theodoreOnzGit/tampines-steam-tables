# TAMPINES Steam Tables
In house steam tables for the Thermo-hydraulic Artificial intelligence 
Multi-Phase INtegrated Emulator System (TAMPINES) solver.


This relies heavily upon the [Rust-steam](https://github.com/marciorvneto/rusteam)
library licensed using the MIT license. 

However, [Rust-steam](https://github.com/marciorvneto/rusteam) is incomplete 
for now. Moreover, it does not use the units of measure library. This 
set of steam-tables is meant to used dimensioned units by default. It will 
also incorporate verification tests from the following reference:

Kretzschmar, H. J., & Wagner, W. (2019). 
International steam tables. Springer Berlin Heidelberg.

Significant portions of code will be copied from the rust-steam package.
Hence, I am putting the rust-steam license here.

# Changelog 

v0.1.1 

Starting the h,s flash algorithms.

Also, copied some openfoam algorithms which will form the basis for which 
the steam tables in tampines are used to solve two phase flow in transient 
scenarios. Since OpenFOAM is licensed under 
GNU GPLv3, tampines-steam-tables will also be licensed under GNU GPL v3.

Note that for points near boundaries, correction factors have not been 
applied for (p,h) and (p,s) flashes. (h,s) flashes have only been 
partly implemented. 

v0.1.0 

Added dielectric constant and surface tension functions.
Didn't yet test across the whole steam table, but it works for the 
small unit tests.

v0.0.9
Implemented thermal conductivity for ps flash. However, for 160 bar 
and 220 bar steam tables, the max error is 30 and 40% respectively.
For the 220 bar steam tables, it is quite near critical point, 
so thermal conductivity equations were not meant to be accurate there.
However, for 160 bar, it is sufficiently far from critical point that this 
shouldn't be the case. But i'm leaving it as such for now, till such time 
there is a better reference for such properties.

Near critical temperatures and pressures, eg. 180 bar about 
357+ degrees C, the speed of sound, isentropic exponent
the specific heat capacity are not accurate to within 1%. Some discrepancies 
are larger than 5-10% esp near critical region for these properties.

The thermal conductivity and dynamic viscosity also have similar-ish degrees 
of uncertainty.


v0.0.8
Implemented thermal conductivity for ph flash. However, for 160 bar 
and 220 bar steam tables, the max error is 30 and 40% respectively.
For the 220 bar steam tables, it is quite near critical point, 
so thermal conductivity equations were not meant to be accurate there.
However, for 160 bar, it is sufficiently far from critical point that this 
shouldn't be the case.

I'm quite puzzled as to why that is the case. But this is a bug that needs 
to be fixed.

v0.0.7

Starting development of an interface for the forward and backward flash.
First using functional programming, then object oriented programming.
OOP not implemented in this version yet.

Firstly, (p,T) flash, then (p,H) flash.
This is done for all steam table values except for 1000 bar. 
There is a problem for (p,T) and (p,H) flashing at 1000 bar 
as the algorithm complains it's out of range. Will take some debugging 
to settle.

Dynamic viscosity added, can reproduce steam table to within 2%.
Thermal conductivity can reproduce steam table values to within 1% 
except for regions near critical point (220 bar) and from 100 bar 
up to 200 bar.
Thermal conductivity yet to be implemented for (p,H) flash,
requires some debugging.




v0.0.6
beginning the addition of backwards eqns

First, pressure and enthalpy (p,h) flash. 
This is applicable for 

region 1, which forward equations are (p,t) flash:
- T(p,h)

region 2, which forward eqns are (p,t) flash:
- T(p,h) 

region 3, which forward equations are (v,t) flash, so it accounts for quality:
- T(p,h)
- V(p,h)

once (p,h) flash is done for regions 1,2 and 3, then you can get T,
for region 1 and 2 or (V,T) for region 3

and then get all your other thermodynamic variables

the ps3 equations (enthalpy to pressure equations) are also added

However, the interface for an overall ph flash or tp flash is not yet 
available.


v0.0.5
Add Region 5 equations (no backwards equations here)

v0.0.4 
Added region 4 vapour liq saturation temp and pressure 
line up to critical point. This includes triple point, 
normal boiling point (100C at 1 atm) and critical point of water.

v0.0.3 
Added region 3, and the saturation temperature and pressure boundary equation 
p23 and b23.
Only forward eqns added. That is (T,P) flash.

v0.0.2

Added region 2, including metastable, dimensioned equations, with verification tests.
Only forward eqns added. That is (T,P) flash.

v0.0.1 

Added region 1 dimensioned equations, with verification tests.
Only forward eqns added. That is (T,P) flash.

## Rust-steam license:

Copyright 2023

Permission is hereby granted, free of charge, to any person obtaining a 
copy of this software and associated documentation files (the “Software”), 
to deal in the Software without restriction, including without limitation 
the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, 
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

