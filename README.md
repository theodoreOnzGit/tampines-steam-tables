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

