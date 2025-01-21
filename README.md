# TAMPINES Steam Tables
In house steam tables for the Thermo-hydraulic Artificial intelligence 
Multi-Phase INtegrated Emulator System (TAMPINES) solver

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

