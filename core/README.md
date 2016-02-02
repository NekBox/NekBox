# NekBox on Github
![https://travis-ci.org/maxhutch/NekBox](https://travis-ci.org/maxhutch/NekBox.svg?branch=box)


## NekBox

NekBox is a version of Nek5000 specialized to box geometries and intended for prototyping new methods.

### Constraints
 - Box geometry
 - Only `P`, `W`, `O`, or `SYM` boundaries
 - PN-PN only
 - A few more...

### Advantages
 - 2x less memory used, all dynamically allocated
 - Spectral coarse solver (very fast)
 - No re2 or map files
 - O(minute) startup
 - Pure Fortran90 (including usr file)
 - 2x smaller codebase
 
## Usage

NekBox uses `makenek` and `rea` files just like Nek5000.  The rea files don't need to include mesh information (never any `re2` files).  The `usr` file is compiled in just like in Nek5000, but is written in Fortran90 and references modules instead of common blocks.  There is no map file.

It is highly recommended to setup NekBox runs using [genrun script](https://www.github.com/maxhutch/nek-tools/).
