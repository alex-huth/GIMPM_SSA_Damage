# GIMPM-SSA-Damage

The Generalized Interpolation Material Point Method (GIMPM) for the 
Shallow Shelf Approximation (SSA) of ice flow with Damage

Developer: Alex Huth (ahuth@princeton.edu)

This repository contains the GIMPM-SSA-Damage model detailed in:

      Huth, A., Duddu, R., Smith, B.E. (2020a). A generalized interpolation material point method for shallow ice shelves. Part I: shallow shelf approximation and ice thickness evolution. 

      Huth, A., Duddu, R., & Smith, B. E. (2020b). A generalized interpolation material point method for shallow ice shelves. Part II: Anisotropic creep damage mechanics and application to a marine ice sheet.
      
(manuscripts currently in review)

In addition to the GIMPM, this code also includes the standard Material Point Method (sMPM)

Damage models included (see Huth et al., 2020b for details):
  - SSA creep damage (Huth et al., 2020b)
  - SSA "zero-stress" damage (Sun et al., 2017)
  - The SSA "zero-stress" damage model + a modification to include necking/mass balance effects (Bassis & Ma, 2015)

## Compilation
See README in `PROG`, which contains the main source code.

## Test cases
The examples from Huth et al., 2020a are found in `test1d`, `test2d`, and `mismip/steady`
The examples from Huth et al., 2020b are found in `mismip/damage`
Each directory contains a README with instructions for running the examples.

## Notes
- An installation and knowledge of Elmer FEM and Elmer/Ice is required
  - https://github.com/ElmerCSC/elmerfem
  - https://elmerfem.org/elmerice/wiki/
- This model has not yet been parallelized