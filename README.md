# fmm-taylor-solver

Compile using
```
gfortran vacuum.f90 -o vacuum ($PATH1)/virtual-casing/src/magneto-static-routs.f90 ($PATH2)/superconductor-type1/src/surf_routs.f90 -L($FMMBIE_INSTALL_DIR) -lfmm3dbie -L($FMM_INSTALL_DIR) -lfmm3d
```