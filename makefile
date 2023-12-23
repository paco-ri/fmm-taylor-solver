vacuum_example: vacuum_example.f90 ../virtual-casing/src/magneto-static-routs.f90 ../superconductor-type1/src/surf_routs.f90 geometry.f90 test_vacuum.f90 vacuum.f90 ../fmm3dbie/src/lap_wrappers/lap_comb_dir.f 
	gfortran vacuum_example.f90 -o vacuum_example ../virtual-casing/src/magneto-static-routs.f90 ../superconductor-type1/src/surf_routs.f90 geometry.f90 test_vacuum.f90 vacuum.f90 ../fmm3dbie/src/lap_wrappers/lap_comb_dir.f -L/home/paco/lib -lfmm3dbie -L/home/paco/lib -lfmm3d -fallow-argument-mismatch -fopenmp -fbacktrace -fcheck=all -g -ffpe-trap=invalid -O

B0_conv: B0_conv.f90 ../virtual-casing/src/magneto-static-routs.f90 ../superconductor-type1/src/surf_routs.f90 geometry.f90 test_vacuum.f90 vacuum.f90
	gfortran B0_conv.f90 -o B0_conv ../virtual-casing/src/magneto-static-routs.f90 ../superconductor-type1/src/surf_routs.f90 geometry.f90 test_vacuum.f90 vacuum.f90 -L/home/paco/lib -lfmm3dbie -L/home/paco/lib -lfmm3d -fallow-argument-mismatch -fopenmp	

vacuum_noblock_example: vacuum_noblock_example.f90 ../virtual-casing/src/magneto-static-routs.f90 ../superconductor-type1/src/surf_routs.f90 geometry.f90 test_vacuum.f90 vacuum.f90 vacuum_noblock.f90
	gfortran vacuum_noblock_example.f90 -o vacuum_noblock_example ../virtual-casing/src/magneto-static-routs.f90 ../superconductor-type1/src/surf_routs.f90 geometry.f90 test_vacuum.f90 vacuum.f90 vacuum_noblock.f90 -L/home/paco/lib -lfmm3dbie -L/home/paco/lib -lfmm3d -fallow-argument-mismatch -fopenmp -fcheck=all -g

Sprimematrix: Sprimematrix.f90 ../virtual-casing/src/magneto-static-routs.f90 ../superconductor-type1/src/surf_routs.f90 geometry.f90 test_vacuum.f90 vacuum.f90
	gfortran Sprimematrix.f90 -o Sprimematrix ../virtual-casing/src/magneto-static-routs.f90 ../superconductor-type1/src/surf_routs.f90 geometry.f90 test_vacuum.f90 vacuum.f90 -L/home/paco/lib -lfmm3dbie -L/home/paco/lib -lfmm3d -fallow-argument-mismatch -fopenmp -fcheck=all -g	
