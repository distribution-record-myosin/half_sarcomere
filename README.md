# Source Code for "ACCELERATION OF BIO-MOLECULE SIMULATION BASED ON OVER-DAMPED LANGEVIN EQUATION BY MONTE-CARLO METHOD"

This repository contains the source code used in the numerical experiments reported in the paper
  "ACCELERATION OF BIO-MOLECULE SIMULATION BASED ON OVER-DAMPED LANGEVIN EQUATION BY MONTE-CARLO METHOD"

by Akihiro Fujii and Takumi Washio, submitted to ....

------------------------------------------------------------
1. Overview
------------------------------------------------------------

This code implements the numerical methods and algorithms proposed in the
above-mentioned paper. It was used to generate all numerical results,
figures, and tables presented in the manuscript.

The implementation is intended to support reproducibility of the results
and to serve as a reference for further research.

------------------------------------------------------------
2. Requirements
------------------------------------------------------------

The code has been tested with the following environment:

- Programming language: Fortran 90
- Required libraries:
  - MPI
  - OpenMP

------------------------------------------------------------
3. Compilation
------------------------------------------------------------

Clone the repository:
```
  git clone https://github.com/USERNAME/REPOSITORY.git
  cd REPOSITORY
```
Edit Makefile:
```
    MPIF90=mpifrtpx # Fortran90 compiler with MPI                               
    F90=frtpx       # Fortran90 complier
    OPTFLG=-Kfast   # optimization flag     
    OMPFLG=-Kopenmp # OpenMP flag for the compiler          
```
    
Compile the source files:
```
   make
```

------------------------------------------------------------
4. Usage
------------------------------------------------------------

To reproduce the results of euler_maruyama method in the paper, 
run:
```
    export OMP_NUM_THREADS=40;                                     │
    ./euler_maruyama      
```


To prepare distribution records in the paper, 
run:
```
    export OMP_NUM_THREADS=10
    mpirun -np 10 ./record_gen
```

To reproduce the results of distribution record based method in the paper, 
run:
```
    export OMP_NUM_THREADS=40
    ./montecalro_w_record
```

------------------------------------------------------------
5. Reproducing the Results in the Paper
------------------------------------------------------------

- trans.csv has the reslts of each output time step.
    9-th column data corresponds to the z-line displacement.

------------------------------------------------------------
6. License
------------------------------------------------------------

This code is released under the MIT License.
See the LICENSE file for details.

------------------------------------------------------------
7. Contact
------------------------------------------------------------

For questions or issues related to this code, please contact:

.....
