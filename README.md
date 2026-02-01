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
  git clone https://github.com/distribution-record-myosin/half_sarcomere.git
  cd half_sarcomere
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
    export OMP_NUM_THREADS=40;                                     
    ./euler_maruyama      
```
In Euler Maruyama method, time step width is 0.5 ns. In order to shorten the execution time, the simulation time length is set to 0.02 sec.

To prepare distribution records in the paper,  
run:
```
    export OMP_NUM_THREADS=10
    mpirun -np 10 ./record_gen
```
Distribution records are stored in the dist_records directory.
The distribution of myosin behavior is measured over 1,000 trials for each initial condition.
Parallelism is implemented using multi-threading for trials and multi-processing for the various initial conditions.


To reproduce the results of distribution record based method in the paper,  
run:
```
    export OMP_NUM_THREADS=40
    ./montecalro_w_record
```
In distribution record based method, time step width is set to be 1000 ns with distribution records.
The simulation time length is set to 1.0 sec.

The number of threads and processes can be adjusted to suit your computing environment.

------------------------------------------------------------
5. Reproducing the Results in the Paper
------------------------------------------------------------

- trans.csv has the reslts of each output time step.
    9-th column data corresponds to the z-line displacement.
  - Be careful not to overwrite the output file, as the file name is the same for both Euler Maruyama method and distribution record based method.
- To make a graph of z-line displacement "z_disp.png" from trans.csv,  
  run: 
```
   python z_disp.py
```

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
