# Source Code for "ACCELERATED MULTISCALE SIMULATION OF MOLECULAR MOTOR DYNAMICS USING A MONTE CARLO METHOD FOR OVERDAMPED LANGEVIN SYSTEMS"

This repository contains the source code used in the numerical experiments reported in the paper
  "ACCELERATED MULTISCALE SIMULATION OF MOLECULAR MOTOR DYNAMICS USING A MONTE CARLO METHOD FOR OVERDAMPED LANGEVIN SYSTEMS"

by Akihiro Fujii, Ryo Yoda and Takumi Washio, submitted to ....

## 1. Overview
This code implements the numerical methods and algorithms proposed in the above-mentioned paper. It was used to generate all numerical results, figures, and tables presented in the manuscript. 

The implementation is intended to support the reproducibility of the results and to serve as a reference for further research.

The repository consists of 3 core program components:
1. **The Euler-Maruyama (EM) method**: Direct stochastic simulation.
2. **The distribution record generation code**: Pre-computes transition probability distributions.
3. **The proposed Monte-Carlo (MC) method**: Accelerated simulation utilizing the pre-computed distribution records.

---

## 2. Requirements

The code has been tested and verified in the following environment:

- **Programming Language:** Fortran 90
- **Parallel Programming Libraries:**
  - MPI
  - OpenMP
- **Post-processing:** Python 3 (for data visualization)

---

## 3. Compilation

### 3.1. Clone the Repository
```bash
git clone https://github.com/distribution-record-myosin/half_sarcomere.git
cd half_sarcomere
```

### 3.2. Edit the Makefile
Open the `Makefile` and configure the compiler and flags according to your environment (the following is an example for the Fujitsu compiler):
```makefile
MPIF90=mpifrtpx # Fortran90 compiler with MPI                               
F90=frtpx       # Fortran90 compiler
OPTFLG=-Kfast   # Optimization flag     
OMPFLG=-Kopenmp # OpenMP flag for the compiler          
```

### 3.3. Compile
Run the following command to compile all source files:
```bash
make
```

> [!IMPORTANT]
> If you modify any simulation settings in `src/param.inc`, you must recompile the codes using `make` before running the simulation again.

---

## 4. Usage

> [!NOTE]
> * The number of threads and processes (`OMP_NUM_THREADS` and `-np`) should be adjusted to suit your computing environment.
> * Simulation length is set to `nt_out`$\times 0.1 \, ms$. Initially, `nt_out=1000`. 
> * When `ConstantVelocitySwitch = .true.` (Force-shortening velocity relationship), the executable requires a single input parameter (e.g., constant shortening velocity on the Z-line) via standard input.


### 4.1. Run Euler-Maruyama (EM) Method
To reproduce the EM simulation results:
```bash
export OMP_NUM_THREADS=40

# initially remove the output files
rm trans.csv data.bin

# If ConstantVelocitySwitch = .false. (Spontaneous oscillation)
./euler_maruyama      

# If ConstantVelocitySwitch = .true. (Force-shortening velocity relationship), constant shortening velocity is set as 100 nm/s
echo "100" | ./euler_maruyama
```
- **Time step size ($\Delta t$):** $0.5\, ns$
- **Output:** Z-line displacement and contraction force per actin are recarded in 9th and 10th columns of `trans.csv`. states of all myosins at each time are stored in `data.bin`. 

### 4.2. Generate Distribution Records for the MC Method
To pre-compute the distribution records required for the proposed MC method:
```bash
export OMP_NUM_THREADS=10
mpirun -np 10 ./record_gen
```
- **Output:** Generated distribution records will be stored in the `dist_records` directory.
- **Details:** The distribution of myosin behavior is measured over 15,000 trials for each initial condition. Parallelism is implemented using multi-threading for trials and multi-processing (MPI) for the various initial conditions.

### 4.3. Run Proposed Monte-Carlo (MC) Method
To reproduce the accelerated MC simulation results using the generated records:
```bash
export OMP_NUM_THREADS=40

# initially remove the output files
rm trans.csv data.bin

# If ConstantVelocitySwitch = .false. (Spontaneous oscillation)
./montecalro_w_record

# If ConstantVelocitySwitch = .true. (Force-shortening velocity relationship), constant shortening velocity is set as 100 nm/s
echo "100" | ./montecalro_w_record
```
- **Time step size ($\Delta t$):** $1,000\, ns$ (with distribution records).
- **Output:** Z-line displacement and contraction force per actin are recarded in 9th and 10th columns of `trans.csv`. states of all myosins at each time are stored in `data.bin`. 

## 5. Reproducing the Results in the Paper

Simulation configurations—such as boundary conditions, the time step size of the MC method, or power-stroke potential function shapes—are specified in `src/param.inc`.


### 5.1 Boundary conditions for the Z-line
You can switch the boundary condition by setting the value of `ConstantVelocitySwitch` in `src/param.inc`:

| Value | Boundary Condition |
| :--- | :--- |
| `.false.` | Spontaneous oscillation |
| `.true.` | Force-shortening velocity relationship |

### 5.2 Shapes of power-stroke potentials
The value of `DELTA_SCALE` in `src/param.inc` changes the parameter $\Delta$ of the power-stroke potential function:

| `DELTA_SCALE` | Potential Parameter |
| :--- | :--- |
| `2` | $\Delta = 2k_BT$ |
| `4` | $\Delta = 4k_BT$ |
| `8` | $\Delta = 8k_BT$ |

### 5.3 Time Step Size of the MC Method

The value of `time` in `src/param.inc` sets the time step size in nano-seconds:

| `time` | Time Step Size ($\Delta t$) |
| :--- | :--- |
| `1000` | $\Delta t = 1,000\text{ ns}$ |
| `5000` | $\Delta t = 5,000\text{ ns}$ |

### 5.4. Data Analysis and Visualization

- **Output File:** After executing either the EM or MC method, the results for each output time step ($0.1\, ms$) are stored in `trans.csv`.
- **Data Structure:** The 9th and 10th columns correspond to the Z-line displacement and contraction force per actin filament, respectively.

  - Be careful not to overwrite the output file, as the file name is the same for both EM method and the MC method.
- To make a graph of half-sarcomere shortening length "hssl.png" from trans.csv,
  run: 
```
   python3 hssl.py
```
- To make a graph "contract_force.png" of contraction force per actin filament from trans.csv
  run: 
```
   python3 contract_force.py
```
- To make a graph "dist.png" of distribution of the states of myosins bound to actin at 100 $\times 0.1\, ms$ from data.bin,  
  run: 
```
   # states of myosins bound to actin are recorded 
   # in dist.csv 
   echo "100" | ./get_dist
   python3 dist.py
```

## 6. License


This code is released under the MIT License.
See the LICENSE file for details.


## 7. Contact

For questions or issues related to this code, please contact:
- **Name / Email:** Akihiro Fujii / fujii@cc.kogakuin.ac.jp
- **Paper Authors:** Akihiro Fujii, Ryo Yoda, and Takumi Washio