This repository contains the HiRep simulation code.

# Getting started

## Dependencies

* GCC or different C-compiler
* MPI implementation, i.e. OpenMPI or MPICH for MPI support
* In order to run make use of CUDA GPU acceleration compile with CUDA 11.x, for multi-GPU using a CUDA-aware MPI implementation


## Compilation

### Clone the directory

```
git clone https://github.com/claudiopica/HiRep
```

### Adjust Make Flags
Adjust the file ```Make/MkFlags``` for the right compilation

* Gauge group SU(NG) or SO(NG)
```
NG=3

#CHOICES ARE GAUGE_SUN AND GAUGE_SON
GAUGE_GROUP = GAUGE_SUN
#GAUGE_GROUP = GAUGE_SON
```

* Representation
```
REPR = REPR_FUNDAMENTAL
#REPR = REPR_SYMMETRIC
#REPR = REPR_ANTISYMMETRIC
#REPR = REPR_ADJOINT
```

* Boundary Conditions

Uncomment the line here, when you want to establish certain boundary conditions into the respective direction.
```
#T => PERIODIC, ANTIPERIODIC, OPEN, THETA
#X => PERIODIC, ANTIPERIODIC, THETA
#Y => PERIODIC, ANTIPERIODIC, THETA
#Z => PERIODIC, ANTIPERIODIC, THETA
```

* Macro parameters

Then a number of macro parameters follow. Here you have to specify if you want to compile for certain boundary conditions by adding the identifier to the ```MACRO``` variable.

```
#MACRO += -DBC_T_THETA
#MACRO += -DBC_T_PERIODIC
MACRO += -DBC_T_ANTIPERIODIC
#MACRO += -DBC_T_OPEN
MACRO += -DBC_X_PERIODIC
MACRO += -DBC_Y_PERIODIC
MACRO += -DBC_Z_PERIODIC

#MACRO += -DBC_XYZ_TWISTED

#MACRO += -DHALFBG_SF
#MACRO += -DBASIC_SF
#MACRO += -DROTATED_SF
```

Specify, whether you want to compile with MPI either with or without GPU acceleration by using 

```
#MACRO += -DWITH_MPI
```

For compilation with GPU acceleration for CUDA GPUs, add the identifier ```-DWITH_GPU``` to ```MACRO```.

```
MACRO += -DWITH_GPU
```

* Compilers, wrappers, preprocessors

A number of example combinations are already given in ```MkFlags```.

For compiling with ```GCC``` and ```OpenMPI``` one would compile with

```
CC = gcc
MPICC = mpicc
LDFLAGS =
INCLUDE =
```

Using Intel compilers and Intel's MPI implementation, one can use for example

```
CC = icc
MPICC = mpiicc
LDFLAGS = -L /opt/local/lib/mpich-devel-gcc7/ -L /opt/local/lib/
INCLUDE = -I /opt/local/include/mpich-devel-gcc7/
```

For CUDA acceleration, use ```nvcc``` and adjust the flag ```-arch``` according to the compute capability of the CUDA capable device.

```
CC = nvcc
CFLAGS = -O2 -Xcompiler '-std=c99 -fgcse-sm -fgcse-las -fgcse-after-reload'
GPUFLAGS = --x cu -arch=sm_70 -Xptxas -v -Xptxas -dlcm=ca -dc
LDFLAGS = -lcuda
```

For compilation with CUDA-aware MPI one needs to pass the MPI wrapper of the MPI implementation to the CUDA preprocessor using the flag ```-ccbin```. For OpenMPI this is ```mpicc```.

```
CC = nvcc
CFLAGS = -ccbin mpicc -Xcompiler '-std=c99 -fgcse-sm -fgcse-las -fgcse-after-reload'
GPUFLAGS = --x cu -arch=sm_70 -Xptxas -v -Xptxas -dlcm=ca -dc
LDFLAGS = -lcuda -lmpi
```

## Run

### Adjust input file

Compile the HiRep library for example in ```LibHR``` by typing ```make```. An example of a C-file that generates a binary to run the HMC can be found in ```HMC```, you can navigate into this directory and type ```make``` to create a binary. It is necessary to specify a number of parameters using an input file, see ```HMC/input_file``` for an example. For basic run variables, one can have a look at the section ```Run control variables```.

```
run name = run1
save freq = 1
meas freq = 1
conf dir = cnfg
gauge start = random 
last conf = +1
```

The "+" in front of ```last conf``` specifies the number of trajectories to be generated after the chosen startup configuration. I.e. if the startup configuration is trajectory number 5 and ```last conf = 6``` then one additional trajectory will be generated. If ```last conf = +6``` then six additional trajectories will be generated.

### Execute Binary

Run the HMC using 

```
$ hmc -i input_file
```

where ```hmc``` is the binary generated from ```hmc.c```.


# Documentation

* Development Notes [html](docs/_build/html/index.html) [pdf](docs/_build/latex/hirep-documentation.pdf)





![https://github.com/claudiopica/HiRep/actions?workflow=no-mpi](https://github.com/claudiopica/HiRep/workflows/no-mpi/badge.svg)
![https://github.com/claudiopica/HiRep/actions?workflow=mpi](https://github.com/claudiopica/HiRep/workflows/mpi/badge.svg)
