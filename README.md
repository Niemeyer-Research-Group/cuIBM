cuIBM - A GPU-based Immersed Boundary Method code
=================================================

This is a fork of the [Barba group](https://github.com/barbagroup)'s
[cuIBM](https://github.com/barbagroup/cuIBM).
This version of cuIBM has been tested on CentOS 6.7 with CUDA 7.5 and cusp 0.5.1.

### New Features:

Installation instructions
-------------------------

### Dependencies

Please ensure that the following dependencies are installed before compiling
cuIBM:

* Git distributed version control system (`git`)
* GNU C++ Compiler(`g++`)
* NVIDIA's CUDA Compiler (`nvcc`)
* CUSP Library (available [here](https://github.com/cusplibrary/cusplibrary))

#### Git (`git`)
Check if git is installed. On Ubuntu, this can be done via the Terminal using the following
command:

    >git
    
If it isn't installed you will get the message "The program 'git' is currently not installed."
To install `git`. 

    > sudo apt-get install git-core

#### GNU C++ Compiler (`g++`)

Install `g++` using the following command:

    > sudo apt-get install g++

Check the version of G++ installed:

    > g++ --version

The default version of g++ on Ubuntu 14.04 is 4.8.4.

#### NVIDIA's CUDA Compiler (`nvcc`)

[Download and install](https://developer.nvidia.com/cuda-downloads) the CUDA
Toolkit. cuIBM has been developed and tested with CUDA 7.5.

Check the version of NVCC installed:

    > nvcc --version

#### CUSP Library

CUSP is a library that provides routines to perform sparse linear algebra
computations on Graphics Processing Units. It is written using the CUDA
programming language and runs code on compatible NVIDIA GPUs.

CUSP is currently hosted on
[GitHub](https://github.com/cusplibrary/cusplibrary). cuIBM has been tested
and works with version 0.5.1, available for download
[here](https://github.com/cusplibrary/cusplibrary/archive/0.5.1.zip).

The instructions here assume that the CUSP library is to be installed in the
folder `/scratch/src/lib`, but any other folder with write permissions can be used.
Create a local copy of the CUSP library using the following commands:

    > mkdir -p $/scratch/src/lib
    > cd /scratch/src/lib

Download cusp 0.5.1 from [here](https://github.com/cusplibrary/cusplibrary/releases/tag/v0.5.1).
Copy the zip to /scratch/src/lib then:

    > unzip 0.5.1.zip

The folder `/scratch/src/lib/cusplibrary-0.5.1` is now created.

### Compiling cuIBM

This version of cuIBM can be found at its [GitHub repository](https://github.com/Niemyer-Research-Group/cuIBM).

Run the following commands to create a local copy of the repository in the
folder `/scratch/src` (or any other folder with appropriate read/write/execute
permissions):

    > cd /scratch/src
    > git clone https://github.com/Niemeyer-Research-Group/cuIBM.git

To compile, set the environment variable `CUSP_DIR` to point to the directory
with the CUSP library. For a `bash` shell, add the following line to the file
`~/.bashrc` or `~/.bash_profile`:

    export CUSP_DIR=/path/to/cusp/folder

which for the present case would be `$HOME/lib/cusplibrary-0.5.1`.

We also recommend setting the environment variable `CUIBM_DIR` to point to the
location of the cuIBM folder. While the code can be compiled and run without
setting this variable, some of the validation scripts provided make use of it.

    export CUIBM_DIR=/path/to/cuIBM/folder

which is `/scratch/src/cuIBM`, as per the above instructions.

Reload the file:

    > source ~/.bashrc

Switch to the cuIBM directory:

    > cd /scratch/src/cuIBM/src

Compile all the files:

    > make

Run a test:

    > make lidDrivenCavityRe100

**IMPORTANT**: If your NVIDIA card supports only compute capability 1.3, buy a new computer.  
Otherwise, edit Line 13 of the file `Makefile.inc` in the cuIBM root directory before
compiling: replace `compute_20` with `compute_13`.

Numerical schemes
-----------------

### FADLUN

### LUO

### Luo_iter

### Temporal discretisation
The convection terms are calculated using Crank-Nicolson and the advection terms are calculated using 2nd order explicit Adams-Bashforth.

### Spatial discretisation

The convection terms are calculated using a conservative symmetric
finite-difference scheme, and the diffusion terms are calculated using a
second-order central difference scheme.

Examples
--------

The following examples are available:

lidDrivenCavityRe + ###
Example: lidDrivenCavityRe100
* `100`: Flow in a lid-driven cavity with Reynolds number
100.
* `1000`: Flow in a lid-driven cavity with Reynolds number
1000.
* `10000`: Flow in a lid-driven cavity with Reynolds number
10000.
cylinderRe + ###
Example : cylinderRe40
* `40`: Flow over a circular cylinder at Reynolds number 40. The
flow eventually reaches a steady state.
* `75`: Flow over a circular cylinder at Reynolds number 75. The
initial flow field has an asymmetric perturbation that triggers instability in
the flow and vortex shedding is observed in the wake.
* `100`: Flow over a circular cylinder at Reynolds number 100. The
initial flow field has an asymmetric perturbation that triggers instability in
the flow and vortex shedding is observed in the wake.
* `150`: Flow over a circular cylinder at Reynolds number 150. The
initial flow field has an asymmetric perturbation that triggers instability in
the flow and vortex shedding is observed in the wake.
* `550`: Initial flow over an impulsively started cylinder at
Reynolds number 550.
* `3000`: Initial flow over an impulsively started cylinder at
Reynolds number 3000.
Static flow with oscillating cylinder.
Impulsivly started oscillating cylinder.
Vorticity induced vibrations

### Run the tests

Post-processing
---------------


Known issues
------------


Contact
-------

Please e-mail [Chris Minar](mailto:minarc@oregonstate.edu) if you have any
questions, suggestions or feedback.
