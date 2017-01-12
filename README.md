cuIBM - A GPU-based Immersed Boundary Method code
=================================================

This is a fork of the [Barba group](https://github.com/barbagroup)'s
[cuIBM](https://github.com/barbagroup/cuIBM).
This version of cuIBM has been tested on Ubuntu 14.04 with CUDA 7.5 and cusp 0.5.1.

### New Features:

Installation instructions
-------------------------
The following is an installation guide that begins at a fresh copy of Ubuntu 14.04 and ends with some post processing.

### Dependencies
Please ensure that the following dependencies are installed before compiling
cuIBM:

* Git distributed version control system (`git`)
* GNU C++ Compiler(`g++`)
* NVIDIA's CUDA Compiler (`nvcc`)
* CUSP Library (available [here](https://github.com/cusplibrary/cusplibrary))
* Python/numpy
* postprocessthing

#### Git (`git`)
Check if git is installed. On Ubuntu, this can be done via the Terminal using the following
command:

    > git
    
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
Toolkit. cuIBM has been developed and tested with CUDA 7.5. Make sure the enviornment variables are set as described in the linux CUDA installation guide. If you use the .deb installer it should automatiaclly update the graphics card driver. If not, you may have to update the graphics driver.

Check the version of NVCC installed:

    > nvcc --version

#### CUSP Library

CUSP is a library that provides routines to perform sparse linear algebra
computations on Graphics Processing Units. It is written using the CUDA
programming language and runs code on compatible NVIDIA GPUs.

Create a local copy of the CUSP library using the following commands:

    > sudo mkdir /scratch
    > chmod o+rwx /scratch
    > cd /scratch
    > mkdir /src
    > cd src
    > mkdir lib
    > cd lib

Download cusp 0.5.1 from [here](https://github.com/cusplibrary/cusplibrary/releases/tag/v0.5.1).
Copy the zip to /scratch/src/lib then:

    > unzip cusplibrary-0.5.1.zip

The folder `/scratch/src/lib/cusplibrary-0.5.1` is now created.

### Python Libraries

If you want to use the post processing scripts, install these python libraries. You only need numpy and matplotlib.

    > sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose

### Compiling cuIBM

This version of cuIBM can be found at its [GitHub repository](https://github.com/Niemeyer-Research-Group/cuIBM).

Run the following commands to create a local copy of the repository in the
folder `/scratch/src` (or any other folder with appropriate read/write/execute
permissions):

    > cd /scratch/src
    > git clone https://github.com/Niemeyer-Research-Group/cuIBM.git

To compile, set the environment variable `CUSP_DIR` to point to the directory
with the CUSP library. For a `bash` shell, add the following lines to the file
`~/.bashrc` or `~/.bash_profile`:

    export CUSP_DIR=/scratch/src/lib/cusplibrary-0.5.1
    export CUIBM_DIR=/scratch/src/cuIBM

Reload the file:

    > source ~/.bashrc

Switch to the cuIBM directory:

    > cd /scratch/src/cuIBM/src

Compile all the files:

    > make

Run a test:

    > make cylinderRe40
    
View the drag profile:

    > cd /scratch/src/cuIBM
    > scripts/validation/validateCylinder.py
    
The output will be in `/scratch/src/cuIBM/validation/cylinder/Re40`

**IMPORTANT**: If your NVIDIA card supports only compute capability 1.3, buy a new computer.  
Otherwise, edit Line 13 of the file `Makefile.inc` in the cuIBM root directory before
compiling: replace `compute_20` with `compute_13`.

### Setting up Nvidia Nsight
The installtion of CUDA comes with an a variant of eclipse called NSight, an IDE designed to help debug CUDA kernels.
These steps will import the code into Nsight.  
1. Click `File->Import...`  
2. Select the folder `Git` and the `Projects from Git` inside that. Click `next`.  
3. Click `Existing local repository` then `next`  
4. Click `add` then change the directory to `/scratch/src/cuIBM` then click `search` and add the git project.  
5. Click `cuIBM` then `next`  
6. Click `Use the New PRoject wizard`  
7. Select `Makefile project with existing code`  
8. Set the projet name to `cuIBM` and the existing code location to `/scratch/src/cuIBM/src`  
9. Select CUDA Toolkit 7.5 in the toolchain box and hit `finish`.  


Immersed Boundary Methods
-------------------------

### Modified Fadlun
Based on the work of Fadlun et al. Modifications to prevent mass destruction at immersed boundary.

### External
Based on the work of Luo et al. Operations take place outside of linear algebra

### Embedded
Based on the work of Luo et al. Operations take place inside the linear algebra


Numerical Schemes
-----------------

### Temporal discretisation
The convection terms are calculated using Crank-Nicolson and the advection terms are calculated using 2nd order explicit Adams-Bashforth.

### Spatial discretisation

The convection terms are calculated using a conservative symmetric
finite-difference scheme, and the diffusion terms are calculated using a
second-order central difference scheme.

Examples
--------

Try some of the followng examples:

## Shorter
### Lid driven cavity
Example: make lidDrivenCavityRe100
Reynolds numbers availible: 100, 1000

### Impulsively started cylinder
Example: make cylinderRe40
Reynolds numbers availible: 40, 550, 3000

## Longer
For these simulations try switching between the external and embedded methods. The embedded method generally takes 10x longer. To do this, navigate to the simParams.yaml file for the simulation you want to modify. For example 'validaiton/osc/VIV/Ured4/simParams.yaml'.
Change the line '  SolverType: xxx'. Set xxx to LUO for the external method or OSC_CYLINDER for the embedded emthod.

### In-line oscillating cylinder driven flow
Example: make oscStaticExternal10

### Vorticity induced vibration
Example: make vivUred4

Post-processing
---------------
There are a bunch of post processing scripts under scripts/validaiton and scripts/python


Contact
-------

Please e-mail [Chris Minar](mailto:minarc@oregonstate.edu) if you have any
questions, suggestions or feedback.
