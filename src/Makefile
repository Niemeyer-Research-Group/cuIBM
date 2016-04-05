# file: Makefile
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Compiles and links cuIBM code.


# compiler: nvcc is NVIDIA's CUDA compiler
CC = /usr/local/apps/cuda/cuda-7.0.28/bin/nvcc $(OSXOPTS)

# compiler options
# -arch=compute_20: compile to use double precision on GPU use _12 -03 for your shitty computer
# -O3: optimization flag
CCFLAGS = -arch=compute_20 -O0#-O3# -Xcompiler -Q
#CCFLAGS = -arch=compute_12 -O3

# variables
RM = rm
MAKE = make

# root directory of the project
# return the absolute path of the directory  where is located the Makefile
# variable MAKEFILE_LIST lists all Makefiles in the working directory
# `lastword` picks the last element of the list
PROJ_ROOT = /scratch/src/cuIBM-FSI

#cusp directory
CUSP_ROOT = /scratch/src/lib/cusplibrary-0.4.0 #hardcode path

# source code directory
SRC_DIR = $(PROJ_ROOT)/src

# directory where object files are stored
BUILD_DIR = $(PROJ_ROOT)/build

# directory where binary executables are stored
BIN_DIR = $(PROJ_ROOT)/bin

# extension of source files
SRC_EXT = .cu

# cuIBM executable
TARGET = $(PROJ_ROOT)/bin/cuIBM

# list all source files in the source directory
SRCS = $(shell find $(SRC_DIR) -type f -name \*$(SRC_EXT))
# absolute path of all object files to be created
OBJS = $(patsubst $(SRC_DIR)/%, $(BUILD_DIR)/%, $(SRCS:$(SRC_EXT)=.o))

# include header files from cuIBM and CUSP library
INC = -I $(SRC_DIR) -I $(CUSP_ROOT)

# path of the YAML static library
EXT_LIBS = $(PROJ_ROOT)/external/lib/libyaml-cpp.a
# include YAML header files
INC += -I $(PROJ_ROOT)/external/yaml-cpp/include


.PHONY: all

all: $(TARGET)

$(TARGET): $(OBJS) $(EXT_LIBS)
	@echo "\nLinking ..."
	@mkdir -p $(BIN_DIR)
	$(CC) $^ -o $@

$(EXT_LIBS):
	@echo "\nCreating static library $@ ..."
	cd /scratch/src/cuIBM-FSI/external; $(MAKE) $(MFLAGS) all

$(BUILD_DIR)/%.o: $(SRC_DIR)/%$(SRC_EXT)
	@mkdir -p $(@D)
	$(CC) $(CCFLAGS) $(INC) -c $< -o $@

################################################################################

.PHONY: doc

DOC_DIR = $(PROJ_ROOT)/doc
DOXYGEN = doxygen

doc:
	@echo "\nGenerating Doxygen documentation ..."
	cd $(DOC_DIR); $(DOXYGEN) Doxyfile

################################################################################

.PHONY: clean cleanexternal cleandoc cleanall

clean:
	@echo "\nCleaning cuIBM ..."
	$(RM) -rf $(BUILD_DIR) $(BIN_DIR)

cleanexternal:
	@echo "\nCleaning external YAML ..."
	cd external; $(MAKE) $(MFLAGS) clean

cleandoc:
	@echo "\nCleaning documentation ..."
	find $(DOC_DIR) ! -name 'Doxyfile' -type f -delete
	find $(DOC_DIR)/* ! -name 'Doxyfile' -type d -delete

cleanall: clean cleanexternal cleandoc

################################################################################

# commands to run cuIBM simulations

cylinderRe40:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM-FSI/validation/cylinder/Re40
	
cylinderRe75:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM-FSI/validation/cylinder/Re75	
	
cylinderRe100:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM-FSI/validation/cylinder/Re100	
	
cylinderRe150:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM-FSI/validation/cylinder/Re150	
	
cylinderRe200:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM-FSI/validation/cylinder/Re200	

cylinderRe550:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM-FSI/validation/cylinder/Re550
	
cylinderRe3000:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM-FSI/validation/cylinder/Re3000

lidDrivenCavityRe100:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM-FSI/validation/lidDrivenCavity/Re100
	
lidDrivenCavityRe1000:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM-FSI/validation/lidDrivenCavity/Re1000
	
lidDrivenCavityRe5000:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM-FSI/validation/lidDrivenCavity/Re5000
	
lidDrivenCavityRe10000:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM-FSI/validation/lidDrivenCavity/Re10000
