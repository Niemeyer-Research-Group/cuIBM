# file: Makefile
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Compiles and links cuIBM code.


# compiler: nvcc is NVIDIA's CUDA compiler
CC = /usr/local/cuda-7.5/bin/nvcc $(OSXOPTS)

# compiler options
# -O3: optimization flag
#-g nvcc debug  makes execution take 20x longer
#-lineinfo nvcc debug
#-Xcompiler -Q passes something to gcc c++ compiler
CCFLAGS = -arch=compute_20 -O3# -lineinfo

# variables
RM = rm
MAKE = make

# root directory of the project
# return the absolute path of the directory  where is located the Makefile
# variable MAKEFILE_LIST lists all Makefiles in the working directory
# `lastword` picks the last element of the list
PROJ_ROOT = /scratch/src/cuIBM
#PROJ_ROOT = $(abspath $(dir $(lastword $(MAKEFILE_LIST))))

#cusp directory
CUSP_ROOT = /scratch/src/lib/cusplibrary-0.5.1 #hardcode path

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
	cd /scratch/src/cuIBM/external; $(MAKE) $(MFLAGS) all

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
test:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/cylinder/test
############################
#OSC CLYINDER IN NO FLOW
############################
oscStatic:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/static

############################
#OSC CYLINDER WITH FLOW
############################
oscflowa:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/flow/a

oscflowb:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/flow/b
	
oscflowc:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/flow/c
	
oscflowd:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/flow/d
	
oscflowe:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/flow/e
	
oscflowf:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/flow/f
	
oscflowg:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/flow/g
	
oscflowh:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/flow/h
	
oscflowi:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/flow/i
	
oscflowj:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/flow/j
	
oscflowk:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/flow/k
	
oscflowl:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/flow/l
	
oscflowm:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/flow/m
	
oscflown:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/flow/n

oscflowo:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/flow/o
	
oscflowp:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/flow/p
############################
#VIV
############################
vivUred3:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/VIV/Ured3
	
vivUred4:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/VIV/Ured4
	
vivUred5:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/VIV/Ured5
	
vivUred6:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/VIV/Ured6
	
vivUred7:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/VIV/Ured7
	
vivUred8:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/osc/VIV/Ured8

############################
#CYLINDER
############################
cylinderRe40:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/cylinder/Re40
	
cylinderRe75:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/cylinder/Re75	
	
cylinderRe100:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/cylinder/Re100	
	
cylinderRe150:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/cylinder/Re150	
	
cylinderRe200:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/cylinder/Re200
	
cylinderRe300:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/cylinder/Re300

cylinderRe550:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/cylinder/Re550
	
cylinderRe3000:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/cylinder/Re3000
	
############################
#LID DRIVEN CAVITY
############################

lidDrivenCavityRe100:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/lidDrivenCavity/Re100
	
lidDrivenCavityRe1000:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/lidDrivenCavity/Re1000
	
lidDrivenCavityRe5000:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/lidDrivenCavity/Re5000
	
lidDrivenCavityRe10000:
	$(PROJ_ROOT)/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/lidDrivenCavity/Re10000
