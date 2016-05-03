/***************************************************************************//**
 * \mainpage cuIBM
 *
 *		A GPU-based Immersed Boundary Method
 *
 * \author Anush Krishnan (anush@bu.edu)
 */


/***************************************************************************//**
 * \file cuIBM.cu
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Main source-file of \c cuIBM.
 */

#include "domain.h"
#include "io/io.h"
#include "solvers/NavierStokes/NavierStokesSolver.h"
#include "solvers/NavierStokes/FSI.h"

int main(int argc, char **argv)
{
	// initialize the computational domain
	domain dom_info;

	// initialize the parameters of the simulation
	parameterDB paramDB;

	// read input .yaml files
	io::readInputs(argc, argv, paramDB, dom_info);

	//print simulation info
	io::printSimulationInfo(paramDB, dom_info);

	// create and initialize the flow solver, I think this can be simplified/streamlined now that there is only one solver
	NavierStokesSolver *solver = 0;
	if (paramDB["simulation"]["FSI"].get<int>() == 0)
		solver = new NavierStokesSolver(&paramDB, &dom_info);
	else
		solver = new FSI(&paramDB, &dom_info);
	solver->initialise();

	//prints to output and files
	io::printDeviceMemoryUsage();
	io::writeInfoFile(paramDB, dom_info);

	// time-step loop
	while (!solver->finished())
	{
		solver->stepTime();
		solver->writeData();
	}

	solver->shutDown();
}
