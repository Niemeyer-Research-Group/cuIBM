/***************************************************************************//**
 * \file  fadlunModified.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \based on code by Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class oscCylinder.
 */

#include "fadlunModified.h"
#include <sys/stat.h>
#include <solvers/NavierStokes/FadlunModified/kernels/projectVelocity.h>

/**
 * \brief Constructor. Copies the database and information about the computational grid.
 *
 * \param pDB database that contains all the simulation parameters
 * \param dInfo information related to the computational grid
 */
fadlunModified::fadlunModified(parameterDB *pDB, domain *dInfo)
{
	paramDB = pDB;
	domInfo = dInfo;
}

/*
 * Initialise the solver
 */
void fadlunModified::initialise()
{
	luo_base::initialise();
	logger.startTimer("initialise");

	////////////////////////////////////////////////////////////////////////////////////////////////
	//cast
	////////////////////////////////////////////////////////////////////////////////////////////////
	fadlunModified::cast();
	std::cout << "fadlun arrays resized and cast!" << std::endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//LHS
	/////////////////////////////////////////////////////////////////////////////////////////////////
	initialiseLHS();
}

/*
 * Initliase the LHS matrices
 */
void fadlunModified::initialiseLHS()
{
	parameterDB  &db = *NavierStokesSolver::paramDB;
	int nx = domInfo->nx,
		ny = domInfo->ny,
		numUV = (nx-1)*ny + (ny-1)*nx;
	LHS1.resize(numUV, numUV, (nx-1)*ny*5 - 2*ny-2*(nx-1)       +        (ny-1)*nx*5 - 2*(ny-1) - 2*nx);
	LHS2.resize(nx*ny, nx*ny, 5*nx*ny - 2*ny-2*nx); //flag this should have some zero terms in it because no nodes are being removing to account for the different stencil at the body
	generateLHS1();
	generateLHS2();
	LHS1.sort_by_row_and_column();
	LHS2.sort_by_row_and_column();

	NavierStokesSolver::PC.generate1(NavierStokesSolver::LHS1, db["velocitySolve"]["preconditioner"].get<preconditionerType>());
	NavierStokesSolver::PC.generate2(NavierStokesSolver::LHS2, db["PoissonSolve"]["preconditioner"].get<preconditionerType>());
	std::cout << "Assembled FADLUN LHS matrices!" << std::endl;
}

/**
 * \brief Writes data into files.
 */
void fadlunModified::writeData()
{
	logger.startTimer("output");

	writeCommon();
	calculateForce();
	if (NavierStokesSolver::timeStep == 0)
		forceFile<<"timestep\tFx\tFxX\tFxY\tFxU\tFy\n";
	forceFile << timeStep*dt << '\t' << B.forceX[0] << '\t'<<fxx<<"\t"<<fxy<<"\t"<<fxu<<"\t" << B.forceY[0] << std::endl;

	logger.stopTimer("output");
}

/**
 * \brief Writes numerical solution at current time-step,
 *        as well as the number of iterations performed in each solver.
 */
void fadlunModified::writeCommon()
{
	NavierStokesSolver::writeCommon();
	parameterDB  &db = *NavierStokesSolver::paramDB;
	int nsave = db["simulation"]["nsave"].get<int>();
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();

	// write body output
	if (timeStep % nsave == 0)
	{
		B.writeToFile(folder, NavierStokesSolver::timeStep);
	}
}

void fadlunModified::_intermediate_velocity()
{
	generateRHS1();
	solveIntermediateVelocity();
}

void fadlunModified::_pressure()
{
	generateRHS2();
	solvePoisson();
}

void fadlunModified::_project_velocity()
{
	logger.startTimer("Velocity Projection");

	const int blocksize = 256;

	dim3 dimGridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 dimBlockU(blocksize, 1);
	dim3 dimGridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 dimBlockV(blocksize, 1);

	kernels::project_velocity_X<<<dimGridU,dimBlockU>>>(u_r, uhat_r, uold_r, pressure_r, ghostTagsP_r, dx_r, dt, nx, ny);
	kernels::project_velocity_Y<<<dimGridV,dimBlockV>>>(u_r, uhat_r, uold_r, pressure_r, ghostTagsP_r, dy_r, dt, nx, ny);

	logger.stopTimer("Velocity Projection");
}

/**
 * \brief Prints timing information and closes the different files.
 */
void fadlunModified::shutDown()
{
	NavierStokesSolver::shutDown();
	forceFile.close();
}
