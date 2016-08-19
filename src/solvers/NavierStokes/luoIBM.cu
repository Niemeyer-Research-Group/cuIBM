/***************************************************************************//**
 * \file  luoIBM.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \based on code by Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class oscCylinder.
 */

#include "luoIBM.h"
#include <sys/stat.h>

/**
 * \brief Constructor. Copies the database and information about the computational grid.
 *
 * \param pDB database that contains all the simulation parameters
 * \param dInfo information related to the computational grid
 */
luoIBM::luoIBM(parameterDB *pDB, domain *dInfo)
{
	paramDB = pDB;
	domInfo = dInfo;
}

/*
 * Initialise the solver
 */
void luoIBM::initialise()
{
	luo_base::initialise();
	logger.startTimer("initialise");

	////////////////////////////////////////////////////////////////////////////////////////////////
	//Cast Luo
	////////////////////////////////////////////////////////////////////////////////////////////////
	luoIBM::cast();

	std::cout << "luoIBM: resized and cast!" << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////////////
	//Initialize Velocity
	////////////////////////////////////////////////////////////////////////////////////////////////
	zeroVelocity();//sets the velocity inside the body to 0
	std::cout << "luoIBM: Inside velocity set to body velocity!" << std::endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//LHS
	/////////////////////////////////////////////////////////////////////////////////////////////////
	initialiseLHS();
	std::cout << "luoIBM: LHS Initialised!" << std::endl;

	logger.stopTimer("initialise");
}

/*
 * Initialise the LHS matricies
 */
void luoIBM::initialiseLHS()
{
	generateLHS1();
	generateLHS2();

	PC.generate(LHS1,LHS2, (*paramDB)["velocitySolve"]["preconditioner"].get<preconditionerType>(), (*paramDB)["PoissonSolve"]["preconditioner"].get<preconditionerType>());
}

/**
 * \brief Writes data into files.
 */
void luoIBM::writeData()
{
	logger.startTimer("output");
	writeCommon();
	logger.stopTimer("output");

	logger.startTimer("calculateForce");
	calculateForce();
	//luoForce();
	logger.stopTimer("calculateForce");

	logger.startTimer("output");
	if (NavierStokesSolver::timeStep == 1)
		forceFile<<"timestep\tFx\tFy\n";
	forceFile << timeStep*dt << '\t' << B.forceX << '\t'<< B.forceY << std::endl;
	logger.stopTimer("output");
}

/**
 * \brief Writes numerical solution at current time-step,
 *        as well as the number of iterations performed in each solver.
 */
void luoIBM::writeCommon()
{
	luo_base::writeCommon();
}

void luoIBM::_intermediate_velocity()
{
	generateRHS1();
	solveIntermediateVelocity();
	weightUhat();
}
void luoIBM::_pressure()
{
	generateRHS2();
	solvePoisson();
	weightPressure();
}

/**
 * \brief Prints timing information and closes the different files.
 */
void luoIBM::shutDown()
{
	luo_base::shutDown();
}
