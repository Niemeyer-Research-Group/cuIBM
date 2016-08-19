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
	NavierStokesSolver::initialiseNoBody();
	logger.startTimer("initialise");

	////////////////////////////////////////////////////////////////////////////////////////////////
	//Cast Luo
	////////////////////////////////////////////////////////////////////////////////////////////////
	luoIBM::cast();

	std::cout << "Luo IBM resized and cast!" << std::endl;
	////////////////////////////////////////////////////////////////////////////////////////////////
	//Initialize Bodies
	////////////////////////////////////////////////////////////////////////////////////////////////
	B.initialise((*paramDB), *domInfo);
	std::cout << "Initialised bodies!" << std::endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//TAG POINTS
	/////////////////////////////////////////////////////////////////////////////////////////////////
	tagPoints();
	std::cout << "Tagged points!" << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////////////
	//Initialize Velocity
	////////////////////////////////////////////////////////////////////////////////////////////////
	zeroVelocity();//sets the velocity inside the body to 0
	std::cout << "Inside velocity set to body velocity!" << std::endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//LHS
	/////////////////////////////////////////////////////////////////////////////////////////////////
	initialiseLHS();
	std::cout << "LHS Initialised!" << std::endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//OUTPUT
	/////////////////////////////////////////////////////////////////////////////////////////////////
	parameterDB  &db = *NavierStokesSolver::paramDB;
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();
	std::stringstream out;
	out << folder << "/forces";
	forceFile.open(out.str().c_str());

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
	std::cout << "Assembled LUO LHS matrices!" << std::endl;
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
	NavierStokesSolver::writeCommon();
	parameterDB  &db = *NavierStokesSolver::paramDB;
	int nsave = db["simulation"]["nsave"].get<int>();
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();

	// write body output
	if (timeStep % nsave == 0)
	{
		B.writeToFile(folder, NavierStokesSolver::timeStep);
	}

	// write the number of iterations for each solve
	iterationsFile << timeStep << '\t' << iterationCount1 << '\t' << iterationCount2 << std::endl;
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
	NavierStokesSolver::shutDown();
	forceFile.close();
}
