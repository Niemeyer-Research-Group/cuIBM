/***************************************************************************//**
 * \file  fadlunModified.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \based on code by Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class oscCylinder.
 */

#include "fadlunModified.h"
#include <sys/stat.h>

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
	NavierStokesSolver::initialiseNoBody();
	NavierStokesSolver::logger.startTimer("initialise");

	int nx = NavierStokesSolver::domInfo->nx,
		ny = NavierStokesSolver::domInfo->ny;

	int numUV = (nx-1)*ny + nx*(ny-1);
	int numP  = nx*ny;
	////////////////////////////////////////////////////////////////////////////////////////////////
	//ARRAYS
	////////////////////////////////////////////////////////////////////////////////////////////////
	//tagpoints, size uv, device
	tags.resize(numUV);//used in lhs1
	tagsOld.resize(numUV);
	tagsPOld.resize(numP);
	tags2.resize(numUV);//used in lhs1
	tagsIn.resize(numUV);//used in lhs1
	distance_from_intersection_to_node.resize(numUV);
	distance_between_nodes_at_IB.resize(numUV);
	uv.resize(numUV);

	//tagpoints, size np, device
	tagsP.resize(numP);//flag
	tagsPOut.resize(numP);//flag
	distance_from_u_to_body.resize(numP);
	distance_from_v_to_body.resize(numP);
	test.resize(numP); //flag

	cusp::blas::fill(tagsOld,-1);
	cusp::blas::fill(tagsPOld,-1);

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

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//LHS
	/////////////////////////////////////////////////////////////////////////////////////////////////
	initialiseLHS();

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
 * Initliase the LHS matrices
 */
void fadlunModified::initialiseLHS()
{
	parameterDB  &db = *NavierStokesSolver::paramDB;
	generateLHS1();
	generateLHS2();

	NavierStokesSolver::PC.generate(NavierStokesSolver::LHS1,NavierStokesSolver::LHS2, db["velocitySolve"]["preconditioner"].get<preconditionerType>(), db["PoissonSolve"]["preconditioner"].get<preconditionerType>());
	std::cout << "Assembled FADLUN LHS matrices!" << std::endl;
}


/**
 * \brief Writes data into files.
 */
void fadlunModified::writeData()
{

	parameterDB  &db = *NavierStokesSolver::paramDB;
	double dt  = db["simulation"]["dt"].get<double>();

	logger.startTimer("output");
	writeCommon();
	logger.stopTimer("output");

	calculateForce();

	logger.startTimer("output");
	if (NavierStokesSolver::timeStep == 0)
		forceFile<<"timestep\tFx\tFxX\tFxY\tFxU\tFy\n";
	forceFile << timeStep*dt << '\t' << B.forceX << '\t'<<fxx<<"\t"<<fxy<<"\t"<<fxu<<"\t" << B.forceY << std::endl;
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

void fadlunModified::stepTime()
{
	generateRHS1();
	NavierStokesSolver::solveIntermediateVelocity();

	generateRHS2();
	NavierStokesSolver::solvePoisson();

	velocityProjection();

	//std::cout<<timeStep<<std::endl;
	NavierStokesSolver::timeStep++;
}

/**
 * \brief Prints timing information and closes the different files.
 */
void fadlunModified::shutDown()
{
	NavierStokesSolver::shutDown();
	forceFile.close();
}

#include "FadlunModified/intermediateVelocity.inl"
#include "FadlunModified/intermediatePressure.inl"
#include "FadlunModified/projectVelocity.inl"
#include "FadlunModified/tagpoints.inl"
#include "FadlunModified/calculateForce.inl"
#include "FadlunModified/checkTags.inl"
