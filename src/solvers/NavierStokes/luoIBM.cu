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
	//ARRAYS
	////////////////////////////////////////////////////////////////////////////////////////////////
	pressureStar.resize(numP);
	ustar.resize(numUV);
	ghostTagsUV.resize(numUV);
	hybridTagsUV.resize(numUV);
	hybridTagsUV2.resize(numUV);
	body_intercept_x.resize(numUV);
	body_intercept_y.resize(numUV);
	image_point_x.resize(numUV);
	image_point_y.resize(numUV);
	body_intercept_p_x.resize(numP);
	body_intercept_p_y.resize(numP);
	body_intercept_p.resize(numP);
	image_point_p_x.resize(numP);
	image_point_p_y.resize(numP);
	distance_from_intersection_to_node.resize(numUV);
	distance_between_nodes_at_IB.resize(numUV);
	uv.resize(numUV);

	//testing
	x1_ip.resize(numUV);
	x2_ip.resize(numUV);
	y1_ip.resize(numUV);
	y2_ip.resize(numUV);
	x1_ip_p.resize(numP);
	x2_ip_p.resize(numP);
	y1_ip_p.resize(numP);
	y2_ip_p.resize(numP);
	image_point_u.resize(numUV);
	x1.resize(numUV);
	x2.resize(numUV);
	x3.resize(numUV);
	x4.resize(numUV);
	y1.resize(numUV);
	y2.resize(numUV);
	y3.resize(numUV);
	y4.resize(numUV);
	q1.resize(numUV);
	q2.resize(numUV);
	q3.resize(numUV);
	q4.resize(numUV);
	x1_p.resize(numP);
	x2_p.resize(numP);
	x3_p.resize(numP);
	x4_p.resize(numP);
	y1_p.resize(numP);
	y2_p.resize(numP);
	y3_p.resize(numP);
	y4_p.resize(numP);
	q1_p.resize(numP);
	q2_p.resize(numP);
	q3_p.resize(numP);
	q4_p.resize(numP);
	a0.resize(numP);
	a1.resize(numP);
	a2.resize(numP);
	a3.resize(numP);
	dudt.resize(numP);
	ududx.resize(numP);
	vdudy.resize(numP);
	dvdt.resize(numP);
	udvdx.resize(numP);
	vdvdy.resize(numP);

	//tagpoints, size nump
	ghostTagsP.resize(numP);
	hybridTagsP.resize(numP);
	distance_from_u_to_body.resize(numP);
	distance_from_v_to_body.resize(numP);

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
		forceFile<<"timestep\told\tPressure\tdudn\tnew\n";
	forceFile << timeStep*dt << '\t' << B.forceX << '\t'<<fxx<<"\t"<<fxy<<"\t"<<fxu<<"\t" << B.forceY << std::endl;
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

/**
 * \brief Calculates the variables at the next time step.
 */
void luoIBM::stepTime()
{
	generateRHS1();
	solveIntermediateVelocity();
	weightUhat();

	generateRHS2();
	solvePoisson();
	weightPressure();

	velocityProjection();

	timeStep++;
	std::cout<<timeStep<<std::endl;
}

/**
 * \brief Prints timing information and closes the different files.
 */
void luoIBM::shutDown()
{
	NavierStokesSolver::shutDown();
	forceFile.close();
}
