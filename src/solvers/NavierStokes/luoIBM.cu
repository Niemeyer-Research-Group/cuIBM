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
	NavierStokesSolver::logger.startTimer("initialise");
	int nx = NavierStokesSolver::domInfo->nx,
		ny = NavierStokesSolver::domInfo->ny;

	int numUV = (nx-1)*ny + nx*(ny-1);
	int numP  = nx*ny;
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

	//interp
	detA.resize(numP);
	alpha.resize(numP);
	b11.resize(numP);
	b12.resize(numP);
	b13.resize(numP);
	b14.resize(numP);
	b21.resize(numP);
	b22.resize(numP);
	b23.resize(numP);
	b24.resize(numP);
	b31.resize(numP);
	b32.resize(numP);
	b33.resize(numP);
	b34.resize(numP);
	b41.resize(numP);
	b42.resize(numP);
	b43.resize(numP);
	b44.resize(numP);
	stencilCoef.resize(numP);
	interpCoef.resize(numP);
	countD.resize(numP);
	countH.resize(numP);
	q1flag.resize(numP);
	q2flag.resize(numP);
	q3flag.resize(numP);
	q4flag.resize(numP);
	index1.resize(numP);
	index2.resize(numP);
	index3.resize(numP);
	index4.resize(numP);

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
	parameterDB  &db = *NavierStokesSolver::paramDB;
	int nx = domInfo->nx,
		ny = domInfo->ny,
		numUV = (nx-1)*ny + (ny-1)*nx;
	LHS1.resize(numUV, numUV, (nx-1)*ny*5 - 2*ny-2*(nx-1)       +        (ny-1)*nx*5 - 2*(ny-1) - 2*nx);
	//LHS2.resize(nx*ny, nx*ny, 5*nx*ny - 2*ny-2*nx + nx*ny*3); //flag
	generateLHS1();
	//generateLHS2();

	//NavierStokesSolver::PC.generate(NavierStokesSolver::LHS1,NavierStokesSolver::LHS2, db["velocitySolve"]["preconditioner"].get<preconditionerType>(), db["PoissonSolve"]["preconditioner"].get<preconditionerType>());
	std::cout << "Assembled LUO LHS matrices!" << std::endl;
}

/**
 * \brief Writes data into files.
 */
void luoIBM::writeData()
{
	parameterDB  &db = *NavierStokesSolver::paramDB;
	double dt  = db["simulation"]["dt"].get<double>();

	logger.startTimer("output");
	writeCommon();
	logger.stopTimer("output");

	logger.startTimer("calculateForce");
	calculateForce();
	//luoForce();
	logger.stopTimer("calculateForce");

	logger.startTimer("output");
	if (NavierStokesSolver::timeStep == 1)
<<<<<<< HEAD
		forceFile<<"timestep\told\tPressure\tdudn\tnew\n";
	forceFile << timeStep*dt << '\t' << B.forceX << '\t'<<fxx<<"\t"<<fxy<<"\t"<<fxu<<"\t" << B.forceY << std::endl;
=======
		forceFile<<"timestep\tFx\tFxX\tFxY\tFxU\tFy\n";
	forceFile << timeStep*dt << '\t' << B.forceX[0] << '\t'<<fxx<<"\t"<<fxy<<"\t"<<fxu<<"\t" << B.forceY[0] << std::endl;

>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder
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
}

/**
 * \brief Calculates the variables at the next time step.
 */
void luoIBM::stepTime()
{
	generateRHS1();
	solveIntermediateVelocity();
	weightUhat();

	preRHS2();
	sizeLHS2();
	generateLHS2();
	generateRHS2();
	LHS2.sort_by_row_and_column();
	//print(LHS2);
	//printLHS();
	PC.generate(LHS1,LHS2, (*paramDB)["velocitySolve"]["preconditioner"].get<preconditionerType>(), (*paramDB)["PoissonSolve"]["preconditioner"].get<preconditionerType>());

	solvePoisson();

	interpPGN();
	velocityProjection();

	std::cout<<timeStep<<std::endl;
	timeStep++;
	if (timeStep == 1000)
	{
<<<<<<< HEAD
		arrayprint(uhat,"uhat","x",-1);
		arrayprint(pressure,"p","p",-1);
		arrayprint(u,"u","x",-1);
		arrayprint(ghostTagsP,"ghostp","p",-1);
		arrayprint(ghostTagsUV,"ghostu","x",-1);
		arrayprint(hybridTagsP,"hybridp","p",-1);
=======
		arrayprint(u,"u","x");
		arrayprint(u,"v","y");
		arrayprint(pressure,"pressure","p");
		divergence();
>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder
	}
}

/**
 * \brief Prints timing information and closes the different files.
 */
void luoIBM::shutDown()
{
	NavierStokesSolver::shutDown();
	forceFile.close();
}

#include "luoIBM/intermediateVelocity.inl"
#include "luoIBM/intermediatePressure.inl"
#include "luoIBM/projectVelocity.inl"
#include "luoIBM/tagpoints.inl"
#include "luoIBM/calculateForce.inl"
#include "luoIBM/testing.inl"
