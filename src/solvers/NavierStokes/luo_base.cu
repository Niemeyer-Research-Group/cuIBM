/***************************************************************************//**
 * \file  luo_base.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \based on code by Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class oscCylinder.
 */

#include "luo_base.h"
#include <sys/stat.h>
#include <solvers/NavierStokes/luo_base/kernels/structure.h>

/**
 * \brief Constructor. Copies the database and information about the computational grid.
 *
 * \param pDB database that contains all the simulation parameters
 * \param dInfo information related to the computational grid
 */
luo_base::luo_base(parameterDB *pDB, domain *dInfo)
{
	paramDB = pDB;
	domInfo = dInfo;
}

/*
 * Initialise the solver
 */
void luo_base::initialise()
{
	NavierStokesSolver::initialiseNoBody();
	logger.startTimer("initialise");

	////////////////////////////////////////////////////////////////////////////////////////////////
	//Cast Luo
	////////////////////////////////////////////////////////////////////////////////////////////////
	luo_base::cast();
	std::cout << "luo_base: resized and cast!" << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////////////
	//Initialize Bodies
	////////////////////////////////////////////////////////////////////////////////////////////////
	B.initialise((*paramDB), *domInfo);
	std::cout << "luo_base: Initialised bodies!" << std::endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//TAG POINTS
	/////////////////////////////////////////////////////////////////////////////////////////////////
	tagPoints();
	std::cout << "luo_base: Tagged points!" << std::endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//OUTPUT
	/////////////////////////////////////////////////////////////////////////////////////////////////
	parameterDB  &db = *NavierStokesSolver::paramDB;
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();
	std::stringstream out;
	out << folder << "/forces";
	forceFile.open(out.str().c_str());
	std::string folder2 = db["inputs"]["caseFolder"].get<std::string>();
	std::stringstream outPosition2;
	outPosition2 << folder2 <<"/midPosition";
	midPositionFile.open(outPosition2.str().c_str());

	logger.stopTimer("initialise");
}

/*
 * Initialise the LHS matricies
 */
void luo_base::initialiseLHS()
{
}

/**
 * \brief Writes data into files.
 */
void luo_base::writeData()
{
	logger.startTimer("output");
	writeCommon();

	if (NavierStokesSolver::timeStep == 1)
		forceFile<<"timestep\tFx\tFy\n";
	forceFile << timeStep*dt << '\t' << B.forceX << '\t'<< B.forceY << std::endl;
	logger.stopTimer("output");
}

/*
 * Calculates new cell indices
 * Calculates new body bounding boxes
 * Tags Points
 * Remakes LHS matricies
 * updates Preconditioners
 */
void luo_base::updateSolver()
{
	logger.startTimer("Bounding Boxes");
	B.calculateBoundingBoxes(*paramDB, *domInfo);//flag this isn't really needed because the body never moves out of the bounding box
	logger.stopTimer("Bounding Boxes");

	tagPoints();
	generateLHS1();//is this needed?
	generateLHS2();

	logger.startTimer("Preconditioner");
	if (iterationCount2 > 100)
	{
		PC.update1(LHS1);
		PC.update2(LHS2);
	}
	logger.stopTimer("Preconditioner");
}

/*
 * Calculates Force
 * Moves body
 */
//flag this could probably be done with B.update
void luo_base::moveBody()
{
	calculateForce();
	//luoForce();

	logger.startTimer("moveBody");
	double	t = dt*timeStep,
			f = B.frequency,
			xCoeff = B.xCoeff,
			uCoeff = B.uCoeff,
			xPhase = B.xPhase,
			uPhase = B.uPhase,
			totalPoints=B.totalPoints,
			xold= B.midX,
			unew,
			xnew;

	//xnew = -1/(2*M_PI)*sin(2*M_PI*f*t);
	//unew = -f*cos(2*M_PI*f*t);
	xnew = xCoeff*sin(2*M_PI*f*t + xPhase);
	unew = uCoeff*cos(2*M_PI*f*t + uPhase);

	B.centerVelocityU = unew;
	B.midX = xnew;

	const int blocksize = 256;
	dim3 grid( int( (totalPoints)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	B.uBk = B.uB;
	kernels::update_body_viv<<<grid,block>>>(B.x_r, B.uB_r, xnew-xold, unew, totalPoints);
	logger.stopTimer("moveBody");
}

/**
 * \brief Writes numerical solution at current time-step,
 *        as well as the number of iterations performed in each solver.
 */
void luo_base::writeCommon()
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
	if (timeStep == 1)
	{
		midPositionFile << "timeStep"
						<< "\t" << "X"
						<< "\t" << "Y"
						<< "\t" << "U"
						<< "\t" << "V" << std::endl; //flag this is formatted quite properly
	}
	midPositionFile << timeStep
					<< "\t" << B.midX
					<< "\t" << B.midY
					<< "\t" << B.centerVelocityU
					<< "\t" << B.centerVelocityV <<std::endl;
}

void luo_base::stepTime()
{
	_intermediate_velocity();

	_pressure();

	_project_velocity();

	_post_step();
}

void luo_base::_intermediate_velocity()
{}
void luo_base::_pressure()
{}
void luo_base::_post_step()
{
	//Release the body after a certain timestep
	if (timeStep >= (*paramDB)["simulation"]["startStep"].get<int>())
	{
		moveBody();
		updateSolver();
		CFL();
	}

	//arrayprint(u,"u","x",-1);
	//arrayprint(ghostTagsUV,"ghostu","x",-1);

	//update time
	timeStep++;
	std::cout<<timeStep<<std::endl;

	//print stuff if its done
	if (timeStep%(*paramDB)["simulation"]["nt"].get<int>() == 0)
	{
		std::cout<<"Maximun CFL: " << cfl_max << std::endl;
		std::cout<<"Expected CFL: " << (*paramDB)["simulation"]["dt"].get<double>()*bc[XMINUS][0]/domInfo->mid_h << std::endl;
		std::cout<<"CFL I: " << cfl_I << std::endl;
		std::cout<<"CFL J: " << cfl_J << std::endl;
		std::cout<<"CFL ts: " << cfl_ts << std::endl;
	}
}

/**
 * \brief Prints timing information and closes the different files.
 */
void luo_base::shutDown()
{
	NavierStokesSolver::shutDown();
	forceFile.close();
	midPositionFile.close();
}
