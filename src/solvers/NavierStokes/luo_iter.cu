/***************************************************************************//**
 * \file  luo_iter.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \based on code by Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class luo_iter.
 */

#include <solvers/NavierStokes/luo_base/kernels/structure.h>
#include "luo_iter.h"
#include <sys/stat.h>

/**
 * \brief Constructor. Copies the database and information about the computational grid.
 *
 * \param pDB database that contains all the simulation parameters
 * \param dInfo information related to the computational grid
 */
luo_iter::luo_iter(parameterDB *pDB, domain *dInfo)
{
	paramDB = pDB;
	domInfo = dInfo;
}

void luo_iter::writeData()
{
	luo_base::writeData();
}

void luo_iter::writeCommon()
{
	luo_base::writeCommon();
}


/*
 * Calculates new cell indices
 * Calculates new body bounding boxes
 * Tags Points
 * Remakes LHS matricies
 * updates Preconditioners
 */
void luo_iter::updateSolver()
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
		PC.update(LHS1, LHS2);
	}
	logger.stopTimer("Preconditioner");
}

/*
 * Calculates Force
 * Moves body
 */
//flag this could probably be done with B.update
void luo_iter::moveBody()
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

/*
 * initialise the simulation
 */
void luo_iter::initialise()
{
	luo_base::initialise();
	luo_iter::cast();

	//output
	parameterDB  &db = *paramDB;
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();
	std::stringstream outPosition;
	outPosition << folder <<"/midPosition";
	midPositionFile.open(outPosition.str().c_str());

	double	dt	= db["simulation"]["dt"].get<double>(),
			nu	= db["flow"]["nu"].get<double>(),
			t = dt*timeStep,
			//D = 0.2,
			//uMax = 1,
			f = B.frequency,
			xCoeff = B.xCoeff,
			uCoeff = B.uCoeff,
			xPhase = B.xPhase,
			uPhase = B.uPhase,
			//KC = uMax/f/D,
			//Re = uMax*D/nu,
			totalPoints=B.totalPoints,
			xold= B.midX,
			unew,
			xnew;

	//xnew = -1/(2*M_PI)*sin(2*M_PI*f*t);
	//unew = -f*cos(2*M_PI*f*t);
	xnew = xCoeff*sin(2*M_PI*f*t + xPhase);
	unew = uCoeff*cos(2*M_PI*f*t + uPhase);

	B.centerVelocityU0 = unew;
	B.midX0 = xnew;
	B.centerVelocityU = unew;
	B.midX = xnew;

	const int blocksize = 256;
	dim3 grid( int( (totalPoints)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	B.uBk = B.uB;
	//update position/velocity for current values
	kernels::update_body_viv<<<grid,block>>>(B.x_r, B.uB_r, xnew-xold, unew, totalPoints);
	//set position/velocity for old values
	kernels::initialise_old<<<grid,block>>>(B.uBk_r,unew,totalPoints);//flag not sure if this should be done or not, as it is it simulates the body being in motion before we actually start, and it is technically more like an impulsivly started motion
																	//it effects du/dt for the calcualtion of the material derivative in the bilinear interp functions, its overall effect is pretty minimal
}

/**
 * \brief Calculates the variables at the next time step.
 */
void luo_iter::stepTime()
{
	int val = 1;
	for (int i=0; i<val; i++)
	{
	generateRHS1();
	solveIntermediateVelocity();
	//weightUhat();

	generateRHS2();
	solvePoisson();
	//weightPressure();

	if (i == val-1)
		uold = u;
	velocityProjection();
	}
	//Release the body after a certain timestep
	if (timeStep >= (*paramDB)["simulation"]["startStep"].get<int>())
	{
		moveBody();
		updateSolver();
		CFL();
	}

	timeStep++;
	std::cout<<timeStep<<std::endl;
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
void luo_iter::shutDown()
{
	luo_iter::shutDown();
	midPositionFile.close();
}
