/***************************************************************************//**
 * \file  oscCylinder.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \based on code by Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class oscCylinder.
 */

#include <solvers/NavierStokes/oscCylinder/kernels/structure.h>
#include "oscCylinder.h"
#include <sys/stat.h>

/**
 * \brief Constructor. Copies the database and information about the computational grid.
 *
 * \param pDB database that contains all the simulation parameters
 * \param dInfo information related to the computational grid
 */
oscCylinder::oscCylinder(parameterDB *pDB, domain *dInfo)
{
	paramDB = pDB;
	domInfo = dInfo;
}

/**
 * \brief Writes data into files.
 */
void oscCylinder::writeData()
{
	parameterDB  &db = *NavierStokesSolver::paramDB;
	double dt  = db["simulation"]["dt"].get<double>();
	NavierStokesSolver::logger.startTimer("output");
	writeCommon();
	if (NavierStokesSolver::timeStep == 0)
		forceFile<<"timestep\tFx\tFxX\tFxY\tFxU\tFy\n";
	forceFile << timeStep*dt << '\t' << B.forceX[0] << '\t'<<fxx<<"\t"<<fxy<<"\t"<<fxu<<"\t" << B.forceY[0] << std::endl;
	logger.stopTimer("output");
}

/**
 * \brief Writes numerical solution at current time-step,
 *        as well as the number of iterations performed in each solver,
 *        the force,
 *        and the middle position of the body (calculated as an average of all the nodes)
 */
void oscCylinder::writeCommon()
{
	fadlunModified::writeCommon();
	midPositionFile << timeStep << '\t' << B.midX << '\t' << B.midY <<std::endl;
}

/*
 * Calculates new cell indices
 * Calculates new body bounding boxes
 * Tags Points
 * Remakes LHS matricies
 * updates Preconditioners
 */
void oscCylinder::updateSolver()
{
	fadlunModified::B.calculateCellIndices(*NavierStokesSolver::domInfo);
	fadlunModified::B.calculateBoundingBoxes(*NavierStokesSolver::paramDB, *NavierStokesSolver::domInfo);
	fadlunModified::tagPoints();
	fadlunModified::generateLHS1();//is this needed?
	fadlunModified::generateLHS2();

	NavierStokesSolver::logger.startTimer("Preconditioner");
	if (NavierStokesSolver::iterationCount2 > 100)
	{
		NavierStokesSolver::PC.update(NavierStokesSolver::LHS1, NavierStokesSolver::LHS2);
	}
	NavierStokesSolver::logger.stopTimer("Preconditioner");
}

/*
 * Calculates Force
 * Moves body
 */
void oscCylinder::moveBody()
{
	parameterDB  &db = *NavierStokesSolver::paramDB;
	fadlunModified::calculateForce();

	double *x_r	= thrust::raw_pointer_cast( &(fadlunModified::B.x[0]) ),
		   *uB_r= thrust::raw_pointer_cast( &(fadlunModified::B.uB[0]) );
	double	dt	= db["simulation"]["dt"].get<double>(),
			nu	= db["flow"]["nu"].get<double>(),
			t = dt*NavierStokesSolver::timeStep,
			D = 1.0,
			uMax = 100*nu/D, //Re
			f = uMax*D/5.0, //KC
			A = uMax/(M_PI*2.0*f), //umax
			totalPoints=fadlunModified::B.totalPoints,
			xold= fadlunModified::B.midX,
			unew,
			xnew;

	//calc new velocity and position
	xnew = A*cos(2*M_PI*f*t);
	unew = -A*2*M_PI*f*sin(2*M_PI*f*t);
	fadlunModified::B.centerVelocityU = unew;
	fadlunModified::B.midX = xnew;

	const int blocksize = 256;
	dim3 grid( int( (totalPoints)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	kernels::update_body_viv<<<grid,block>>>(x_r, uB_r, xnew-xold, unew, totalPoints);

}

/*
 * initialise the simulation
 */
void oscCylinder::initialise()
{
	fadlunModified::initialise();

	//output
	parameterDB  &db = *NavierStokesSolver::paramDB;
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();
	std::stringstream outPosition;
	outPosition << folder <<"/midPosition";
	midPositionFile.open(outPosition.str().c_str());

	double *x_r	= thrust::raw_pointer_cast( &(fadlunModified::B.x[0]) ),
		   *uB_r= thrust::raw_pointer_cast( &(fadlunModified::B.uB[0]) );
	double	dt	= db["simulation"]["dt"].get<double>(),
			nu	= db["flow"]["nu"].get<double>(),
			t = dt*NavierStokesSolver::timeStep,
			D = 1.0,
			uMax = 100*nu/D, //Re
			f = uMax*D/5.0, //KC
			A = uMax/(M_PI*2.0*f), //umax
			totalPoints=fadlunModified::B.totalPoints,
			xold= fadlunModified::B.midX,
			unew,
			xnew;

	//calc new velocity and position
	xnew = A*cos(2*M_PI*f*t);
	unew = -A*2*M_PI*f*sin(2*M_PI*f*t);
	fadlunModified::B.centerVelocityU = unew;
	fadlunModified::B.midX = xnew;

	const int blocksize = 256;
	dim3 grid( int( (totalPoints)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	kernels::update_body_viv<<<grid,block>>>(x_r, uB_r, xnew-xold, unew, totalPoints);
}

/**
 * \brief Calculates the variables at the next time step.
 */
void oscCylinder::stepTime()
{
	fadlunModified::generateRHS1();
	//arrayprint(tags,"tags","x");
	//arrayprint(tagsP,"tagsP","p");
	//arrayprint(distance_from_u_to_body,"dub","p");
	//arrayprint(distance_from_v_to_body,"dvb","p");
	NavierStokesSolver::solveIntermediateVelocity();
	//arrayprint(uhat,"uhat","x");

	fadlunModified::generateRHS2();
	//arrayprint(rhs2,"rhs2","p");
	NavierStokesSolver::solvePoisson();

	fadlunModified::velocityProjection();
	//arrayprint(u,"u","x");

	//Release the body after a certain timestep
	if (NavierStokesSolver::timeStep >= (*NavierStokesSolver::paramDB)["simulation"]["startStep"].get<int>())
	{
		moveBody();
		updateSolver();
	}
	std::cout<<timeStep<<std::endl;
	NavierStokesSolver::timeStep++;
	if (timeStep%(*paramDB)["simulation"]["nsave"].get<int>() == 0)
	{
		//arrayprint(u,"u","x");
	}
}

/**
 * \brief Prints timing information and closes the different files.
 */
void oscCylinder::shutDown()
{
	fadlunModified::shutDown();
	midPositionFile.close();
}
