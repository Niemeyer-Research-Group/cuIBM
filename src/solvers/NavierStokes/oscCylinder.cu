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
	parameterDB  &db = *paramDB;
	double dt  = db["simulation"]["dt"].get<double>();
	logger.startTimer("output");
	writeCommon();
	if (timeStep == 0)
		forceFile<<"timestep\tFx\tFxX\tFxY\tFxU\tFy\n";
	forceFile << timeStep*dt << '\t' << B.forceX << '\t'<<fxx<<"\t"<<fxy<<"\t"<<fxu<<"\t" << B.forceY << std::endl;//flag writing from b.forcex takes forever
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
	luoIBM::writeCommon();
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
	B.calculateCellIndices(*domInfo);
	B.calculateBoundingBoxes(*paramDB, *domInfo);
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
void oscCylinder::moveBody()
{
	parameterDB  &db = *paramDB;
	calculateForce();

	double *x_r	= thrust::raw_pointer_cast( &(B.x[0]) ),
		   *uB_r= thrust::raw_pointer_cast( &(B.uB[0]) );
	double	dt	= db["simulation"]["dt"].get<double>(),
			nu	= db["flow"]["nu"].get<double>(),
			t = dt*timeStep,
			f = 1,
			totalPoints=B.totalPoints,
			xold= B.midX,
			unew,
			xnew;

	xnew = -1/(2*M_PI)*sin(2*M_PI*f*t);
	unew = -f*cos(2*M_PI*f*t);

	B.centerVelocityU = unew;
	B.midX = xnew;

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
	luoIBM::initialise();

	//output
	parameterDB  &db = *paramDB;
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();
	std::stringstream outPosition;
	outPosition << folder <<"/midPosition";
	midPositionFile.open(outPosition.str().c_str());

	double *x_r	= thrust::raw_pointer_cast( &(B.x[0]) ),
		   *uB_r= thrust::raw_pointer_cast( &(B.uB[0]) );
	double	dt	= db["simulation"]["dt"].get<double>(),
			nu	= db["flow"]["nu"].get<double>(),
			t = dt*timeStep,
			D = 0.2,
			uMax = 1,
			f = 1,
			KC = uMax/f/D,
			Re = uMax*D/nu,
			totalPoints=B.totalPoints,
			xold= B.midX,
			unew,
			xnew;

	std::cout<<"dt\tf\tuMax\tD\tnu\tRe\tKC\n";
	std::cout<<dt<<"\t"<<f<<"\t"<<uMax<<"\t"<<D<<"\t"<<nu<<"\t"<<Re<<"\t"<<KC<<"\n";
	xnew = -1/(2*M_PI)*sin(2*M_PI*f*t);
	unew = -f*cos(2*M_PI*f*t);
	B.centerVelocityU = unew;
	B.midX = xnew;

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
	generateRHS1();
	solveIntermediateVelocity();
	weightUhat();

	generateRHS2();
	solvePoisson();
	weightPressure();

	velocityProjection();
	//Release the body after a certain timestep
	if (timeStep >= (*paramDB)["simulation"]["startStep"].get<int>())
	{
		moveBody();
		updateSolver();
	}

	std::cout<<timeStep<<std::endl;
	timeStep++;

	if (timeStep%(*paramDB)["simulation"]["nsave"].get<int>() == 0)
	{
		//arrayprint(pressure,"p","p");
		//arrayprint(u,"u","x");
	}
}

/**
 * \brief Prints timing information and closes the different files.
 */
void oscCylinder::shutDown()
{
	luoIBM::shutDown();
	midPositionFile.close();
}

#include "oscCylinder/intermediateVelocity.inl"
