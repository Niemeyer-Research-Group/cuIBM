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
	luoIBM::writeData();
	/*parameterDB  &db = *paramDB;
	double dt  = db["simulation"]["dt"].get<double>();
	logger.startTimer("output");
	writeCommon();
	if (timeStep == 0)
<<<<<<< HEAD
		forceFile<<"timestep\tFx\tFy\n";
	//forceFile << timeStep*dt << '\t' << B.forceX << '\t' << B.forceY << std::endl;
	forceFile	<< timeStep*dt << "\t"
				<< B.forceX << "\t"
				<< fxx << "\t"
				<< fxy << "\t"
				<< fxu << "\t"
				<< B.forceY << std::endl;
	logger.stopTimer("output");
=======
		forceFile<<"timestep\tFx\tFxX\tFxY\tFxU\tFy\n";
	forceFile << timeStep*dt << '\t' << B.forceX[0] << '\t'<<fxx<<"\t"<<fxy<<"\t"<<fxu<<"\t" << B.forceY[0] << std::endl;
	logger.stopTimer("output");*/
>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder
}

/**
 * \brief Writes numerical solution at current time-step,
 *        as well as the number of iterations performed in each solver,
 *        the force,
 *        and the middle position of the body (calculated as an average of all the nodes)
 */
void oscCylinder::writeCommon()
{
<<<<<<< HEAD
	luoIBM::writeCommon();
	if (timeStep == 1)
	{
		midPositionFile << "timeStep"
						<< std::setw(16) << "X"
						<< std::setw(16) << "Y"
						<< std::setw(16) << "U"
						<< std::setw(16) << "V" << std::endl; //flag this is formatted quite properly
	}
	midPositionFile << timeStep << "\t"
					<< std::setw(16) << std::setprecision(4) << std::fixed << B.midX << "\t"
					<< std::setw(16) << std::setprecision(4) << std::fixed << B.midY << "\t"
					<< std::setw(16) << std::setprecision(4) << std::fixed << B.centerVelocityU << "\t"
					<< std::setw(16) << std::setprecision(4) << std::fixed << B.centerVelocityV <<std::endl;
=======
	/*luoIBM::writeCommon();
	midPositionFile << timeStep << '\t' << B.midX << '\t' << B.midY <<std::endl;*/
>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder
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
<<<<<<< HEAD
	logger.startTimer("Bounding Boxes");
	B.calculateBoundingBoxes(*paramDB, *domInfo);//flag this isn't really needed because the body never moves out of the bounding box
	logger.stopTimer("Bounding Boxes");

	tagPoints();
	generateLHS1();//is this needed?
	//generateLHS2();
=======
	/*B.calculateCellIndices(*domInfo);
	B.calculateBoundingBoxes(*paramDB, *domInfo);
	tagPoints();
	generateLHS1();//is this needed?
	generateLHS2();

	logger.startTimer("Preconditioner");
	if (iterationCount2 > 100)
	{
		PC.update(LHS1, LHS2);
	}
	logger.stopTimer("Preconditioner");*/
>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder
}

/*
 * Calculates Force
 * Moves body
 */
//flag this could probably be done with B.update
void oscCylinder::moveBody()
{
<<<<<<< HEAD
	parameterDB  &db = *paramDB;
	calculateForce();
	//luoForce();
=======
	/*parameterDB  &db = *paramDB;
	luoIBM::calculateForce();
>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder

	logger.startTimer("moveBody");
	double *x_r	= thrust::raw_pointer_cast( &(B.x[0]) ),
		   *uB_r= thrust::raw_pointer_cast( &(B.uB[0]) );
	double	dt	= db["simulation"]["dt"].get<double>(),
			nu	= db["flow"]["nu"].get<double>(),
			t = dt*timeStep,
<<<<<<< HEAD
			f = B.frequency,
			xCoeff = B.xCoeff,
			uCoeff = B.uCoeff,
			xPhase = B.xPhase,
			uPhase = B.uPhase,
=======
			D = 1.0,
			uMax = 100*nu/D, //Re
			f = uMax*D/5.0, //KC
			A = uMax/(M_PI*2.0*f), //umax
>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder
			totalPoints=B.totalPoints,
			xold= B.midX,
			unew,
			xnew;

<<<<<<< HEAD
	//xnew = -1/(2*M_PI)*sin(2*M_PI*f*t);
	//unew = -f*cos(2*M_PI*f*t);
	xnew = xCoeff*sin(2*M_PI*f*t + xPhase);
	unew = uCoeff*cos(2*M_PI*f*t + uPhase);

=======
	//calc new velocity and position
	xnew = A*cos(2*M_PI*f*t);
	unew = -A*2*M_PI*f*sin(2*M_PI*f*t);
>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder
	B.centerVelocityU = unew;
	B.midX = xnew;

	const int blocksize = 256;
	dim3 grid( int( (totalPoints)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	B.uBk = B.uB;
	kernels::update_body_viv<<<grid,block>>>(x_r, uB_r, xnew-xold, unew, totalPoints);
<<<<<<< HEAD
	logger.stopTimer("moveBody");
=======
	 */
>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder
}

/*
 * initialise the simulation
 */
void oscCylinder::initialise()
{
<<<<<<< HEAD
	luoIBM::initialise();
	cfl.resize(domInfo->nx*domInfo->ny);
	distance.resize((domInfo->nx-1)*domInfo->ny + (domInfo->ny-1)*domInfo->nx);
	cfl_max = 0;
=======
	/*luoIBM::initialise();
>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder

	//output
	parameterDB  &db = *paramDB;
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();
	std::stringstream outPosition;
	outPosition << folder <<"/midPosition";
	midPositionFile.open(outPosition.str().c_str());

	double *x_r	= thrust::raw_pointer_cast( &(B.x[0]) ),
		   *uB_r= thrust::raw_pointer_cast( &(B.uB[0]) ),
		   *uB0_r=thrust::raw_pointer_cast( &(B.uBk[0]) );
	double	dt	= db["simulation"]["dt"].get<double>(),
			nu	= db["flow"]["nu"].get<double>(),
			t = dt*timeStep,
<<<<<<< HEAD
			//D = 0.2,
			//uMax = 1,
			f = B.frequency,
			xCoeff = B.xCoeff,
			uCoeff = B.uCoeff,
			xPhase = B.xPhase,
			uPhase = B.uPhase,
			//KC = uMax/f/D,
			//Re = uMax*D/nu,
=======
			D = 1.0,
			uMax = 100*nu/D, //Re
			f = uMax*D/5.0, //KC
			A = uMax/(M_PI*2.0*f), //umax
>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder
			totalPoints=B.totalPoints,
			xold= B.midX,
			unew,
			xnew;

<<<<<<< HEAD
	//xnew = -1/(2*M_PI)*sin(2*M_PI*f*t);
	//unew = -f*cos(2*M_PI*f*t);
	xnew = xCoeff*sin(2*M_PI*f*t + xPhase);
	unew = uCoeff*cos(2*M_PI*f*t + uPhase);

	B.centerVelocityU0 = unew;
	B.midX0 = xnew;
=======
	//calc new velocity and position
	xnew = A*cos(2*M_PI*f*t);
	unew = -A*2*M_PI*f*sin(2*M_PI*f*t);
>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder
	B.centerVelocityU = unew;
	B.midX = xnew;

	const int blocksize = 256;
	dim3 grid( int( (totalPoints)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
<<<<<<< HEAD
	B.uBk = B.uB;
	//update position/velocity for current values
	kernels::update_body_viv<<<grid,block>>>(x_r, uB_r, xnew-xold, unew, totalPoints);
	//set position/velocity for old values
	kernels::initialise_old<<<grid,block>>>(uB0_r,unew,totalPoints);//flag not sure if this should be done or not, as it is it simulates the body being in motion before we actually start, and it is technically more like an impulsivly started motion
																	//it effects du/dt for the calcualtion of the material derivative in the bilinear interp functions, its overall effect is pretty minimal
=======
	kernels::update_body_viv<<<grid,block>>>(x_r, uB_r, xnew-xold, unew, totalPoints);*/
>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder
}

/**
 * \brief Calculates the variables at the next time step.
 */
void oscCylinder::stepTime()
{
	luoIBM::stepTime();
	/*generateRHS1();
	//arrayprint(tags,"tags","x");
	//arrayprint(tagsP,"tagsP","p");
	//arrayprint(distance_from_u_to_body,"dub","p");
	//arrayprint(distance_from_v_to_body,"dvb","p");
	solveIntermediateVelocity();
	//arrayprint(uhat,"uhat","x");

	preRHS2();
	sizeLHS2();
	generateLHS2();
	generateRHS2();
<<<<<<< HEAD
	logger.startTimer("sort LHS2");
	LHS2.sort_by_row_and_column();
	logger.stopTimer("sort LHS2");
	//print(LHS2);
	//printLHS();

=======
	//arrayprint(rhs2,"rhs2","p");
>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder
	solvePoisson();

	interpPGN();
	velocityProjection();
<<<<<<< HEAD

	std::cout<<timeStep<<std::endl;
	timeStep++;
	if (timeStep == 1000)
	{
		arrayprint(uhat,"uhat","x",-1);
		arrayprint(pressure,"p","p",-1);
		arrayprint(u,"u","x",-1);
		arrayprint(ghostTagsP,"ghostp","p",-1);
		arrayprint(ghostTagsUV,"ghostu","x",-1);
		arrayprint(hybridTagsP,"hybridp","p",-1);
	}
=======
	//arrayprint(u,"u","x");

>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder
	//Release the body after a certain timestep
	if (timeStep >= (*paramDB)["simulation"]["startStep"].get<int>())
	{
		moveBody();
		updateSolver();
		CFL();
	}
<<<<<<< HEAD

	//print at end
	if (timeStep%(*paramDB)["simulation"]["nt"].get<int>() == 0)
	{
		//arrayprint(pressure,"p","p");
		//arrayprint(u,"u","x",-1);
		std::cout<<"Maximun CFL: " << cfl_max << std::endl;
		std::cout<<"Expected CFL: " << (*paramDB)["simulation"]["dt"].get<double>()*bc[XMINUS][0]/domInfo->mid_h << std::endl;
		std::cout<<"CFL I: " << cfl_I << std::endl;
		std::cout<<"CFL J: " << cfl_J << std::endl;
		std::cout<<"CFL ts: " << cfl_ts << std::endl;
	}
=======
	std::cout<<timeStep<<std::endl;
	timeStep++;
	if (timeStep%(*paramDB)["simulation"]["nsave"].get<int>() == 0)
	{
		//arrayprint(u,"u","x");
	}*/
>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder
}

/**
 * \brief Prints timing information and closes the different files.
 */
void oscCylinder::shutDown()
{
	luoIBM::shutDown();
	//midPositionFile.close();
}
<<<<<<< HEAD

#include "oscCylinder/intermediateVelocity.inl"
#include "oscCylinder/CFL.inl"
=======
>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder
