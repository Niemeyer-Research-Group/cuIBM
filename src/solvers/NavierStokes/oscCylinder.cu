/***************************************************************************//**
 * \file  oscCylinder.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \based on code by Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class oscCylinder.
 */

#include <solvers/NavierStokes/kernels/structure.h>
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

void oscCylinder::writeData()
{
	parameterDB  &db = *NavierStokesSolver::paramDB;
	double dt  = db["simulation"]["dt"].get<double>();
	NavierStokesSolver::logger.startTimer("output");
	NavierStokesSolver::writeCommon();
	if (NavierStokesSolver::timeStep == 0)
		forceFile<<"timestep\tFx\tFxX\tFxY\tFxU\tFy\n";
	forceFile << timeStep*dt << '\t' << B.forceX[0] << '\t'<<fxx<<"\t"<<fxy<<"\t"<<fxu<<"\t" << B.forceY[0] << std::endl;
	logger.stopTimer("output");
}

void oscCylinder::updateSolver()
{
	NavierStokesSolver::B.calculateCellIndices(*NavierStokesSolver::domInfo);
	NavierStokesSolver::B.calculateBoundingBoxes(*NavierStokesSolver::paramDB, *NavierStokesSolver::domInfo);
	NavierStokesSolver::tagPoints();
	NavierStokesSolver::generateLHS1();//is this needed?
	NavierStokesSolver::generateLHS2();

	NavierStokesSolver::logger.startTimer("Preconditioner");
	if (NavierStokesSolver::iterationCount2 > 100)
	{
		NavierStokesSolver::PC.update(NavierStokesSolver::LHS1, NavierStokesSolver::LHS2);
	}
	NavierStokesSolver::logger.stopTimer("Preconditioner");
}

void oscCylinder::moveBody()
{
	parameterDB  &db = *NavierStokesSolver::paramDB;
	NavierStokesSolver::calculateForce();

	double *x_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::B.x[0]) ),
		   *uB_r= thrust::raw_pointer_cast( &(NavierStokesSolver::B.uB[0]) );
	double	dt	= db["simulation"]["dt"].get<double>(),
			nu	= db["flow"]["nu"].get<double>(),
			t = dt*NavierStokesSolver::timeStep,
			D = 1.0,
			uMax = 100*nu/D, //Re
			f = uMax*D/5.0, //KC
			A = uMax/(M_PI*2.0*f), //umax
			totalPoints=NavierStokesSolver::B.totalPoints,
			xold= NavierStokesSolver::B.midX,
			unew,
			xnew;

	//calc new velocity and position
	xnew = A*cos(2*M_PI*f*t);
	unew = -A*2*M_PI*f*sin(2*M_PI*f*t);
	NavierStokesSolver::B.centerVelocityU = unew;
	NavierStokesSolver::B.midX = xnew;

	const int blocksize = 256;
	dim3 grid( int( (totalPoints)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	kernels::update_body_viv<<<grid,block>>>(x_r, uB_r, xnew-xold, unew, totalPoints);

}

void oscCylinder::initialise()
{
	parameterDB  &db = *NavierStokesSolver::paramDB;

	double *x_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::B.x[0]) ),
		   *uB_r= thrust::raw_pointer_cast( &(NavierStokesSolver::B.uB[0]) );
	double	dt	= db["simulation"]["dt"].get<double>(),
			nu	= db["flow"]["nu"].get<double>(),
			t = dt*NavierStokesSolver::timeStep,
			D = 1.0,
			uMax = 100*nu/D, //Re
			f = uMax*D/5.0, //KC
			A = uMax/(M_PI*2.0*f), //umax
			totalPoints=NavierStokesSolver::B.totalPoints,
			xold= NavierStokesSolver::B.midX,
			unew,
			xnew;

	//calc new velocity and position
	xnew = A*cos(2*M_PI*f*t);
	unew = -A*2*M_PI*f*sin(2*M_PI*f*t);
	NavierStokesSolver::B.centerVelocityU = unew;
	NavierStokesSolver::B.midX = xnew;

	const int blocksize = 256;
	dim3 grid( int( (totalPoints)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	kernels::update_body_viv<<<grid,block>>>(x_r, uB_r, xnew-xold, unew, totalPoints);
}

void oscCylinder::stepTime()
{
	if (timeStep == 0)
		initialise();

	NavierStokesSolver::generateRHS1();
	//arrayprint(tags,"tags","x");
	//arrayprint(tagsP,"tagsP","p");
	//arrayprint(distance_from_u_to_body,"dub","p");
	//arrayprint(distance_from_v_to_body,"dvb","p");
	NavierStokesSolver::solveIntermediateVelocity();
	//arrayprint(uhat,"uhat","x");

	NavierStokesSolver::generateRHS2();
	//arrayprint(rhs2,"rhs2","p");
	NavierStokesSolver::solvePoisson();

	NavierStokesSolver::velocityProjection();
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

