/***************************************************************************//**
 * \file  FSI.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \based on code by Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class FSI.
 */

#include <solvers/NavierStokes/oscCylinder/kernels/structure.h>
#include "FSI.h"
#include <solvers/NavierStokes/NavierStokes/kernels/intermediatePressure.h>
#include <sys/stat.h>

/**
 * \brief Constructor. Copies the database and information about the computational grid.
 *
 * \param pDB database that contains all the simulation parameters
 * \param dInfo information related to the computational grid
 *//*
FSI::FSI(parameterDB *pDB, domain *dInfo)
{
	paramDB = pDB;
	domInfo = dInfo;
}

void FSI::writeData()
{
	parameterDB  &db = *NavierStokesSolver::paramDB;
	double dt  = db["simulation"]["dt"].get<double>();
	NavierStokesSolver::logger.startTimer("output");
	NavierStokesSolver::writeCommon();
	forceFile << timeStep*dt << '\t' << B.forceX[0] << '\t' << B.forceY[0] << std::endl;
	logger.stopTimer("output");
}

void FSI::updateSolver()
{
	NavierStokesSolver::B.calculateBoundingBoxes(*NavierStokesSolver::paramDB, *NavierStokesSolver::domInfo);
	NavierStokesSolver::tagPoints();
	NavierStokesSolver::generateLHS1();//is this needed?
	NavierStokesSolver::generateLHS2();

	NavierStokesSolver::logger.startTimer("Preconditioner");
	if (NavierStokesSolver::iterationCount2 > 20)
	{
		//NavierStokesSolver::PC.update(NavierStokesSolver::LHS1,NavierStokesSolver::LHS2);
		//NavierStokesSolver::PC1->update(NavierStokesSolver::LHS1);
		//NavierStokesSolver::PC2->update(NavierStokesSolver::LHS2);
	}
	NavierStokesSolver::logger.stopTimer("Preconditioner");
}

void FSI::moveBody()
{
	parameterDB  &db = *NavierStokesSolver::paramDB;
	NavierStokesSolver::calculateForce();

	double *y_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::B.y[0]) ),
		   *vB_r= thrust::raw_pointer_cast( &(NavierStokesSolver::B.vB[0]) );
	double	Cy	= NavierStokesSolver::B.forceY[0]*2.0,
			U	= NavierStokesSolver::bc[XMINUS][0],
			Mred= 2.0,
			Ured= db["flow"]["Ured"].get<double>(),
			dt	= db["simulation"]["dt"].get<double>(),
			totalPoints=NavierStokesSolver::B.totalPoints,
			vold= NavierStokesSolver::B.centerVelocityV,
			yold= NavierStokesSolver::B.midY,
			vnew,
			ynew;

	double a = dt*M_PI*M_PI*4/(Ured*Ured),
		   b = dt*dt*2*M_PI*M_PI/(Ured*Ured);

	//calc updated velocity
	vnew = (vold - a*(yold+ dt/2*vold) + dt*Cy/2/Mred)/(1+b);
	ynew = yold + dt/2*(vnew + vold);
	NavierStokesSolver::B.centerVelocityV = vnew;
	NavierStokesSolver::B.midY = ynew;

	const int blocksize = 256;
	dim3 grid( int( (totalPoints)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	kernels::update_body_viv<<<grid,block>>>(y_r, vB_r, ynew-yold, vnew, totalPoints);
}

void FSI::moveBodySC()
{
	parameterDB  &db = *NavierStokesSolver::paramDB;
	NavierStokesSolver::calculateForce();

	double *y_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::B.y[0]) ),
		   *vB_r= thrust::raw_pointer_cast( &(NavierStokesSolver::B.vB[0]) );
	double	Cy	= NavierStokesSolver::B.forceY[0]*2.0,
			U	= NavierStokesSolver::bc[XMINUS][0],
			Mred= 2.0,
			Ured= db["simulation"]["Ured"].get<double>(),
			dt	= db["simulation"]["dt"].get<double>(),
			totalPoints=NavierStokesSolver::B.totalPoints,
			vold= NavierStokesSolver::B.centerVelocityV0,
			yold= NavierStokesSolver::B.midY0,
			vnew,
			ynew,
			relaxation_coeff = 0.75;

	double a = dt*M_PI*M_PI*4/(Ured*Ured),
		   b = dt*dt*2*M_PI*M_PI/(Ured*Ured);

	//calc updated velocity
	vnew = (vold - a*(yold+ dt/2*vold) + dt*Cy/2/Mred)/(1+b);
	vnew = relaxation_coeff * vnew + (1-relaxation_coeff) * NavierStokesSolver::B.centerVelocityV; //relax velocity
	ynew = yold + dt/2*(vnew + vold);

	NavierStokesSolver::B.centerVelocityV = vnew;
	NavierStokesSolver::B.midY = ynew;

	const int blocksize = 256;
	dim3 grid( int( (totalPoints)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	kernels::update_body_viv<<<grid,block>>>(y_r, vB_r, ynew-yold, vnew, totalPoints);
}

void FSI::stepTime()
{
	LC();
}

void FSI::LC()
{
	NavierStokesSolver::generateRHS1();
	NavierStokesSolver::solveIntermediateVelocity();

	NavierStokesSolver::generateRHS2();
	NavierStokesSolver::solvePoisson();

	NavierStokesSolver::velocityProjection();

	//Release the body after a certain timestep
	if (NavierStokesSolver::timeStep >= (*NavierStokesSolver::paramDB)["simulation"]["startStep"].get<int>())
	{
		std::cout<<"5.1\n\n";
		moveBody();
		updateSolver();
	}
	NavierStokesSolver::timeStep++;

	if (NavierStokesSolver::timeStep > 140)
	{
		//arrayprint(rhs1,"rhs1","x");
		//arrayprint(uhat,"uhat","x");
		//arrayprint(uhat,"vhat","y");

		//arrayprint(rhs2,"rhs2","p");
		//arrayprint(pressure,"pressure","p");

		//arrayprint(u,"u","x");
		//arrayprint(u,"v","y");
	}
}

void FSI::SC()
{
	NavierStokesSolver::B.centerVelocityV0 = NavierStokesSolver::B.centerVelocityV;
	NavierStokesSolver::B.midY0 = NavierStokesSolver::B.midY;
	int count = 0;
	do
	{
		NavierStokesSolver::B.forceYk[0] = NavierStokesSolver::B.forceY[0];
		NavierStokesSolver::generateRHS1();
		NavierStokesSolver::solveIntermediateVelocity();

		NavierStokesSolver::generateRHS2();
		NavierStokesSolver::solvePoisson();

		NavierStokesSolver::velocityProjection();
		//Release the body after a certain timestep
		if (NavierStokesSolver::timeStep >= (*NavierStokesSolver::paramDB)["simulation"]["startStep"].get<int>())
		{
			moveBodySC();
			updateSolver();
		}
		count += 1;
	}
	while (fabs(NavierStokesSolver::B.forceY[0]- NavierStokesSolver::B.forceYk[0]) > 0.0001);
	if (count > 1)
		std::cout<<count<<"\n";
	std::cout<<NavierStokesSolver::timeStep<<"\n";

	NavierStokesSolver::timeStep++;

}*/
/*
void FSI::callTest()
{
	int nx = NavierStokesSolver::domInfo ->nx,
		ny = NavierStokesSolver::domInfo ->ny;
	double	dt = (*NavierStokesSolver::paramDB)["simulation"]["dt"].get<double>();
	const int blocksize = 256;
	dim3 grid( int( (nx*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	double	*test_r = thrust::raw_pointer_cast( &(NavierStokesSolver::test[0]) ),
			*uhat_r = thrust::raw_pointer_cast( &(NavierStokesSolver::uhat[0]) ),
			*pressure_r = thrust::raw_pointer_cast( &(NavierStokesSolver::pressure[0]) ),
			*dx_r = thrust::raw_pointer_cast( &(NavierStokesSolver::domInfo->dx[0]) );


	kernels::testkernel<<<grid,block>>>(test_r, uhat_r, pressure_r, dx_r, dt, nx, ny);
}*/
