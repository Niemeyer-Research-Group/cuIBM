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
{luo_base::writeData();}

void luo_iter::writeCommon()
{luo_base::writeCommon();}

void luo_iter::initialise()
{
	luo_base::initialise();
	luo_iter::cast();

	parameterDB  &db = *paramDB;

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
			uold,
			xnew;

	//xnew = -1/(2*M_PI)*sin(2*M_PI*f*t);
	//unew = -f*cos(2*M_PI*f*t);
	xnew = xCoeff*sin(2*M_PI*f*t + xPhase);
	unew = uCoeff*cos(2*M_PI*f*t + uPhase);
	uold = uCoeff*cos(2*M_PI*f*(t-dt)+uPhase);
	std::cout<<unew<<"\n";
	std::cout<<uold<<"\n";

	B.centerVelocityU0 = unew;
	B.midX0 = xnew;
	B.centerVelocityU = unew;
	B.midX = xnew;

	const int blocksize = 256;
	dim3 grid( int( (totalPoints)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	//update position/velocity for current values
	kernels::update_body_viv<<<grid,block>>>(B.x_r, B.uB_r, xnew-xold, unew, totalPoints);
	//set position/velocity for old values
	kernels::initialise_old<<<grid,block>>>(B.uBk_r,uold,totalPoints);
}

void luo_iter::_intermediate_velocity()
{
	intermediate_velocity_setup();
	/*int index = 0;
	for (int i = 0; i < numUV*5; i++)
	{
		if (ghostTagsUV[LHS1.row_indices[i]]>0)
		{
			if (LHS1.row_indices[i]>index)
			{
				std::cout<<"\n";
				index = LHS1.row_indices[i];
			}
			std::cout<<LHS1.row_indices[i]<<"\t";
			std::cout<<LHS1.column_indices[i]<<"\t";
			std::cout<<LHS1.values[i]<<std::endl;
		}
		if (LHS1.row_indices[i]>70532)
			break;
	}*/
	solveIntermediateVelocity();
	testInterpX();
}

void luo_iter::_pressure()
{
	poisson_setup();
	solvePoisson();
	/*int index = 0, ip, I;
	for (int i = 0; i < numUV*5; i++)
	{
		ip = LHS2.row_indices[i];
		I = ip%nx;
		if (hybridTagsP[LHS2.row_indices[i]]>0 && I <160)
		{
			if (LHS2.row_indices[i] >index)
			{
				std::cout<<"\n";
				index = LHS2.row_indices[i];
			}
			std::cout<<LHS2.row_indices[i]<<"\t";
			std::cout<<LHS2.column_indices[i]<<"\t";
			std::cout<<LHS2.values[i]<<std::endl;
		}
		if (LHS2.row_indices[i]>70000)
			break;
	}*/

	arrayprint(rhs2,"rhs2","p",-1);
	arrayprint(ghostTagsP,"ghostp","p",-1);
	arrayprint(hybridTagsP,"hybridP","p",-1);
	arrayprint(pressure, "p","p",-1);
}
