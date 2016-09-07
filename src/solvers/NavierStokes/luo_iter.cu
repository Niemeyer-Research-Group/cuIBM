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
}

void luo_iter::_intermediate_velocity()
{
	intermediate_velocity_setup();
	solveIntermediateVelocity();
}

void luo_iter::_pressure()
{
	poisson_setup();
	solvePoisson();
}
