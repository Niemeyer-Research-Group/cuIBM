/***************************************************************************//**
 * \file
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernels that setup the intermediate velocity solve
 */

#include <solvers/NavierStokes/luo_base.h>

#include <solvers/NavierStokes/FadlunModified/kernels/intermediateVelocity.h> //updateboundary

void luo_base::updateRobinBoundary()
{
	double 	Uinf = 1, //need a better way to enforce these, ie read from yaml file
			Vinf = 1;

	const int blocksize = 256;
	dim3 dimGridBCX( int(ny/blocksize) + 1, 1);
	dim3 dimGridBCY( int(nx/blocksize) + 1, 1);
	dim3 dimBlockBC(blocksize, 1);

	kernels::updateBoundaryX<<<dimGridBCX,dimBlockBC>>>(u_r, xp_r, dx_r, dt, Uinf, nx, ny);
	kernels::updateBoundaryY<<<dimGridBCY,dimBlockBC>>>(u_r, xp_r, dx_r, dt, Vinf, nx, ny);
}
