/***************************************************************************//**
 * \file
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernels that setup the intermediate velocity solve
 */
#include <solvers/NavierStokes/NavierStokesSolver.h>

#include <solvers/NavierStokes/NavierStokes/kernels/CFL.h>

void NavierStokesSolver::CFL()
{
	logger.startTimer("CFL");

	const int blocksize = 256;

	dim3 grid( int( (nx*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	kernels::calculateCFL<<<grid,block>>>(cfl_r, u_r, dx_r, dy_r, nx, ny, dt);

	thrust::device_vector<double>::iterator iter = thrust::max_element(cfl.begin(),cfl.end());
	unsigned int position = iter - cfl.begin();
	double max_val = *iter;
	if (max_val > cfl_max)
	{
		cfl_max = max_val;
		cfl_I = position%nx;
		cfl_J = int(position/nx);
		cfl_ts = timeStep;
	}
	logger.stopTimer("CFL");
}

/*
void NavierStokesSolver::calcDistance()
{

	int width_i = B.numCellsXHost, //flag this value is only moved to the host once (in B.initialise) if the body is moving too much this could break
		height_j=B.numCellsYHost;  //this is done because we need the value on the host to calculate the grid size, but copying it to the host every TS is expensive

	const int blocksize = 256;

	dim3 grid( int( (width_i*height_j-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	kernels::testDistance<<<grid,block>>>(distance_r, ghostTagsUV_r, ghostTagsP_r, xu_r, xv_r, yu_r, yv_r, B.midX, B.midY,
											B.startI_r, B.startJ_r, width_i, nx, ny);
	arrayprint(distance,"distance","p",-1);
}
*/
