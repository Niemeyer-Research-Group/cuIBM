/***************************************************************************//**
 * \file CFL.inl
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernels that setup the intermediate velocity solve
 */

#include <solvers/NavierStokes/oscCylinder/kernels/CFL.h>

void oscCylinder::CFL()
{
	logger.startTimer("CFL");
	double	*u_r = thrust::raw_pointer_cast ( &(u[0]) ),
			*cfl_r = thrust::raw_pointer_cast ( &(cfl[0]) ),
			*dx_r = thrust::raw_pointer_cast ( &(domInfo->dx[0]) ),
			*dy_r = thrust::raw_pointer_cast ( &(domInfo->dy[0]) );
			
	double	dt = (*paramDB)["simulation"]["dt"].get<double>();
		
	int nx = domInfo ->nx,
		ny = domInfo ->ny;
	
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

