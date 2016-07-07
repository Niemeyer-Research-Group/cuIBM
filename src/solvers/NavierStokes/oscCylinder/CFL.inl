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

void oscCylinder::calcDistance()
{
	double	*distance_r = thrust::raw_pointer_cast ( &(distance[0]) ),
			*xu_r = thrust::raw_pointer_cast ( &(domInfo->xu[0]) ),
			*xv_r = thrust::raw_pointer_cast ( &(domInfo->xv[0]) ),
			*yu_r = thrust::raw_pointer_cast ( &(domInfo->yu[0]) ),
			*yv_r = thrust::raw_pointer_cast ( &(domInfo->yv[0]) );
	
	int	*ghostTagsUV_r = thrust::raw_pointer_cast ( &(ghostTagsUV[0]) ),
		*ghostTagsP_r  = thrust::raw_pointer_cast ( &(ghostTagsP[0]) );
			
	int nx = domInfo ->nx,
		ny = domInfo ->ny,
		width_i = B.numCellsXHost, //flag this value is only moved to the host once (in B.initialise) if the body is moving too much this could break 
		height_j=B.numCellsYHost;  //this is done because we need the value on the host to calculate the grid size, but copying it to the host every TS is expensive
	
	int *i_start_r = thrust::raw_pointer_cast ( &(B.startI[0]) ),
		*j_start_r = thrust::raw_pointer_cast ( &(B.startJ[0]) );
	double dt = (*paramDB)["simulation"]["dt"].get<double>();

	const int blocksize = 256;
	dim3 grid( int( (width_i*height_j-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	
	kernels::testDistance<<<grid,block>>>(distance_r, ghostTagsUV_r, ghostTagsP_r, xu_r, xv_r, yu_r, yv_r, B.midX, B.midY,
											i_start_r, j_start_r, width_i, nx, ny);
	arrayprint(distance,"distance","p",-1);
}
