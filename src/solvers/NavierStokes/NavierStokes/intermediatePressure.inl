/***************************************************************************//**
 * \file intermediatePressure.inl
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernals that setup the prerequisites to solve the poission equation
 */


#include <solvers/NavierStokes/NavierStokes/kernels/intermediatePressure.h>
#include <solvers/NavierStokes/NavierStokes/kernels/LHS2.h>

void NavierStokesSolver::generateRHS2()
{
	logger.startTimer("RHS2");

	int nx = domInfo ->nx,
		ny = domInfo ->ny;

	double	*rhs2_r	= thrust::raw_pointer_cast( &(rhs2[0]) ),
			*uhat_r	= thrust::raw_pointer_cast( &(uhat[0]) ),
			*dx_r	= thrust::raw_pointer_cast( &(domInfo->dx[0]) ),
			*dy_r	= thrust::raw_pointer_cast( &(domInfo->dy[0]) ),
			*ym_r	= thrust::raw_pointer_cast( &(bc[YMINUS][0]) ),
			*yp_r	= thrust::raw_pointer_cast( &(bc[YPLUS][0]) ),
			*xm_r	= thrust::raw_pointer_cast( &(bc[XMINUS][0]) ),
			*xp_r	= thrust::raw_pointer_cast( &(bc[XPLUS][0]) );

	const int blocksize = 256;
	
	double	dt = (*paramDB)["simulation"]["dt"].get<double>();

	dim3 grid( int( (nx*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	kernels::intermediatePressureNoBody<<<grid,block>>>(rhs2_r, uhat_r, ym_r, yp_r, xm_r, xp_r, dx_r, dy_r, nx, ny);
	pressure_old = pressure;
	logger.stopTimer("RHS2");
}

void NavierStokesSolver::generateLHS2()
{
	logger.startTimer("LHS2");
	int nx = domInfo ->nx,
		ny = domInfo ->ny;

	double	dt 		= (*paramDB)["simulation"]["dt"].get<double>(),
			*dx_r	= thrust::raw_pointer_cast( &(domInfo->dx[0]) ),
			*dy_r	= thrust::raw_pointer_cast( &(domInfo->dy[0]) ),
			*val_r	= thrust::raw_pointer_cast( &(LHS2.values[0]) );

	int		*row_r	= thrust::raw_pointer_cast( &(LHS2.row_indices[0]) ),
			*col_r	= thrust::raw_pointer_cast( &(LHS2.column_indices[0]) );

	const int blocksize = 256;
	dim3 grid( int( (nx*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	cusp::blas::fill(LHS2.values,0);
	kernels::LHS2_mid_nobody<<<grid,block>>>(row_r, col_r, val_r, dx_r, dy_r, nx, ny, dt);
	kernels::LHS2_BC<<<grid,block>>>(row_r, col_r, val_r, dx_r, dy_r, nx,ny,dt);
	logger.stopTimer("LHS2");
}
