/***************************************************************************//**
 * \file intermediatePressure.inl
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernals that setup the prerequisites to solve the poission equation
 */

#include <solvers/NavierStokes/luoIBM/kernels/intermediatePressure.h>
#include <solvers/NavierStokes/luoIBM/kernels/LHS2.h>
#include <solvers/NavierStokes/NavierStokes/kernels/LHS2.h>

void luoIBM::generateRHS2()
{
	NavierStokesSolver::logger.startTimer("RHS2");

	int nx = NavierStokesSolver::domInfo ->nx,
		ny = NavierStokesSolver::domInfo ->ny;

	double	*rhs2_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::rhs2[0]) ),
			*uhat_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::uhat[0]) ),
			*dx_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::domInfo->dx[0]) ),
			*dy_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::domInfo->dy[0]) ),
			*ym_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::bc[YMINUS][0]) ),
			*yp_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::bc[YPLUS][0]) ),
			*xm_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::bc[XMINUS][0]) ),
			*xp_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::bc[XPLUS][0]) ),
			*distance_from_u_to_body_r = thrust::raw_pointer_cast( &(distance_from_u_to_body[0]) ),
			*distance_from_v_to_body_r = thrust::raw_pointer_cast( &(distance_from_v_to_body[0]) );

	int		*ghostTagsP_r= thrust::raw_pointer_cast( &(ghostTagsP[0]) ),
			*ghostTagsUV_r=thrust::raw_pointer_cast( &(ghostTagsUV[0]) ),
			*hybridTagsP_r=thrust::raw_pointer_cast( &(hybridTagsP[0]) );
	const int blocksize = 256;

	parameterDB  &db = *NavierStokesSolver::paramDB;
	double	dt = db["simulation"]["dt"].get<double>();

	dim3 grid( int( (nx*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	kernels::intermediatePressure_luo<<<grid,block>>>(rhs2_r, uhat_r, ghostTagsP_r, hybridTagsP_r, ghostTagsUV_r, distance_from_u_to_body_r, distance_from_v_to_body_r, ym_r, yp_r, xm_r, xp_r, dx_r, dy_r, nx, ny);
	NavierStokesSolver::pressure_old = NavierStokesSolver::pressure;
	NavierStokesSolver::logger.stopTimer("RHS2");
}

void luoIBM::generateLHS2()
{
	NavierStokesSolver::logger.startTimer("LHS2");
	int nx = domInfo ->nx,
		ny = domInfo ->ny;
	parameterDB  &db = *NavierStokesSolver::paramDB;
	
	double	dt 		= db["simulation"]["dt"].get<double>(),
			*dx_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::domInfo->dx[0]) ),
			*dy_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::domInfo->dy[0]) ),
			*val_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::LHS2.values[0]) ),
			*dxu_r	= thrust::raw_pointer_cast( &(distance_from_u_to_body[0]) ),
			*dyv_r	= thrust::raw_pointer_cast( &(distance_from_v_to_body[0]) );

	int		*row_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::LHS2.row_indices[0]) ),
			*col_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::LHS2.column_indices[0]) ),
			*ghostTagsP_r= thrust::raw_pointer_cast( &(ghostTagsP[0]) ),
			*hybridTagsP_r= thrust::raw_pointer_cast( &(hybridTagsP[0]) );

	const int blocksize = 256;
	dim3 grid( int( (nx*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	cusp::blas::fill(LHS2.values,0);
	kernels::LHS2_mid_luo<<<grid,block>>>(row_r, col_r, val_r, dxu_r, dyv_r, ghostTagsP_r, hybridTagsP_r, dx_r, dy_r, nx, ny, dt);
	kernels::LHS2_BC<<<grid,block>>>(row_r, col_r, val_r, dx_r, dy_r, nx,ny,dt);
	logger.stopTimer("LHS2");
}
