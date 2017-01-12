/***************************************************************************//**
 * \file
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernals that setup the prerequisites to solve the poission equation
 */

#include <solvers/NavierStokes/fadlunModified.h>

#include <solvers/NavierStokes/FadlunModified/kernels/intermediatePressure.h>
#include <solvers/NavierStokes/FadlunModified/kernels/LHS2.h>
#include <solvers/NavierStokes/NavierStokes/kernels/LHS2.h>

void fadlunModified::generateRHS2()
{
	const int blocksize = 256;

	dim3 grid( int( (nx*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	kernels::intermediatePressure<<<grid,block>>>(rhs2_r, uhat_r, ghostTagsP_r, hybridTagsP_r, ghostTagsUV_r, distance_from_u_to_body_r, distance_from_v_to_body_r, ym_r, yp_r, xm_r, xp_r, dx_r, dy_r, nx, ny);
	NavierStokesSolver::pressure_old = NavierStokesSolver::pressure;
}

void fadlunModified::generateLHS2()
{
	const int blocksize = 256;

	dim3 grid( int( (nx*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	cusp::blas::fill(LHS2.values,0);

	kernels::LHS2_mid<<<grid,block>>>(LHS2_row_r, LHS2_col_r, LHS2_val_r, distance_from_u_to_body_r, distance_from_v_to_body_r, ghostTagsP_r, hybridTagsP_r, dx_r, dy_r, nx, ny, dt);
	kernels::LHS2_BC<<<grid,block>>>(LHS2_row_r, LHS2_col_r, LHS2_val_r, dx_r, dy_r, nx,ny,dt);
}
