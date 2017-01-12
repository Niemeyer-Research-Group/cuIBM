/***************************************************************************//**
 * \file intermediatePressure.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief Functions that call kernels to solve the Poisson equation
 */

#include <solvers/NavierStokes/NavierStokesSolver.h>

#include <solvers/NavierStokes/NavierStokes/kernels/intermediatePressure.h>
#include <solvers/NavierStokes/NavierStokes/kernels/LHS2.h>

void NavierStokesSolver::generateRHS2()
{
	logger.startTimer("Poisson Setup");

	const int blocksize = 256;

	dim3 grid( int( (nx*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	kernels::intermediatePressureNoBody<<<grid,block>>>(rhs2_r, uhat_r, ym_r, yp_r, xm_r, xp_r, dx_r, dy_r, nx, ny);

	logger.stopTimer("Poisson Setup");
}

void NavierStokesSolver::generateLHS2()
{
	const int blocksize = 256;

	dim3 grid( int( (nx*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	cusp::blas::fill(LHS2.values,0);

	kernels::LHS2_mid_nobody<<<grid,block>>>(LHS2_row_r, LHS2_col_r, LHS2_val_r, dx_r, dy_r, nx, ny, dt);
	kernels::LHS2_BC<<<grid,block>>>(LHS2_row_r, LHS2_col_r, LHS2_val_r, dx_r, dy_r, nx,ny,dt);
}
