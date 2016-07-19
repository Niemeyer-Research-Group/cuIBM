/***************************************************************************//**
 * \file intermediateVelocity.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief Functions that call kernels to solve for the intermediate velocity
 */

#include <solvers/NavierStokes/NavierStokesSolver.h>

#include <solvers/NavierStokes/NavierStokes/kernels/intermediateVelocity.h>
#include <solvers/NavierStokes/NavierStokes/kernels/LHS1.h>
#include <solvers/NavierStokes/NavierStokes/kernels/N.h>
#include <solvers/NavierStokes/NavierStokes/kernels/L.h>

void NavierStokesSolver::generateRHS1()
{
	//set correct grid and block size
	const int blocksize = 256;

	dim3 dimGridUV( int( ((nx-1)*ny+nx*(ny-1)-0.5)/blocksize ) + 1, 1);
	dim3 dimBlockUV(blocksize, 1);
	dim3 dimGridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 dimBlockU(blocksize, 1);
	dim3 dimGridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 dimBlockV(blocksize, 1);

	//boundary update goes here

	Nold = N;
	generateN();

	generateL();

	generateBC1();

	kernels::generateRHS<<<dimGridUV,dimBlockUV>>>(rhs1_r, L_r, Nold_r, N_r, u_r, bc1_r, dt, nx, ny);
}

/**
 * \brief Fills the Convection matrix N
 * \ param nx number of pressure nodes in the x direction
 * \ param ny number of pressure nodes in the y direction
 * \ param dx array containing the X-direction cell widths
 * \ param dy array containing the Y-direction cell heights
 * \ param N  convection matrix (stored as an array)
 * \ param u  velocity matrix (stored as an array)
 */
void NavierStokesSolver::generateN()
{
	logger.startTimer("Advection Terms");

	const int blocksize = 256;

	dim3 gridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 blockU(blocksize, 1);
	dim3 gridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 blockV(blocksize, 1);

	kernels::Nmidx<<<gridU, blockU>>>(N_r,u_r,dx_r,dy_r,nx,ny);
	kernels::Nbcx<<<gridU, blockU>>>(N_r, u_r, dx_r, dy_r, ym_r, yp_r, xm_r, xp_r, nx, ny);
	kernels::Nmidy<<<gridV, blockV>>>(N_r,u_r,dx_r,dy_r,nx,ny);
	kernels::Nbcy<<<gridV, blockV>>>(N_r, u_r, dx_r, dy_r, ym_r, yp_r, xm_r, xp_r, nx, ny);

	logger.stopTimer("Advection Terms");
}

/**
 * \brief Fills the Laplacian matrix L
 * \ param nx number of pressure nodes in the x direction
 * \ param ny number of pressure nodes in the y direction
 * \ param nu viscosity: effectivly the inverse Reynolds number
 * \ param dx array containing the X-direction cell widths
 * \ param dy array containing the Y-direction cell heights
 * \ param L  laplacian matrix (stored as an array)
 * \ param u  velocity matrix (stored as an array)
 */
void NavierStokesSolver::generateL()
{
	logger.startTimer("Diffusion Terms");

	const int blocksize = 256;

	dim3 gridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 blockU(blocksize, 1);
	dim3 gridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 blockV(blocksize, 1);

	kernels::Lmidx<<<gridU, blockU>>>(L_r,u_r,dx_r,dy_r,nx,ny,nu);
	kernels::Lbcx<<<gridU, blockU>>>(L_r, u_r, dx_r, dy_r, ym_r, yp_r, xm_r, xp_r, nx, ny, nu);
	kernels::Lmidy<<<gridV, blockV>>>(L_r,u_r,dx_r,dy_r,nx,ny,nu);
	kernels::Lbcy<<<gridV, blockV>>>(L_r, u_r, dx_r, dy_r, ym_r, yp_r, xm_r, xp_r, nx, ny, nu);

	logger.stopTimer("Diffusion Terms");
}

void NavierStokesSolver::generateLHS1()
{
	logger.startTimer("LHS1");

	const int blocksize = 256;

	dim3 gridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 blockU(blocksize, 1);
	dim3 gridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 blockV(blocksize, 1);

	kernels::LHS_mid_X_nobody<<<gridU,blockU>>>(LHS1_row_r, LHS1_col_r, LHS1_val_r, dx_r, dy_r, dt, nu, nx, ny);
	kernels::LHS_mid_Y_nobody<<<gridV,blockV>>>(LHS1_row_r, LHS1_col_r, LHS1_val_r, dx_r, dy_r, dt, nu, nx, ny);
	kernels::LHS_BC_X<<<gridU,blockU>>>(LHS1_row_r, LHS1_col_r, LHS1_val_r, dx_r, dy_r, dt, nu, nx, ny);
	kernels::LHS_BC_Y<<<gridV,blockV>>>(LHS1_row_r, LHS1_col_r, LHS1_val_r, dx_r, dy_r, dt, nu, nx, ny);

	logger.stopTimer("LHS1");
}

void NavierStokesSolver::generateBC1()
{
	logger.startTimer("BC1");

	const int blocksize = 256;

	dim3 gridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 blockU(blocksize, 1);
	dim3 gridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 blockV(blocksize, 1);

	kernels::bc1X<<<gridU, blockU>>>(u_r, bc1_r, ym_r, yp_r, xm_r, xp_r, dx_r, dy_r, nu, dt, nx, ny);
	kernels::bc1Y<<<gridV, blockV>>>(u_r, bc1_r, ym_r, yp_r, xm_r, xp_r, dx_r, dy_r, nu, dt, nx, ny);

	logger.stopTimer("BC1");
}
