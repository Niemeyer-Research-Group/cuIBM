/***************************************************************************//**
 * \file
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernels that setup the intermediate velocity solve
 */

#include <solvers/NavierStokes/luoIBM.h>

#include <solvers/NavierStokes/NavierStokes/kernels/intermediateVelocity.h> //generaterhs1
#include <solvers/NavierStokes/luoIBM/kernels/intermediateVelocity.h>//updaterhs1_luo, zeroInside
#include <solvers/NavierStokes/luoIBM/kernels/LHS1.h> //generatelhs_luo _mid
#include <solvers/NavierStokes/NavierStokes/kernels/LHS1.h> //lhs_bc
#include <solvers/NavierStokes/FadlunModified/kernels/intermediateVelocity.h> //updateboundary
#include <solvers/NavierStokes/luo_base/kernels/biLinearInterpolation.h> //interpolate
#include <solvers/NavierStokes/luoIBM/kernels/weight.h>//weighting function

void luoIBM::generateRHS1()
{
	logger.startTimer("RHS1 Setup");

	//set correct grid and block size
	const int blocksize = 256;

	dim3 dimGridUV( int( ((nx-1)*ny+nx*(ny-1)-0.5)/blocksize ) + 1, 1);
	dim3 dimBlockUV(blocksize, 1);
	dim3 dimGridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 dimBlockU(blocksize, 1);
	dim3 dimGridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 dimBlockV(blocksize, 1);

	logger.stopTimer("RHS1 Setup");
	//update right boundary of the domain for the convective boundary condition
	updateRobinBoundary();
	//set Nold to N
	Nold = N;

	//interpolate u and v to image points, then again to ghost nodes
	preRHS1Interpolation();

	//calculate explicit advection terms
	generateN();

	//calculate explicit diffusion terms
	generateL();

	//calculate boundary terms
	generateBC1();

	logger.startTimer("RHS1 Setup");
	//sum rhs components
	kernels::generateRHS<<<dimGridUV,dimBlockUV>>>(rhs1_r, L_r, Nold_r, N_r, u_r, bc1_r, dt, nx, ny);
	logger.stopTimer("RHS1 Setup");
}

void luoIBM::generateLHS1()
{
	NavierStokesSolver::logger.startTimer("LHS1");

	const int blocksize = 256;
	dim3 gridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 blockU(blocksize, 1);
	dim3 gridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 blockV(blocksize, 1);

	kernels::LHS1_mid_luo_X<<<gridU,blockU>>>(LHS1_row_r, LHS1_col_r, LHS1_val_r, ghostTagsUV_r, dx_r, dy_r, dt, nu, nx, ny);//flag lhs1 mid luo is the same as lhs1 mid nobody from NS
	kernels::LHS1_mid_luo_Y<<<gridV,blockV>>>(LHS1_row_r, LHS1_col_r, LHS1_val_r, ghostTagsUV_r, dx_r, dy_r, dt, nu, nx, ny);

	kernels::LHS_BC_X<<<gridU,blockU>>>(LHS1_row_r, LHS1_col_r, LHS1_val_r, dx_r, dy_r, dt, nu, nx, ny);
	kernels::LHS_BC_Y<<<gridV,blockV>>>(LHS1_row_r, LHS1_col_r, LHS1_val_r, dx_r, dy_r, dt, nu, nx, ny);

	NavierStokesSolver::logger.stopTimer("LHS1");
}

void luoIBM::weightUhat()
{
	logger.startTimer("weightUhat");

	const int blocksize = 256;
	dim3 grid( int( (B.numCellsXHost*B.numCellsYHost-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	kernels::weightX<<<grid,block>>>(uhat_r, ustar_r, ghostTagsUV_r, hybridTagsUV_r, yu_r, xu_r,
									body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
									B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny);
	kernels::weightY<<<grid,block>>>(uhat_r, ustar_r, ghostTagsUV_r, hybridTagsUV_r, yv_r, xv_r,
									body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
									B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny);
	logger.stopTimer("weightUhat");
}

void luoIBM::preRHS1Interpolation()
{
	logger.startTimer("RHS1 Interpolation");

	const int blocksize = 256;

	dim3 grid( int( (B.numCellsXHost*B.numCellsYHost-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	kernels::interpolateVelocityToGhostNodeX<<<grid,block>>>(u_r, ghostTagsUV_r, B.x_r, B.y_r, B.uB_r, yu_r, xu_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny,
													x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r,image_point_u_r);
	kernels::interpolateVelocityToGhostNodeY<<<grid,block>>>(u_r, ghostTagsUV_r, B.x_r, B.y_r, B.vB_r, yv_r, xv_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny,
													x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r,image_point_u_r);
	zeroVelocity();
	//interpolate velocity to hybrid node
	//zero ustar
	kernels::interpolateVelocityToHybridNodeX<<<grid,block>>>(u_r, ustar_r, hybridTagsUV_r, detA_r,
																B.x_r, B.y_r, B.uB_r, yu_r, xu_r,
																body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
																B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny,
																index1_r, index2_r, index3_r, index4_r,
																b11_r, b12_r, b13_r, b14_r,
																b21_r, b22_r, b23_r, b24_r,
																b31_r, b32_r, b33_r, b34_r,
																b41_r, b42_r, b43_r, b44_r,
																x1_r, x2_r ,x3_r ,x4_r,
																y1_r, y2_r, y3_r, y4_r,
																q1_r, q2_r, q3_r, q4_r, image_point_u_r);
	kernels::interpolateVelocityToHybridNodeY<<<grid,block>>>(u_r, ustar_r, hybridTagsUV_r, detA_r,
																B.x_r, B.y_r, B.vB_r, yv_r, xv_r,
																body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
																B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny,
																index1_r, index2_r, index3_r, index4_r,
																b11_r, b12_r, b13_r, b14_r,
																b21_r, b22_r, b23_r, b24_r,
																b31_r, b32_r, b33_r, b34_r,
																b41_r, b42_r, b43_r, b44_r,
																x1_r, x2_r ,x3_r ,x4_r,
																y1_r, y2_r, y3_r, y4_r,
																q1_r, q2_r, q3_r, q4_r, image_point_u_r);
	//testInterpX();
	//testInterpY();
	logger.stopTimer("RHS1 Interpolation");
}

void luoIBM::zeroVelocity()
{
	const int blocksize = 256;

	dim3 grid( int( (nx*(ny-1) + (nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	kernels::zeroInside<<<grid,block>>>(ghostTagsUV_r, u_r, nx*(ny-1) + (nx-1)*ny);
}
