/***************************************************************************//**
 * \file
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernals that setup the prerequisites to solve the poission equation
 */
#include <solvers/NavierStokes/luoIBM.h>

#include <solvers/NavierStokes/luoIBM/kernels/intermediatePressure.h>
#include <solvers/NavierStokes/luoIBM/kernels/LHS2.h>
#include <solvers/NavierStokes/NavierStokes/kernels/LHS2.h>
#include <solvers/NavierStokes/luoIBM/kernels/weight.h>//weighting function
#include <solvers/NavierStokes/luoIBM/kernels/intermediateVelocity.h>
#include <solvers/NavierStokes/luoIBM/kernels/biLinearInterpolation.h>

void luoIBM::generateRHS2()
{
	NavierStokesSolver::logger.startTimer("RHS2");

	const int blocksize = 256;

	dim3 grid( int( (nx*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	preRHS2Interpolation();

	kernels::intermediatePressure_luo<<<grid,block>>>(rhs2_r, uhat_r, ym_r, yp_r, xm_r, xp_r, dx_r, dy_r, nx, ny);
	NavierStokesSolver::logger.stopTimer("RHS2");
}

void luoIBM::generateLHS2()
{
	NavierStokesSolver::logger.startTimer("LHS2");

	const int blocksize = 256;

	dim3 grid( int( (nx*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	cusp::blas::fill(LHS2.values,0);

	kernels::LHS2_mid_luo<<<grid,block>>>(LHS2_row_r, LHS2_col_r, LHS2_val_r, dx_r, dy_r, nx, ny, dt);
	kernels::LHS2_BC<<<grid,block>>>(LHS2_row_r, LHS2_col_r, LHS2_val_r, dx_r, dy_r, nx,ny,dt);
	logger.stopTimer("LHS2");
}

void luoIBM::weightPressure()
{
	logger.startTimer("weightPressure");

	const int blocksize = 256;

	dim3 grid( int( (B.numCellsXHost*B.numCellsYHost-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	cusp::blas::fill(pressureStar,0);//flag not needed

	kernels::interpolatePressureToHybridNode<<<grid,block>>>(pressure_r, pressureStar_r, u_r, hybridTagsP_r, B.x_r, B.y_r,
									B.uB_r, B.uBk_r, B.vB_r, B.vBk_r, yu_r, yv_r, xu_r, xv_r,
									body_intercept_p_x_r, body_intercept_p_y_r, image_point_p_x_r, image_point_p_y_r,
									B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny, dt,
									dudt_r,ududx_r,vdudy_r,dvdt_r,udvdx_r,vdvdy_r,
									a0_r, a1_r, a2_r, a3_r,
									x1_p_r, x2_p_r, x3_p_r, x4_p_r, y1_p_r, y2_p_r, y3_p_r, y4_p_r, q1_p_r, q2_p_r, q3_p_r, q4_p_r, timeStep);
	kernels::weightP<<<grid,block>>>(pressure_r, pressureStar_r, ghostTagsP_r, hybridTagsP_r, yu_r, xv_r,
									body_intercept_p_x_r, body_intercept_p_y_r, image_point_p_x_r, image_point_p_y_r,
									B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny);
	kernels::interpolatePressureToGhostNode<<<grid,block>>>(pressure_r, u_r, ghostTagsP_r, B.x_r, B.y_r,
									B.uB_r, B.uBk_r, B.vB_r, B.vBk_r, yu_r, yv_r, xu_r, xv_r,
									body_intercept_p_x_r, body_intercept_p_y_r, image_point_p_x_r, image_point_p_y_r,  body_intercept_p_r,
									B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny, dt,
									dudt_r,ududx_r,vdudy_r,dvdt_r,udvdx_r,vdvdy_r,
									a0_r, a1_r, a2_r, a3_r,
									x1_p_r, x2_p_r, x3_p_r, x4_p_r, y1_p_r, y2_p_r, y3_p_r, y4_p_r, q1_p_r, q2_p_r, q3_p_r, q4_p_r);
	testInterpP();
	logger.stopTimer("weightPressure");
}

void luoIBM::preRHS2Interpolation()
{
	const int blocksize = 256;

	dim3 grid( int( (B.numCellsXHost*B.numCellsYHost-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	//interpolate uhat to image point and ghost node
	kernels::interpolateVelocityToGhostNodeX<<<grid,block>>>(uhat_r, ghostTagsUV_r, B.x_r, B.y_r, B.uB_r, yu_r, xu_r,
																body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
																B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny,
																x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r,image_point_u_r);
	kernels::interpolateVelocityToGhostNodeY<<<grid,block>>>(uhat_r, ghostTagsUV_r, B.x_r, B.y_r, B.vB_r, yv_r, xv_r,
																body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
																B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny,
																x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r,image_point_u_r);
}
