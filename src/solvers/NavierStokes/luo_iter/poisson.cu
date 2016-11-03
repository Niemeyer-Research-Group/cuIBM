/***************************************************************************//**
 * \file
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief
 */
#include <solvers/NavierStokes/luo_iter.h>

#include <solvers/NavierStokes/luoIBM/kernels/intermediatePressure.h> //intermediatepressure_luo
#include <solvers/NavierStokes/luo_iter/kernels/weight.h> //alpha_p
#include <solvers/NavierStokes/luo_base/kernels/biLinearInterpolation.h> //interpolatepressureto_node
#include <solvers/NavierStokes/luo_iter/kernels/pressure.h> //size_lhs2, update_rhs2
#include <solvers/NavierStokes/luo_iter/kernels/LHS2.h> //lhs2_mid_iter
#include <solvers/NavierStokes/NavierStokes/kernels/LHS2.h> //lhs2 BC

void luo_iter::poisson_setup()
{
	logger.startTimer("Poisson Setup");

	const int blocksize = 256;
	dim3 grid( int( (numP-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	//generate RHS2
	kernels::intermediatePressure_luo<<<grid,block>>>(rhs2_r, uhat_r, ym_r, yp_r, xm_r, xp_r, dx_r, dy_r, nx, ny);

	//calculate alpha
	poisson_alpha();

	//interpolation setup
	poisson_interpolation_setup();

	//size lhs
	poisson_size_lhs();

	//calculate lhs
	poisson_calculate_lhs();

	//sort
	LHS2.sort_by_row_and_column();

	//update preconditioner
	if (timeStep == 0)
		PC.generate2(LHS2, (*paramDB)["PoissonSolve"]["preconditioner"].get<preconditionerType>());
	else
		PC.update2(LHS2);
	//update rhs
	poisson_update_rhs();
	logger.stopTimer("Poisson Setup");
}

void luo_iter::poisson_alpha()
{
	const int blocksize = 256;

	dim3 grid( int( (numP-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	cusp::blas::fill(alpha,0);//I don't think this is needed but it cleans up the alpha array

	kernels::alpha_p<<<grid,block>>>(alpha_r, ghostTagsP_r, hybridTagsP_r, yu_r, xv_r,
										body_intercept_p_x_r, body_intercept_p_y_r, image_point_p_x_r, image_point_p_y_r,
										B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny);
}

void luo_iter::poisson_interpolation_setup()
{
	const int blocksize = 256;
	dim3 grid( int( (B.numCellsXHost*B.numCellsYHost-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	cusp::blas::fill(ns_rhs,0);
	cusp::blas::fill(interp_rhs,0);
	cusp::blas::fill(index1,0);//index and qcoef fills are not needed but make it cleaner
	cusp::blas::fill(index2,0);
	cusp::blas::fill(index3,0);
	cusp::blas::fill(index4,0);
	cusp::blas::fill(q1coef,0);
	cusp::blas::fill(q2coef,0);
	cusp::blas::fill(q3coef,0);
	cusp::blas::fill(q4coef,0);

	kernels::interpolatePressureToGhostNode<<<grid,block>>>(pressure_r, false, u_r, ghostTagsP_r, B.x_r, B.y_r, dpdn_r,
																B.uB_r, B.uBk_r, B.vB_r, B.vBk_r, yu_r, yv_r, xu_r, xv_r,
																body_intercept_p_x_r, body_intercept_p_y_r, image_point_p_x_r, image_point_p_y_r,  body_intercept_p_r,
																B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny, dt,
																index1_r, index2_r, index3_r, index4_r,
																q1coef_r, q2coef_r, q3coef_r, q4coef_r,
																a0_r, a1_r, a2_r, a3_r,
																x1_p_r, x2_p_r, x3_p_r, x4_p_r,
																y1_p_r, y2_p_r, y3_p_r, y4_p_r,
																q1_p_r, q2_p_r, q3_p_r, q4_p_r);

	kernels::interpolatePressureToHybridNode<<<grid,block>>>(pressure_r, pressureStar_r, u_r, hybridTagsP_r, B.x_r, B.y_r,
																B.uB_r, B.uBk_r, B.vB_r, B.vBk_r, yu_r, yv_r, xu_r, xv_r,
																body_intercept_p_x_r, body_intercept_p_y_r, image_point_p_x_r, image_point_p_y_r,
																B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny, dt,
																index1_r, index2_r, index3_r, index4_r,
																q1coef_r, q2coef_r, q3coef_r, q4coef_r,
																a0_r, a1_r, a2_r, a3_r,
																x1_p_r, x2_p_r, x3_p_r, x4_p_r,
																y1_p_r, y2_p_r, y3_p_r, y4_p_r,
																q1_p_r, q2_p_r, q3_p_r, q4_p_r);
	//testInterpP();
}

void luo_iter::poisson_size_lhs()
{
	cusp::blas::fill(count, 0);

	dim3 grid( 1, 1);
	dim3 block(1, 1);
	kernels::size_LHS2<<<grid,block>>>(hybridTagsP_r, count_r, B.startI_r, B.startJ_r, B.numCellsXHost, B.numCellsYHost, nx, ny);

	thrust::device_vector<int>::iterator iter = thrust::max_element(count.begin(),count.end());
	unsigned int position = iter - count.begin();
	double max_val = *iter;

	LHS2.resize(numP,numP, nx*ny*5 - 2*nx-ny*2 + max_val);

	LHS2_row_r = thrust::raw_pointer_cast( &(LHS2.row_indices[0]) );
	LHS2_col_r = thrust::raw_pointer_cast( &(LHS2.column_indices[0]) );
	LHS2_val_r = thrust::raw_pointer_cast( &(LHS2.values[0]) );
}

void luo_iter::poisson_calculate_lhs()
{
	const int blocksize = 256;

	dim3 grid( int( (numP-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	cusp::blas::fill(LHS2.values,0);

	kernels::LHS2_mid_iter<<<grid,block>>>(LHS2_row_r, LHS2_col_r, LHS2_val_r, dx_r, dy_r, nx, ny, dt,
											count_r, ns_rhs_r, interp_rhs_r, hybridTagsP_r, ghostTagsP_r,
											alpha_r, dpdn_r,
											index1_r, index2_r, index3_r, index4_r,
											q1coef_r, q2coef_r, q3coef_r, q4coef_r,
											q1_p_r, q2_p_r, q3_p_r, q4_p_r, timeStep);
	kernels::LHS2_BC<<<grid,block>>>(LHS2_row_r, LHS2_col_r, LHS2_val_r, dx_r, dy_r, nx,ny,dt);
}

void luo_iter::poisson_update_rhs()
{
	const int blocksize = 256;
	dim3 grid( int( (numP-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	kernels::update_rhs2<<<grid, block>>>(rhs2_r, ns_rhs_r, interp_rhs_r, nx, ny);
}
