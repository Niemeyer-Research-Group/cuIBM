/***************************************************************************//**
 * \file
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernels that setup the intermediate velocity solve
 */
#include <solvers/NavierStokes/luo_iter.h>

//kernels
#include <solvers/NavierStokes/NavierStokes/kernels/intermediateVelocity.h> //rhs init
#include <solvers/NavierStokes/luo_base/kernels/intermediateVelocity.h> //set inside velocity
#include <solvers/NavierStokes/luo_iter/kernels/intermediateVelocity.h> //size lhs1, update rhs1
#include <solvers/NavierStokes/luo_base/kernels/biLinearInterpolation.h> //interpolate
#include <solvers/NavierStokes/luo_iter/kernels/weight.h> //alpha
#include <solvers/NavierStokes/NavierStokes/kernels/LHS1.h> //lhs_bc
#include <solvers/NavierStokes/luo_iter/kernels/LHS1.h> //lhs mid

void luo_iter::intermediate_velocity_setup()
{
	logger.startTimer("Intermediate Velocity Setup");

	//set correct grid and block size
	const int blocksize = 256;

	dim3 grid( int( ((nx-1)*ny+nx*(ny-1)-0.5)/blocksize ) + 1, 1);
	dim3 block(blocksize, 1);
	dim3 grid_inside( int( (numUV-0.5)/blocksize ) +1, 1);
	dim3 block_inside(blocksize, 1);

	//update right boundary of the domain for the convective boundary condition
	updateRobinBoundary();

	//set Nold to N
	Nold = N;

	//set inside velocity to the velocity of the body
	kernels::setInsideVelocity<<<grid_inside,block_inside>>>(ghostTagsUV_r, u_r, B.uB_r, B.vB_r, nx, ny);

	//calculate explicit advection terms
	generateN();

	//calculate explicit diffusion terms
	generateL();

	//calculate boundary terms
	generateBC1();

	//sum rhs components
	kernels::generateRHS<<<grid,block>>>(rhs1_r, L_r, Nold_r, N_r, u_r, bc1_r, dt, nx, ny);

	//calculate alpha
	intermediate_velocity_alpha();

	//interpolation setup
	intermediate_velocity_interpolation_setup();

	//size lhs
	intermediate_velocity_size_lhs();

	//calculate lhs
	intermediate_velocity_calculate_lhs();

	//sort
	LHS1.sort_by_row_and_column();

	//update preconditioner
	if (timeStep == 0)
		PC.generate1(LHS1, (*paramDB)["velocitySolve"]["preconditioner"].get<preconditionerType>());
	else
		PC.update1(LHS1);

	//update rhs
	intermediate_velocity_update_rhs();

	logger.stopTimer("Intermediate Velocity Setup");
}

void luo_iter::intermediate_velocity_alpha()
{
	const int blocksize = 256;

	dim3 grid( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	cusp::blas::fill(alpha,0);//I don't think this is needed but it cleans up the alpha array

	kernels::alpha_u<<<grid,block>>>(alpha_r, ghostTagsUV_r, hybridTagsUV_r, yu_r, xu_r,
										body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
										B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny);

	kernels::alpha_v<<<grid,block>>>(alpha_r, ghostTagsUV_r, hybridTagsUV_r, yv_r, xv_r,
										body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
										B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny);
}

void luo_iter::intermediate_velocity_interpolation_setup()
{
	const int blocksize = 256;

	dim3 grid( int( (B.numCellsXHost*B.numCellsYHost-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	cusp::blas::fill(ns_rhs,0);
	cusp::blas::fill(interp_rhs,0);
	//interpolate velocity to image point and ghost node
	kernels::interpolateVelocityToGhostNodeX<<<grid,block>>>(u_r, false, ghostTagsUV_r, B.x_r, B.y_r, B.uB_r, yu_r, xu_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny,
													index1_r, index2_r, index3_r, index4_r,
													q1coef_r, q2coef_r, q3coef_r, q4coef_r,
													x1_r, x2_r, x3_r, x4_r,
													y1_r, y2_r, y3_r, y4_r,
													q1_r, q2_r, q3_r, q4_r, image_point_u_r);
	kernels::interpolateVelocityToGhostNodeY<<<grid,block>>>(u_r, false, ghostTagsUV_r, B.x_r, B.y_r, B.vB_r, yv_r, xv_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny,
													index1_r, index2_r, index3_r, index4_r,
													q1coef_r, q2coef_r, q3coef_r, q4coef_r,
													x1_r, x2_r, x3_r, x4_r,
													y1_r, y2_r, y3_r, y4_r,
													q1_r, q2_r, q3_r, q4_r, image_point_u_r);
	kernels::interpolateVelocityToHybridNodeX<<<grid,block>>>(u_r, ustar_r, hybridTagsUV_r,
																B.x_r, B.y_r, B.uB_r, yu_r, xu_r,
																body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
																B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny,
																index1_r, index2_r, index3_r, index4_r,
																q1coef_r, q2coef_r, q3coef_r, q4coef_r,
																x1_r, x2_r ,x3_r ,x4_r,
																y1_r, y2_r, y3_r, y4_r,
																q1_r, q2_r, q3_r, q4_r, image_point_u_r);
	kernels::interpolateVelocityToHybridNodeY<<<grid,block>>>(u_r, ustar_r, hybridTagsUV_r,
																B.x_r, B.y_r, B.vB_r, yv_r, xv_r,
																body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
																B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny,
																index1_r, index2_r, index3_r, index4_r,
																q1coef_r, q2coef_r, q3coef_r, q4coef_r,
																x1_r, x2_r ,x3_r ,x4_r,
																y1_r, y2_r, y3_r, y4_r,
																q1_r, q2_r, q3_r, q4_r, image_point_u_r);
}

void luo_iter::intermediate_velocity_size_lhs()
{
	cusp::blas::fill(count, 0);

	dim3 grid( 1, 1);
	dim3 block(1, 1);
	kernels::size_LHS1<<<grid,block>>>(hybridTagsUV_r, count_r, B.startI_r, B.startJ_r, B.numCellsXHost, B.numCellsYHost, nx, ny);

	thrust::device_vector<int>::iterator iter = thrust::max_element(count.begin(),count.end());
	unsigned int position = iter - count.begin();
	double max_val = *iter;

	LHS1.resize(numUV,numUV, (nx-1)*ny*5 - 2*(nx-1)-ny*2 + nx*(ny-1)*5 - 2*nx - 2*(ny-1) + max_val);

	LHS1_row_r = thrust::raw_pointer_cast( &(LHS1.row_indices[0]) );
	LHS1_col_r = thrust::raw_pointer_cast( &(LHS1.column_indices[0]) );
	LHS1_val_r = thrust::raw_pointer_cast( &(LHS1.values[0]) );
}

void luo_iter::intermediate_velocity_calculate_lhs()
{
	const int blocksize = 256;
	dim3 gridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 blockU(blocksize, 1);
	dim3 gridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 blockV(blocksize, 1);

	kernels::LHS1_mid_iter_X<<<gridU,blockU>>>(LHS1_row_r, LHS1_col_r, LHS1_val_r, dx_r, dy_r, dt, nu, nx, ny,
													hybridTagsUV_r, ghostTagsUV_r, ns_rhs_r, interp_rhs_r, count_r,
													index1_r, index2_r, index3_r, index4_r,
													xu_r, yu_r, alpha_r, B.uB_r,
													q1coef_r, q2coef_r, q3coef_r, q4coef_r,
													q1_r, q2_r, q3_r, q4_r
													);

	kernels::LHS1_mid_iter_Y<<<gridV,blockV>>>(LHS1_row_r, LHS1_col_r, LHS1_val_r, dx_r, dy_r, dt, nu, nx, ny,
													hybridTagsUV_r, ghostTagsUV_r, ns_rhs_r, interp_rhs_r, count_r,
													index1_r, index2_r, index3_r, index4_r,
													xv_r, yv_r, alpha_r, B.vB_r,
													q1coef_r, q2coef_r, q3coef_r, q4coef_r,
													q1_r, q2_r, q3_r, q4_r
													);

	kernels::LHS_BC_X<<<gridU,blockU>>>(LHS1_row_r, LHS1_col_r, LHS1_val_r, dx_r, dy_r, dt, nu, nx, ny);
	kernels::LHS_BC_Y<<<gridV,blockV>>>(LHS1_row_r, LHS1_col_r, LHS1_val_r, dx_r, dy_r, dt, nu, nx, ny);
}

void luo_iter::intermediate_velocity_update_rhs()
{
	const int blocksize = 256;
	dim3 gridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	dim3 gridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);

	kernels::update_rhs1_x<<<gridU, block>>>(rhs1_r, ns_rhs_r, interp_rhs_r, nx, ny);
	kernels::update_rhs1_y<<<gridU, block>>>(rhs1_r, ns_rhs_r, interp_rhs_r, nx, ny);
}

//testing
void luo_iter::generateRHS1()
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

	//calculate lhs
	intermediate_velocity_calculate_lhs();

	//sort
	LHS1.sort_by_row_and_column();

	//update preconditioner
	if (timeStep == 0)
		PC.generate1(LHS1, (*paramDB)["velocitySolve"]["preconditioner"].get<preconditionerType>());
	else
		PC.update1(LHS1);

	intermediate_velocity_update_rhs();
}
#include <solvers/NavierStokes/luoIBM/kernels/weight.h>//weighting function
void luo_iter::weightUhat()
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

void luo_iter::preRHS1Interpolation()
{
	logger.startTimer("RHS1 Interpolation");

	const int blocksize = 256;

	dim3 grid( int( (B.numCellsXHost*B.numCellsYHost-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	kernels::interpolateVelocityToGhostNodeX<<<grid,block>>>(u_r, true, ghostTagsUV_r, B.x_r, B.y_r, B.uB_r, yu_r, xu_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny,
													index1_r, index2_r, index3_r, index4_r,
													q1coef_r, q2coef_r, q3coef_r, q4coef_r,
													x1_r, x2_r, x3_r, x4_r,
													y1_r, y2_r, y3_r, y4_r,
													q1_r, q2_r, q3_r, q4_r, image_point_u_r);
	kernels::interpolateVelocityToGhostNodeY<<<grid,block>>>(u_r, true, ghostTagsUV_r, B.x_r, B.y_r, B.vB_r, yv_r, xv_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny,
													index1_r, index2_r, index3_r, index4_r,
													q1coef_r, q2coef_r, q3coef_r, q4coef_r,
													x1_r, x2_r, x3_r, x4_r,
													y1_r, y2_r, y3_r, y4_r,
													q1_r, q2_r, q3_r, q4_r, image_point_u_r);
	dim3 grid_inside( int( (numUV-0.5)/blocksize ) +1, 1);
	dim3 block_inside(blocksize, 1);
	kernels::setInsideVelocity<<<grid_inside,block_inside>>>(ghostTagsUV_r, u_r, B.uB_r, B.vB_r, nx, ny);
	//interpolate velocity to hybrid node
	//zero ustar
	kernels::interpolateVelocityToHybridNodeX<<<grid,block>>>(u_r, ustar_r, hybridTagsUV_r,
																B.x_r, B.y_r, B.uB_r, yu_r, xu_r,
																body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
																B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny,
																index1_r, index2_r, index3_r, index4_r,
																q1coef_r, q2coef_r, q3coef_r, q4coef_r,
																x1_r, x2_r ,x3_r ,x4_r,
																y1_r, y2_r, y3_r, y4_r,
																q1_r, q2_r, q3_r, q4_r, image_point_u_r);
	kernels::interpolateVelocityToHybridNodeY<<<grid,block>>>(u_r, ustar_r, hybridTagsUV_r,
																B.x_r, B.y_r, B.vB_r, yv_r, xv_r,
																body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
																B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny,
																index1_r, index2_r, index3_r, index4_r,
																q1coef_r, q2coef_r, q3coef_r, q4coef_r,
																x1_r, x2_r ,x3_r ,x4_r,
																y1_r, y2_r, y3_r, y4_r,
																q1_r, q2_r, q3_r, q4_r, image_point_u_r);
	//testInterpX();
	//testInterpY();
	logger.stopTimer("RHS1 Interpolation");
}
