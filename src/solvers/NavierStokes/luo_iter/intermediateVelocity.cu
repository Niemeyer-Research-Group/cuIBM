/***************************************************************************//**
 * \file
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernels that setup the intermediate velocity solve
 */
#include <solvers/NavierStokes/luo_iter.h>

//kernels
#include <solvers/NavierStokes/NavierStokes/kernels/intermediateVelocity.h> //generateRHS
#include <solvers/NavierStokes/luo_base/kernels/intermediateVelocity.h>
#include <solvers/NavierStokes/luo_base/kernels/biLinearInterpolation.h> //interpolate

void luo_iter::intermediate_velocity_setup()
{
	logger.startTimer("Intermediate Velocity Setup");

	//set correct grid and block size
	const int blocksize = 256;

	dim3 grid( int( ((nx-1)*ny+nx*(ny-1)-0.5)/blocksize ) + 1, 1);
	dim3 block(blocksize, 1);
	dim3 grid_inside( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
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

	//interpolate u and v to image points, then again to ghost nodes
	intermediate_velocity_interpolation_setup();//here

	//size lhs
	intermediate_velocity_size_lhs();

	//calculate lhs
	intermediate_velocity_calculate_lhs();

	//update rhs
	intermediate_velocity_update_rhs();

	logger.stopTimer("Intermediate Velocity Setup");
}

void luo_iter::intermediate_velocity_alpha()
{

}

void luo_iter::intermediate_velocity_interpolation_setup()
{
	const int blocksize = 256;

	dim3 grid( int( (B.numCellsXHost*B.numCellsYHost-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	//interpolate velocity to image point and ghost node
	kernels::interpolateVelocityToGhostNodeX<<<grid,block>>>(u_r, ghostTagsUV_r, B.x_r, B.y_r, B.uB_r, yu_r, xu_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny,
													x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r,image_point_u_r);
	kernels::interpolateVelocityToGhostNodeY<<<grid,block>>>(u_r, ghostTagsUV_r, B.x_r, B.y_r, B.vB_r, yv_r, xv_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny,
													x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r,image_point_u_r);

	cusp::blas::fill(ustar,0.0);//flag not needed?
	kernels::interpolateVelocityToHybridNodeX<<<grid,block>>>(u_r, ustar_r, hybridTagsUV_r, B.x_r, B.y_r, B.uB_r, yu_r, xu_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny,
													x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r, image_point_u_r);
	kernels::interpolateVelocityToHybridNodeY<<<grid,block>>>(u_r, ustar_r, hybridTagsUV_r, B.x_r, B.y_r, B.vB_r, yv_r, xv_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													B.startI_r, B.startJ_r, B.numCellsXHost, nx, ny,
													x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r, image_point_u_r);
	//testInterpX();
	//testInterpY();
}


void luo_iter::intermediate_velocity_size_lhs()
{

}

void luo_iter::intermediate_velocity_calculate_lhs()
{

}

void luo_iter::intermediate_velocity_update_rhs()
{

}
