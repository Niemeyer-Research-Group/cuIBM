/***************************************************************************//**
 * \file intermediatePressure.inl
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernals that setup the prerequisites to solve the poission equation
 */

#include <solvers/NavierStokes/luoIBM/kernels/intermediatePressure.h>
#include <solvers/NavierStokes/luoIBM/kernels/LHS2.h>
#include <solvers/NavierStokes/NavierStokes/kernels/LHS2.h>
#include <solvers/NavierStokes/luoIBM/kernels/weight.h>//weighting function
#include <solvers/NavierStokes/luoIBM/kernels/intermediateVelocity.h>

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
			*xp_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::bc[XPLUS][0]) );

	const int blocksize = 256;

	dim3 grid( int( (nx*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	
	preRHS2Interpolation();
	
	kernels::intermediatePressure_luo<<<grid,block>>>(rhs2_r, uhat_r, ym_r, yp_r, xm_r, xp_r, dx_r, dy_r, nx, ny);
	NavierStokesSolver::pressure_old = NavierStokesSolver::pressure;//flag not used anywhere
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
			*val_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::LHS2.values[0]) );

	int		*row_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::LHS2.row_indices[0]) ),
			*col_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::LHS2.column_indices[0]) );

	const int blocksize = 256;
	dim3 grid( int( (nx*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	cusp::blas::fill(LHS2.values,0);
	kernels::LHS2_mid_luo<<<grid,block>>>(row_r, col_r, val_r, dx_r, dy_r, nx, ny, dt);
	kernels::LHS2_BC<<<grid,block>>>(row_r, col_r, val_r, dx_r, dy_r, nx,ny,dt);
	logger.stopTimer("LHS2");
}

void luoIBM::weightPressure()
{
	double	*u_r 		= thrust::raw_pointer_cast ( &(u[0]) ),
			*pressure_r 	= thrust::raw_pointer_cast ( &(pressure[0]) ),
			*pressureStar_r	= thrust::raw_pointer_cast ( &(pressureStar[0]) ),
			*yu_r		= thrust::raw_pointer_cast ( &(domInfo->yu[0]) ),
			*xu_r		= thrust::raw_pointer_cast ( &(domInfo->xu[0]) ),
			*yv_r		= thrust::raw_pointer_cast ( &(domInfo->yv[0]) ),
			*xv_r		= thrust::raw_pointer_cast ( &(domInfo->xv[0]) ),
			*body_intercept_p_x_r = thrust::raw_pointer_cast( &(body_intercept_p_x[0]) ),
			*body_intercept_p_y_r = thrust::raw_pointer_cast( &(body_intercept_p_y[0]) ),
			*image_point_p_x_r = thrust::raw_pointer_cast( &(image_point_p_x[0]) ),
			*image_point_p_y_r = thrust::raw_pointer_cast( &(image_point_p_y[0]) ),
			*uB0_r		= thrust::raw_pointer_cast ( &(B.uBk[0]) ),
			*vB0_r		= thrust::raw_pointer_cast ( &(B.vBk[0]) ),
			*uB_r		= thrust::raw_pointer_cast ( &(B.uB[0]) ),
			*vB_r		= thrust::raw_pointer_cast ( &(B.vB[0]) ),
			*bx_r		= thrust::raw_pointer_cast ( &(B.x[0]) ),//not sure if these are on the host or not
			*by_r		= thrust::raw_pointer_cast ( &(B.y[0]) );

	double 	*x1_p_r = thrust::raw_pointer_cast ( &(x1_p[0]) ),
			*x2_p_r = thrust::raw_pointer_cast ( &(x2_p[0]) ),
			*x3_p_r = thrust::raw_pointer_cast ( &(x3_p[0]) ),
			*x4_p_r = thrust::raw_pointer_cast ( &(x4_p[0]) ),
			*y1_p_r = thrust::raw_pointer_cast ( &(y1_p[0]) ),
			*y2_p_r = thrust::raw_pointer_cast ( &(y2_p[0]) ),
			*y3_p_r = thrust::raw_pointer_cast ( &(y3_p[0]) ),
			*y4_p_r = thrust::raw_pointer_cast ( &(y4_p[0]) ),
			*q1_p_r = thrust::raw_pointer_cast ( &(q1_p[0]) ),
			*q2_p_r = thrust::raw_pointer_cast ( &(q2_p[0]) ),
			*q3_p_r = thrust::raw_pointer_cast ( &(q3_p[0]) ),
			*q4_p_r = thrust::raw_pointer_cast ( &(q4_p[0]) ),
			*a0_r	= thrust::raw_pointer_cast ( &(a0[0]) ),
			*a1_r	= thrust::raw_pointer_cast ( &(a1[0]) ),
			*a2_r	= thrust::raw_pointer_cast ( &(a2[0]) ),
			*a3_r	= thrust::raw_pointer_cast ( &(a3[0]) );

	int 	*ghostTagsP_r		= thrust::raw_pointer_cast ( &(ghostTagsP[0]) ),
			*hybridTagsP_r		= thrust::raw_pointer_cast ( &(hybridTagsP[0]) );

	int nx = domInfo ->nx,
		ny = domInfo ->ny,
		totalPoints = B.totalPoints,
		i_start = B.startI[0],
		j_start = B.startJ[0],
		width_i = B.numCellsX[0],
		height_j = B.numCellsY[0],
		i_end = i_start + width_i,
		j_end = j_start + height_j;

	const int blocksize = 256;
	dim3 grid( int( (width_i*height_j-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	kernels::interpolatePressureToHybridNode<<<grid,block>>>(pressure_r, pressureStar_r, u_r, hybridTagsP_r, bx_r, by_r,
									uB_r, uB0_r, vB_r, vB0_r, yu_r, yv_r, xu_r, xv_r,
									body_intercept_p_x_r, body_intercept_p_y_r, image_point_p_x_r, image_point_p_y_r,
									i_start, j_start, i_end, j_end, nx, ny, totalPoints,
									a0_r, a1_r, a2_r, a3_r,
									x1_p_r, x2_p_r, x3_p_r, x4_p_r, y1_p_r, y2_p_r, y3_p_r, y4_p_r, q1_p_r, q2_p_r, q3_p_r, q4_p_r, timeStep);
	kernels::weightP<<<grid,block>>>(pressure_r, pressureStar_r, ghostTagsP_r, hybridTagsP_r, yu_r, xv_r,
									body_intercept_p_x_r, body_intercept_p_y_r, image_point_p_x_r, image_point_p_y_r,
									i_start, j_start, i_end, j_end, nx, ny);
	kernels::interpolatePressureToGhostNode<<<grid,block>>>(pressure_r, u_r, ghostTagsP_r, bx_r, by_r,
									uB_r, uB0_r, vB_r, vB0_r, yu_r, yv_r, xu_r, xv_r,
									body_intercept_p_x_r, body_intercept_p_y_r, image_point_p_x_r, image_point_p_y_r,
									i_start, j_start, i_end, j_end, nx, ny, totalPoints,
									a0_r, a1_r, a2_r, a3_r,
									x1_p_r, x2_p_r, x3_p_r, x4_p_r, y1_p_r, y2_p_r, y3_p_r, y4_p_r, q1_p_r, q2_p_r, q3_p_r, q4_p_r);
	//testInterpP();
}

void luoIBM::preRHS2Interpolation()
{
	double	*u_r 		= thrust::raw_pointer_cast ( &(u[0]) ),
			*uhat_r 	= thrust::raw_pointer_cast ( &(uhat[0]) ),
			*ustar_r	= thrust::raw_pointer_cast ( &(ustar[0]) ),
			*bx_r		= thrust::raw_pointer_cast ( &(B.x[0]) ),
			*by_r		= thrust::raw_pointer_cast ( &(B.y[0]) ),
			*uB_r		= thrust::raw_pointer_cast ( &(B.uB[0]) ),
			*vB_r		= thrust::raw_pointer_cast ( &(B.vB[0]) ),
			*yu_r		= thrust::raw_pointer_cast ( &(domInfo->yu[0]) ),
			*xu_r		= thrust::raw_pointer_cast ( &(domInfo->xu[0]) ),
			*yv_r		= thrust::raw_pointer_cast ( &(domInfo->yv[0]) ),
			*xv_r		= thrust::raw_pointer_cast ( &(domInfo->xv[0]) ),
			*body_intercept_x_r = thrust::raw_pointer_cast( &(body_intercept_x[0]) ),
			*body_intercept_y_r = thrust::raw_pointer_cast( &(body_intercept_y[0]) ),
			*image_point_x_r = thrust::raw_pointer_cast( &(image_point_x[0]) ),
			*image_point_y_r = thrust::raw_pointer_cast( &(image_point_y[0]) );

	double	*x1_r = thrust::raw_pointer_cast ( &(x1[0]) ),
			*x2_r = thrust::raw_pointer_cast ( &(x2[0]) ),
			*x3_r = thrust::raw_pointer_cast ( &(x3[0]) ),
			*x4_r = thrust::raw_pointer_cast ( &(x4[0]) ),
			*y1_r = thrust::raw_pointer_cast ( &(y1[0]) ),
			*y2_r = thrust::raw_pointer_cast ( &(y2[0]) ),
			*y3_r = thrust::raw_pointer_cast ( &(y3[0]) ),
			*y4_r = thrust::raw_pointer_cast ( &(y4[0]) ),
			*q1_r = thrust::raw_pointer_cast ( &(q1[0]) ),
			*q2_r = thrust::raw_pointer_cast ( &(q2[0]) ),
			*q3_r = thrust::raw_pointer_cast ( &(q3[0]) ),
			*q4_r = thrust::raw_pointer_cast ( &(q4[0]) ),
			*ip_u_r = thrust::raw_pointer_cast( &(ip_u[0]) );

	int 	*ghostTagsUV_r		= thrust::raw_pointer_cast ( &(ghostTagsUV[0]) ),
			*hybridTagsUV_r		= thrust::raw_pointer_cast ( &(hybridTagsUV[0]) );

	parameterDB  &db = *paramDB;
	double	nu = db["flow"]["nu"].get<double>();
	double	dt = db["simulation"]["dt"].get<double>();

	int nx = domInfo ->nx,
		ny = domInfo ->ny,
		totalPoints = B.totalPoints,
		i_start = B.startI[0],
		j_start = B.startJ[0],
		width_i = B.numCellsX[0],
		height_j = B.numCellsY[0],
		i_end = i_start + width_i,
		j_end = j_start + height_j;

	const int blocksize = 256;
	dim3 grid( int( (width_i*height_j-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	//interpolate uhat to image point and ghost node
	kernels::interpolateVelocityToGhostNodeX<<<grid,block>>>(uhat_r, ghostTagsUV_r, bx_r, by_r, uB_r, yu_r, xu_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													i_start, j_start, i_end, j_end, nx, ny, totalPoints,
													x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r,ip_u_r);
	kernels::interpolateVelocityToGhostNodeY<<<grid,block>>>(uhat_r, ghostTagsUV_r, bx_r, by_r, vB_r, yv_r, xv_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													i_start, j_start, i_end, j_end, nx, ny, totalPoints,
													x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r,ip_u_r);
}