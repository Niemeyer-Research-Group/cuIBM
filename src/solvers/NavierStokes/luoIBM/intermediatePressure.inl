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
#include <solvers/NavierStokes/NavierStokes/kernels/intermediatePressure.h>

void luoIBM::generateRHS2()
{
	NavierStokesSolver::logger.startTimer("RHS2");

	int nx = NavierStokesSolver::domInfo ->nx,
		ny = NavierStokesSolver::domInfo ->ny;
	
	//variables related to the grid
	double	*uB_r		= thrust::raw_pointer_cast ( &(B.uB[0]) ),
			*vB_r		= thrust::raw_pointer_cast ( &(B.vB[0]) ),	
			*uB0_r		= thrust::raw_pointer_cast ( &(B.uBk[0]) ),
			*vB0_r		= thrust::raw_pointer_cast ( &(B.vBk[0]) ),
			*bx_r		= thrust::raw_pointer_cast ( &(B.x[0]) ),
			*by_r		= thrust::raw_pointer_cast ( &(B.y[0]) ),
			*dx_r	= thrust::raw_pointer_cast( &(domInfo->dx[0]) ),
			*dy_r	= thrust::raw_pointer_cast( &(domInfo->dy[0]) ),
			*yu_r		= thrust::raw_pointer_cast ( &(domInfo->yu[0]) ),
			*xv_r		= thrust::raw_pointer_cast ( &(domInfo->xv[0]) ),
			*xu_r		= thrust::raw_pointer_cast ( &(domInfo->xu[0]) ),
			*yv_r		= thrust::raw_pointer_cast ( &(domInfo->yv[0]) ),
			*ym_r	= thrust::raw_pointer_cast( &(bc[YMINUS][0]) ),
			*yp_r	= thrust::raw_pointer_cast( &(bc[YPLUS][0]) ),
			*xm_r	= thrust::raw_pointer_cast( &(bc[XMINUS][0]) ),
			*xp_r	= thrust::raw_pointer_cast( &(bc[XPLUS][0]) );

	//other variables
	double	*alpha_r= thrust::raw_pointer_cast( &(alpha[0]) ),
			*detA_r = thrust::raw_pointer_cast( &(detA[0]) ),
			*rhs2_r	= thrust::raw_pointer_cast( &(rhs2[0]) ),
			*uhat_r	= thrust::raw_pointer_cast( &(uhat[0]) ),
			*u_r 		= thrust::raw_pointer_cast ( &(u[0]) );
	
	//variables used for interpolation
	double 	*body_intercept_p_x_r = thrust::raw_pointer_cast( &(body_intercept_p_x[0]) ),
			*body_intercept_p_y_r = thrust::raw_pointer_cast( &(body_intercept_p_y[0]) ),
			*image_point_p_x_r = thrust::raw_pointer_cast( &(image_point_p_x[0]) ),
			*image_point_p_y_r = thrust::raw_pointer_cast( &(image_point_p_y[0]) ),
			*body_intercept_x_r = thrust::raw_pointer_cast( &(body_intercept_x[0]) ),
			*body_intercept_y_r = thrust::raw_pointer_cast( &(body_intercept_y[0]) ),
			*image_point_x_r = thrust::raw_pointer_cast( &(image_point_x[0]) ),
			*image_point_y_r = thrust::raw_pointer_cast( &(image_point_y[0]) ),
			*x1_p_r = thrust::raw_pointer_cast ( &(x1_p[0]) ),
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
			*b11_r	= thrust::raw_pointer_cast( &(b11[0]) ),
			*b12_r	= thrust::raw_pointer_cast( &(b12[0]) ),
			*b13_r	= thrust::raw_pointer_cast( &(b13[0]) ),
			*b14_r	= thrust::raw_pointer_cast( &(b14[0]) ),
			*b21_r	= thrust::raw_pointer_cast( &(b21[0]) ),
			*b22_r	= thrust::raw_pointer_cast( &(b22[0]) ),
			*b23_r	= thrust::raw_pointer_cast( &(b23[0]) ),
			*b24_r	= thrust::raw_pointer_cast( &(b24[0]) ),
			*b31_r	= thrust::raw_pointer_cast( &(b31[0]) ),
			*b32_r	= thrust::raw_pointer_cast( &(b32[0]) ),
			*b33_r	= thrust::raw_pointer_cast( &(b33[0]) ),
			*b34_r	= thrust::raw_pointer_cast( &(b34[0]) ),
			*b41_r	= thrust::raw_pointer_cast( &(b41[0]) ),
			*b42_r	= thrust::raw_pointer_cast( &(b42[0]) ),
			*b43_r	= thrust::raw_pointer_cast( &(b43[0]) ),
			*b44_r	= thrust::raw_pointer_cast( &(b44[0]) ),
			*interpCoef_r	= thrust::raw_pointer_cast( &(interpCoef[0]) ),
			*stencilCoef_r	= thrust::raw_pointer_cast( &(stencilCoef[0]) );
			
	//variables related to the uhat interpolation
	double	*dudt_r	= thrust::raw_pointer_cast ( &(dudt[0]) ),
			*dvdt_r	= thrust::raw_pointer_cast ( &(dvdt[0]) ),
			*ududx_r	= thrust::raw_pointer_cast ( &(ududx[0]) ),
			*vdudy_r	= thrust::raw_pointer_cast ( &(vdudy[0]) ),
			*udvdx_r	= thrust::raw_pointer_cast ( &(udvdx[0]) ),
			*vdvdy_r	= thrust::raw_pointer_cast ( &(vdvdy[0]) ),
			*x1_r = thrust::raw_pointer_cast ( &(x1[0]) ),
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
			*image_point_u_r = thrust::raw_pointer_cast( &(image_point_u[0]) );

	int 	*ghostTagsP_r		= thrust::raw_pointer_cast ( &(ghostTagsP[0]) ),
			*hybridTagsP_r		= thrust::raw_pointer_cast ( &(hybridTagsP[0]) ),
			*ghostTagsUV_r		= thrust::raw_pointer_cast ( &(ghostTagsUV[0]) ),
			*hybridTagsUV_r		= thrust::raw_pointer_cast ( &(hybridTagsUV[0]) ),
			*index1_r = thrust::raw_pointer_cast ( &(index1[0]) ),
			*index2_r = thrust::raw_pointer_cast ( &(index2[0]) ),
			*index3_r = thrust::raw_pointer_cast ( &(index3[0]) ),
			*index4_r = thrust::raw_pointer_cast ( &(index4[0]) );
			
	bool	*q1flag_r = thrust::raw_pointer_cast ( &(q1flag[0]) ),
			*q2flag_r = thrust::raw_pointer_cast ( &(q2flag[0]) ),
			*q3flag_r = thrust::raw_pointer_cast ( &(q3flag[0]) ),
			*q4flag_r = thrust::raw_pointer_cast ( &(q4flag[0]) );
	
	double	dt = (*paramDB)["simulation"]["dt"].get<double>();

	int *i_start_r = thrust::raw_pointer_cast ( &(B.startI[0]) ),
		*j_start_r = thrust::raw_pointer_cast ( &(B.startJ[0]) ),
		width = B.numCellsXHost, //flag this value is only moved to the host once (in B.initialise) if the body is moving too much this could break 
		height=B.numCellsYHost;
	
	const int blocksize = 256;
	dim3 gridsmall( int( (width*height-0.5)/blocksize ) +1, 1);
	dim3 gridbig( int( (nx*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	//rhsnormal
	kernels::intermediatePressureNoBody<<<gridbig,block>>>(rhs2_r, uhat_r, ym_r, yp_r, xm_r, xp_r, dx_r, dy_r, nx, ny);
	//rhsinterp
	kernels::intermediatePressure_luo<<<gridbig,block>>>(rhs2_r, detA_r, hybridTagsP_r, alpha_r, stencilCoef_r,
															xv_r, yu_r,
															b11_r, b12_r, b13_r, b14_r, b21_r, b22_r, b23_r, b24_r,
															b31_r, b32_r, b33_r, b34_r, b41_r, b42_r, b43_r, b44_r,
															q1_p_r, q2_p_r, q3_p_r, q4_p_r,
															q1flag_r, q2flag_r, q3flag_r, q4flag_r,
															nx, ny);
	NavierStokesSolver::logger.stopTimer("RHS2");
}

void luoIBM::generateLHS2()
{
	NavierStokesSolver::logger.startTimer("LHS2");
	int nx = domInfo ->nx,
		ny = domInfo ->ny;
	parameterDB  &db = *NavierStokesSolver::paramDB;
	
	double	dt 		= db["simulation"]["dt"].get<double>(),
			*dx_r	= thrust::raw_pointer_cast( &(domInfo->dx[0]) ),
			*dy_r	= thrust::raw_pointer_cast( &(domInfo->dy[0]) ),
			*val_r	= thrust::raw_pointer_cast( &(LHS2.values[0]) ),
			*alpha_r= thrust::raw_pointer_cast( &(alpha[0]) ),
			*detA_r = thrust::raw_pointer_cast( &(detA[0]) ),
			*xv_r	= thrust::raw_pointer_cast( &(domInfo->xv[0]) ),
			*yu_r	= thrust::raw_pointer_cast( &(domInfo->yu[0]) ),
			*b11_r	= thrust::raw_pointer_cast( &(b11[0]) ),
			*b12_r	= thrust::raw_pointer_cast( &(b12[0]) ),
			*b13_r	= thrust::raw_pointer_cast( &(b13[0]) ),
			*b14_r	= thrust::raw_pointer_cast( &(b14[0]) ),
			*b21_r	= thrust::raw_pointer_cast( &(b21[0]) ),
			*b22_r	= thrust::raw_pointer_cast( &(b22[0]) ),
			*b23_r	= thrust::raw_pointer_cast( &(b23[0]) ),
			*b24_r	= thrust::raw_pointer_cast( &(b24[0]) ),
			*b31_r	= thrust::raw_pointer_cast( &(b31[0]) ),
			*b32_r	= thrust::raw_pointer_cast( &(b32[0]) ),
			*b33_r	= thrust::raw_pointer_cast( &(b33[0]) ),
			*b34_r	= thrust::raw_pointer_cast( &(b34[0]) ),
			*b41_r	= thrust::raw_pointer_cast( &(b41[0]) ),
			*b42_r	= thrust::raw_pointer_cast( &(b42[0]) ),
			*b43_r	= thrust::raw_pointer_cast( &(b43[0]) ),
			*b44_r	= thrust::raw_pointer_cast( &(b44[0]) ),
			*interpCoef_r	= thrust::raw_pointer_cast( &(interpCoef[0]) ),
			*stencilCoef_r	= thrust::raw_pointer_cast( &(stencilCoef[0]) );

	int		*row_r	= thrust::raw_pointer_cast( &(LHS2.row_indices[0]) ),
			*col_r	= thrust::raw_pointer_cast( &(LHS2.column_indices[0]) ),
			*ghostTagsP_r		= thrust::raw_pointer_cast ( &(ghostTagsP[0]) ),
			*hybridTagsP_r		= thrust::raw_pointer_cast ( &(hybridTagsP[0]) ),
			*index1_r = thrust::raw_pointer_cast ( &(index1[0]) ),
			*index2_r = thrust::raw_pointer_cast ( &(index2[0]) ),
			*index3_r = thrust::raw_pointer_cast ( &(index3[0]) ),
			*index4_r = thrust::raw_pointer_cast ( &(index4[0]) ),
			*countD_r	= thrust::raw_pointer_cast( &(countD[0]) );
	
	bool	*q1flag_r = thrust::raw_pointer_cast ( &(q1flag[0]) ),
			*q2flag_r = thrust::raw_pointer_cast ( &(q2flag[0]) ),
			*q3flag_r = thrust::raw_pointer_cast ( &(q3flag[0]) ),
			*q4flag_r = thrust::raw_pointer_cast ( &(q4flag[0]) );
	
	const int blocksize = 256;
	dim3 grid( int( (nx*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	//set LHS to 0
	cusp::blas::fill(LHS2.row_indices,0);
	cusp::blas::fill(LHS2.column_indices,0);
	cusp::blas::fill(LHS2.values,0);
	kernels::LHS2_mid_luo<<<grid,block>>>(row_r, col_r, val_r, dx_r, dy_r, nx, ny, dt, countD_r, stencilCoef_r, interpCoef_r,
											detA_r, hybridTagsP_r, ghostTagsP_r, alpha_r,
											xv_r, yu_r,
											b11_r, b12_r, b13_r, b14_r, b21_r, b22_r, b23_r, b24_r,
											b31_r, b32_r, b33_r, b34_r, b41_r, b42_r, b43_r, b44_r,
											q1flag_r, q2flag_r, q3flag_r, q4flag_r,
											index1_r,index2_r,index3_r,index4_r);
	kernels::LHS2_BC<<<grid,block>>>(row_r, col_r, val_r, dx_r, dy_r, nx,ny,dt);
	logger.stopTimer("LHS2");
	
	logger.startTimer("Preconditioner");
	//if (timeStep == 0 || iterationCount2 > 100)
	{
		//seems like we need to call this every turn due to the pressure interpolation?
		PC.generate(LHS1,LHS2, (*paramDB)["velocitySolve"]["preconditioner"].get<preconditionerType>(), (*paramDB)["PoissonSolve"]["preconditioner"].get<preconditionerType>());
	}
	/*else if (iterationCount2 > 100)
	{
		//PC.update(LHS1, LHS2);
	}*/
	logger.stopTimer("Preconditioner");
}

void luoIBM::sizeLHS2()
{
	NavierStokesSolver::logger.startTimer("LHS2");
	int nx = domInfo ->nx,
		ny = domInfo ->ny;

	int		*hybridTagsP_r		= thrust::raw_pointer_cast ( &(hybridTagsP[0]) ),
			*index1_r = thrust::raw_pointer_cast ( &(index1[0]) ),
			*index2_r = thrust::raw_pointer_cast ( &(index2[0]) ),
			*index3_r = thrust::raw_pointer_cast ( &(index3[0]) ),
			*index4_r = thrust::raw_pointer_cast ( &(index4[0]) ),
			*countD_r	= thrust::raw_pointer_cast( &(countD[0]) );
	
	int *i_start_r = thrust::raw_pointer_cast ( &(B.startI[0]) ),
		*j_start_r = thrust::raw_pointer_cast ( &(B.startJ[0]) ),
		width = B.numCellsXHost, 
		height=B.numCellsYHost;

	//calculate indices for extra nodes
	cusp::blas::fill(countD,0);
	const int blocksize = 256;
	dim3 grid(1,1);
	dim3 block(blocksize, 1);

	kernels::hybridPressureNodeCount<<<grid,block>>>(countD_r, index1_r, index2_r, index3_r, index4_r, hybridTagsP_r,
								i_start_r, j_start_r, width, height, nx, ny);
	//find max value
	thrust::device_vector<int>::iterator iter = thrust::max_element(countD.begin(),countD.end());
	unsigned int position = iter - countD.begin();
	int max_val = *iter;

	LHS2.resize(nx*ny, nx*ny, nx*ny*5-nx*2-ny*2+max_val);
	logger.stopTimer("LHS2");
}

void luoIBM::interpPGN()
{
	logger.startTimer("P interp");
	double	*u_r 		= thrust::raw_pointer_cast ( &(u[0]) ),
			*pressure_r 	= thrust::raw_pointer_cast ( &(pressure[0]) ),
			*yu_r		= thrust::raw_pointer_cast ( &(domInfo->yu[0]) ),
			*xu_r		= thrust::raw_pointer_cast ( &(domInfo->xu[0]) ),
			*yv_r		= thrust::raw_pointer_cast ( &(domInfo->yv[0]) ),
			*xv_r		= thrust::raw_pointer_cast ( &(domInfo->xv[0]) ),
			*body_intercept_p_x_r = thrust::raw_pointer_cast( &(body_intercept_p_x[0]) ),
			*body_intercept_p_y_r = thrust::raw_pointer_cast( &(body_intercept_p_y[0]) ),
			*body_intercept_p_r = thrust::raw_pointer_cast( &(body_intercept_p[0]) ),
			*image_point_p_x_r = thrust::raw_pointer_cast( &(image_point_p_x[0]) ),
			*image_point_p_y_r = thrust::raw_pointer_cast( &(image_point_p_y[0]) ),
			*uB0_r		= thrust::raw_pointer_cast ( &(B.uBk[0]) ),
			*vB0_r		= thrust::raw_pointer_cast ( &(B.vBk[0]) ),
			*uB_r		= thrust::raw_pointer_cast ( &(B.uB[0]) ),
			*vB_r		= thrust::raw_pointer_cast ( &(B.vB[0]) ),
			*bx_r		= thrust::raw_pointer_cast ( &(B.x[0]) ),
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
			*a3_r	= thrust::raw_pointer_cast ( &(a3[0]) ),
			*dudt_r	= thrust::raw_pointer_cast ( &(dudt[0]) ),
			*dvdt_r	= thrust::raw_pointer_cast ( &(dvdt[0]) ),
			*ududx_r	= thrust::raw_pointer_cast ( &(ududx[0]) ),
			*vdudy_r	= thrust::raw_pointer_cast ( &(vdudy[0]) ),
			*udvdx_r	= thrust::raw_pointer_cast ( &(udvdx[0]) ),
			*vdvdy_r	= thrust::raw_pointer_cast ( &(vdvdy[0]) );

	int 	*ghostTagsP_r		= thrust::raw_pointer_cast ( &(ghostTagsP[0]) );

	int nx = domInfo ->nx,
		ny = domInfo ->ny,
		width_i = B.numCellsXHost, //flag this value is only moved to the host once (in B.initialise) if the body is moving too much this could break 
		height_j=B.numCellsYHost;  //this is done because we need the value on the host to calculate the grid size, but copying it to the host every TS is expensive
	int *i_start_r = thrust::raw_pointer_cast ( &(B.startI[0]) ),
		*j_start_r = thrust::raw_pointer_cast ( &(B.startJ[0]) );
	double dt = (*paramDB)["simulation"]["dt"].get<double>();

	const int blocksize = 256;
	dim3 grid( int( (width_i*height_j-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	kernels::interpolatePressureToGhostNode<<<grid,block>>>(pressure_r, u_r, ghostTagsP_r, bx_r, by_r,
									uB_r, uB0_r, vB_r, vB0_r, yu_r, yv_r, xu_r, xv_r,
									body_intercept_p_x_r, body_intercept_p_y_r, image_point_p_x_r, image_point_p_y_r,  body_intercept_p_r,
									i_start_r, j_start_r, width_i, nx, ny, dt, B.totalPoints,
									dudt_r,ududx_r,dvdt_r,vdudy_r,udvdx_r,vdvdy_r,
									a0_r, a1_r, a2_r, a3_r,
									x1_p_r, x2_p_r, x3_p_r, x4_p_r, y1_p_r, y2_p_r, y3_p_r, y4_p_r, q1_p_r, q2_p_r, q3_p_r, q4_p_r);
	logger.stopTimer("P interp");
}

void luoIBM::preRHS2()
{
	NavierStokesSolver::logger.startTimer("RHS2");

	int nx = NavierStokesSolver::domInfo ->nx,
		ny = NavierStokesSolver::domInfo ->ny;

	//variables related to the grid
	double	*uB_r		= thrust::raw_pointer_cast ( &(B.uB[0]) ),
			*vB_r		= thrust::raw_pointer_cast ( &(B.vB[0]) ),	
			*uB0_r		= thrust::raw_pointer_cast ( &(B.uBk[0]) ),
			*vB0_r		= thrust::raw_pointer_cast ( &(B.vBk[0]) ),
			*bx_r		= thrust::raw_pointer_cast ( &(B.x[0]) ),
			*by_r		= thrust::raw_pointer_cast ( &(B.y[0]) ),
			*dx_r	= thrust::raw_pointer_cast( &(domInfo->dx[0]) ),
			*dy_r	= thrust::raw_pointer_cast( &(domInfo->dy[0]) ),
			*yu_r		= thrust::raw_pointer_cast ( &(domInfo->yu[0]) ),
			*xv_r		= thrust::raw_pointer_cast ( &(domInfo->xv[0]) ),
			*xu_r		= thrust::raw_pointer_cast ( &(domInfo->xu[0]) ),
			*yv_r		= thrust::raw_pointer_cast ( &(domInfo->yv[0]) ),
			*ym_r	= thrust::raw_pointer_cast( &(bc[YMINUS][0]) ),
			*yp_r	= thrust::raw_pointer_cast( &(bc[YPLUS][0]) ),
			*xm_r	= thrust::raw_pointer_cast( &(bc[XMINUS][0]) ),
			*xp_r	= thrust::raw_pointer_cast( &(bc[XPLUS][0]) );

	//other variables
	double	*alpha_r= thrust::raw_pointer_cast( &(alpha[0]) ),
			*detA_r = thrust::raw_pointer_cast( &(detA[0]) ),
			*rhs2_r	= thrust::raw_pointer_cast( &(rhs2[0]) ),
			*uhat_r	= thrust::raw_pointer_cast( &(uhat[0]) ),
			*u_r 		= thrust::raw_pointer_cast ( &(u[0]) );

	//variables used for interpolation
	double 	*body_intercept_p_x_r = thrust::raw_pointer_cast( &(body_intercept_p_x[0]) ),
			*body_intercept_p_y_r = thrust::raw_pointer_cast( &(body_intercept_p_y[0]) ),
			*image_point_p_x_r = thrust::raw_pointer_cast( &(image_point_p_x[0]) ),
			*image_point_p_y_r = thrust::raw_pointer_cast( &(image_point_p_y[0]) ),
			*body_intercept_x_r = thrust::raw_pointer_cast( &(body_intercept_x[0]) ),
			*body_intercept_y_r = thrust::raw_pointer_cast( &(body_intercept_y[0]) ),
			*image_point_x_r = thrust::raw_pointer_cast( &(image_point_x[0]) ),
			*image_point_y_r = thrust::raw_pointer_cast( &(image_point_y[0]) ),
			*x1_p_r = thrust::raw_pointer_cast ( &(x1_p[0]) ),
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
			*b11_r	= thrust::raw_pointer_cast( &(b11[0]) ),
			*b12_r	= thrust::raw_pointer_cast( &(b12[0]) ),
			*b13_r	= thrust::raw_pointer_cast( &(b13[0]) ),
			*b14_r	= thrust::raw_pointer_cast( &(b14[0]) ),
			*b21_r	= thrust::raw_pointer_cast( &(b21[0]) ),
			*b22_r	= thrust::raw_pointer_cast( &(b22[0]) ),
			*b23_r	= thrust::raw_pointer_cast( &(b23[0]) ),
			*b24_r	= thrust::raw_pointer_cast( &(b24[0]) ),
			*b31_r	= thrust::raw_pointer_cast( &(b31[0]) ),
			*b32_r	= thrust::raw_pointer_cast( &(b32[0]) ),
			*b33_r	= thrust::raw_pointer_cast( &(b33[0]) ),
			*b34_r	= thrust::raw_pointer_cast( &(b34[0]) ),
			*b41_r	= thrust::raw_pointer_cast( &(b41[0]) ),
			*b42_r	= thrust::raw_pointer_cast( &(b42[0]) ),
			*b43_r	= thrust::raw_pointer_cast( &(b43[0]) ),
			*b44_r	= thrust::raw_pointer_cast( &(b44[0]) );

	//variables related to the uhat interpolation
	double	*dudt_r	= thrust::raw_pointer_cast ( &(dudt[0]) ),
			*dvdt_r	= thrust::raw_pointer_cast ( &(dvdt[0]) ),
			*ududx_r	= thrust::raw_pointer_cast ( &(ududx[0]) ),
			*vdudy_r	= thrust::raw_pointer_cast ( &(vdudy[0]) ),
			*udvdx_r	= thrust::raw_pointer_cast ( &(udvdx[0]) ),
			*vdvdy_r	= thrust::raw_pointer_cast ( &(vdvdy[0]) ),
			*x1_r = thrust::raw_pointer_cast ( &(x1[0]) ),
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
			*image_point_u_r = thrust::raw_pointer_cast( &(image_point_u[0]) );

	int 	*ghostTagsP_r		= thrust::raw_pointer_cast ( &(ghostTagsP[0]) ),
			*hybridTagsP_r		= thrust::raw_pointer_cast ( &(hybridTagsP[0]) ),
			*ghostTagsUV_r		= thrust::raw_pointer_cast ( &(ghostTagsUV[0]) ),
			*hybridTagsUV_r		= thrust::raw_pointer_cast ( &(hybridTagsUV[0]) ),
			*index1_r = thrust::raw_pointer_cast ( &(index1[0]) ),
			*index2_r = thrust::raw_pointer_cast ( &(index2[0]) ),
			*index3_r = thrust::raw_pointer_cast ( &(index3[0]) ),
			*index4_r = thrust::raw_pointer_cast ( &(index4[0]) );

	bool	*q1flag_r = thrust::raw_pointer_cast ( &(q1flag[0]) ),
			*q2flag_r = thrust::raw_pointer_cast ( &(q2flag[0]) ),
			*q3flag_r = thrust::raw_pointer_cast ( &(q3flag[0]) ),
			*q4flag_r = thrust::raw_pointer_cast ( &(q4flag[0]) );

	double	dt = (*paramDB)["simulation"]["dt"].get<double>();

	int *i_start_r = thrust::raw_pointer_cast ( &(B.startI[0]) ),
		*j_start_r = thrust::raw_pointer_cast ( &(B.startJ[0]) ),
		width = B.numCellsXHost, //flag this value is only moved to the host once (in B.initialise) if the body is moving too much this could break 
		height=B.numCellsYHost;

	const int blocksize = 256;
	dim3 gridsmall( int( (width*height-0.5)/blocksize ) +1, 1);
	dim3 gridbig( int( (nx*ny-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	//interpolate uhat to the inside of the body
	kernels::interpolateVelocityToGhostNodeX<<<gridsmall,block>>>(uhat_r, ghostTagsUV_r, bx_r, by_r, uB_r, yu_r, xu_r,
																body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
																i_start_r, j_start_r, width, nx, ny,
																x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r,image_point_u_r);
	kernels::interpolateVelocityToGhostNodeY<<<gridsmall,block>>>(uhat_r, ghostTagsUV_r, bx_r, by_r, vB_r, yv_r, xv_r,
																body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
																i_start_r, j_start_r, width, nx, ny,
																x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r,image_point_u_r);
	//get alpha
	kernels::alpha_<<<gridsmall,block>>>(alpha_r, ghostTagsP_r, hybridTagsP_r, yu_r, xv_r, 
											body_intercept_p_x_r, body_intercept_p_y_r, 
											i_start_r, j_start_r, width, nx, ny);
	//get prereqs for interpolation (detA, index, )
	kernels::interpolate_P_HN_setup<<<gridsmall,block>>>(detA_r, hybridTagsP_r, bx_r, by_r,
															uB_r, uB0_r, vB_r, vB0_r,
															yu_r, xv_r,
															body_intercept_p_x_r, body_intercept_p_y_r, image_point_p_x_r, image_point_p_y_r,
															i_start_r, j_start_r, width, nx, ny, dt, B.totalPoints,
															b11_r, b12_r, b13_r, b14_r, b21_r, b22_r, b23_r, b24_r,
															b31_r, b32_r, b33_r, b34_r, b41_r, b42_r, b43_r, b44_r,
															q1_p_r, q2_p_r, q3_p_r, q4_p_r,
															q1flag_r, q2flag_r, q3flag_r, q4flag_r,
															index1_r,index2_r,index3_r,index4_r,
															x1_p_r, x2_p_r, x3_p_r, x4_p_r,
															y1_p_r, y2_p_r, y3_p_r, y4_p_r,
															dudt_r, ududx_r, vdudy_r, dvdt_r, udvdx_r, vdvdy_r);
	NavierStokesSolver::logger.stopTimer("RHS2");
}