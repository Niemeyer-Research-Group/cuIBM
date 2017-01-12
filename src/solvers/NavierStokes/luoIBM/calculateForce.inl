/***************************************************************************//**
 * \file calculateForce.inl
 * \author Anush Krishnan (anush@bu.edu),
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke kernels that will calculate the force on the immersed body
 */

#include <solvers/NavierStokes/FadlunModified/kernels/calculateForce.h>
#include <solvers/NavierStokes/luoIBM/kernels/calculateForce.h>

/**
 * \brief Calculates forces acting on an immersed body (on the device).
 *
 * Uses the control volume approach explained by Lai and Peskin (2000).
 * This is a general method that can be used with any immersed boundary method.
 * It uses only the velocity and pressure fields to calculate the forces, and
 * does not involve any body forces on the immersed boundary.
 * Currently works only for one body.
 */
void luoIBM::calculateForce()
{
	int  nx = NavierStokesSolver::domInfo->nx,
	     ny = NavierStokesSolver::domInfo->ny;
	parameterDB  &db = *NavierStokesSolver::paramDB;
	double	dt = db["simulation"]["dt"].get<double>(),
			nu = db["flow"]["nu"].get<double>();

	double	*u_r		= thrust::raw_pointer_cast(&NavierStokesSolver::u[0]),
			*uold_r		= thrust::raw_pointer_cast(&NavierStokesSolver::uold[0]),
			*pressure_r	= thrust::raw_pointer_cast(&NavierStokesSolver::pressure[0]),
			*dx		= thrust::raw_pointer_cast(&(NavierStokesSolver::domInfo->dx[0])),
			*dy		= thrust::raw_pointer_cast(&(NavierStokesSolver::domInfo->dy[0]));

	int		*ghostTagsUV_r	= thrust::raw_pointer_cast( &(ghostTagsUV[0]) );

	// Calculating drag
	cusp::array1d<double, cusp::device_memory>
		FxX(B.numCellsY[0]),
		FxY(B.numCellsX[0]+1),
		FxU((B.numCellsX[0]+1)*B.numCellsY[0]);

	double *FxX_r = thrust::raw_pointer_cast(&FxX[0]),
	     *FxY_r = thrust::raw_pointer_cast(&FxY[0]),
	     *FxU_r = thrust::raw_pointer_cast(&FxU[0]);

	const int blockSize = 256;
	dim3 dimGrid( int((B.numCellsX[0]+B.numCellsY[0]+1-0.5)/blockSize)+1, 1 );
	dim3 dimBlock(blockSize, 1);

	kernels::dragLeftRight <<<dimGrid, dimBlock>>> (FxX_r, u_r, pressure_r, nu, dx, dy, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	kernels::dragBottomTop <<<dimGrid, dimBlock>>> (FxY_r, u_r, nu, dx, dy, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	dim3 dimGridX( int( ( (B.numCellsX[0]+1)*B.numCellsY[0]-0.5 )/blockSize )+1, 1 );
	kernels::dragUnsteady <<<dimGridX, dimBlock>>> (FxU_r, u_r, uold_r, ghostTagsUV_r, dx, dy, dt, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	fxx = thrust::reduce(FxX.begin(), FxX.end());
	fxy = thrust::reduce(FxY.begin(), FxY.end());
	fxu = thrust::reduce(FxU.begin(), FxU.end());
	B.forceX[0] =  fxx + fxy + fxu;

	
	// Calculating lift
	cusp::array1d<double, cusp::device_memory>
		FyX(B.numCellsY[0]+1),
		FyY(B.numCellsX[0]),
		FyU((B.numCellsX[0]+1)*B.numCellsY[0]);

	double	*FyX_r = thrust::raw_pointer_cast(&FyX[0]),
			*FyY_r = thrust::raw_pointer_cast(&FyY[0]),
			*FyU_r = thrust::raw_pointer_cast(&FyU[0]);

	kernels::liftLeftRight <<<dimGrid, dimBlock>>> (FyX_r, u_r, nu, dx, dy, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	kernels::liftBottomTop <<<dimGrid, dimBlock>>> (FyY_r, u_r, pressure_r, nu, dx, dy, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	dim3 dimGridY( int( ( B.numCellsX[0]*(B.numCellsY[0]+1)-0.5 )/blockSize )+1, 1 );
	kernels::liftUnsteady <<<dimGridY, dimBlock>>> (FyU_r, u_r, uold_r, ghostTagsUV_r, dx, dy, dt, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	B.forceY[0] = thrust::reduce(FyX.begin(), FyX.end()) + thrust::reduce(FyY.begin(), FyY.end()) + thrust::reduce(FyU.begin(), FyU.end());
}

void luoIBM::luoForce()
{
	double	*force_pressure_r = thrust::raw_pointer_cast ( &(B.force_pressure[0])),
			*force_dudn_r = thrust::raw_pointer_cast ( &(B.force_dudn[0]) ),
			*force_dvdn_r = thrust::raw_pointer_cast ( &(B.force_dvdn[0]) ),
			*force_x_r	= thrust::raw_pointer_cast ( &(B.force_x[0]) ),
			*force_y_r	= thrust::raw_pointer_cast ( &(B.force_y[0]) ),
			*u_r		= thrust::raw_pointer_cast ( &(u[0]) ),
			*bx_r		= thrust::raw_pointer_cast ( &(B.x[0]) ),
			*by_r		= thrust::raw_pointer_cast ( &(B.y[0]) ),
			*uB_r		= thrust::raw_pointer_cast ( &(B.uB[0]) ),
			*vB_r		= thrust::raw_pointer_cast ( &(B.vB[0]) ),
			*body_intercept_p_r = thrust::raw_pointer_cast ( &(body_intercept_p[0]) ),
			*body_intercept_p_x_r = thrust::raw_pointer_cast ( &(body_intercept_p_x[0]) ),
			*body_intercept_p_y_r = thrust::raw_pointer_cast ( &(body_intercept_p_y[0]) ),
			*xv_r = thrust::raw_pointer_cast ( &(domInfo->xv[0]) ),
			*xu_r = thrust::raw_pointer_cast ( &(domInfo->xu[0]) ),
			*yu_r = thrust::raw_pointer_cast ( &(domInfo->yu[0]) ),
			*yv_r = thrust::raw_pointer_cast ( &(domInfo->yv[0]) );
			
		
	int nx = domInfo ->nx,
		ny = domInfo ->ny,
		width_i = B.numCellsXHost, //flag this value is only moved to the host once (in B.initialise) if the body is moving too much this could break 
		height_j = B.numCellsYHost;  //this is done because we need the value on the host to calculate the grid size, but copying it to the host every TS is expensive
	int *i_start_r = thrust::raw_pointer_cast ( &(B.startI[0]) ),
		*j_start_r = thrust::raw_pointer_cast ( &(B.startJ[0]) ),
		*ghostTagsP_r=thrust::raw_pointer_cast ( &(ghostTagsP[0]) );
	
	double nu = (*paramDB)["flow"]["nu"].get<double>();
	const int blocksize = 256;
	dim3 grid( int( (B.totalPoints-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	
	kernels::force_pressure<<<grid,block>>>(force_pressure_r, body_intercept_p_r,
												body_intercept_p_x_r, body_intercept_p_y_r,
												bx_r, by_r, xv_r, yu_r, ghostTagsP_r,
												i_start_r, j_start_r, width_i, height_j, B.totalPoints, nx, ny, B.midX, B.midY);
	kernels::force_velocity_x<<<grid,block>>>(force_dudn_r, uB_r, u_r,
												bx_r, by_r, xu_r, yu_r,
												i_start_r, j_start_r, width_i, height_j, B.totalPoints, nx, ny, B.midX, B.midY, domInfo->mid_h);
	kernels::force_velocity_x<<<grid,block>>>(force_dudn_r, vB_r, u_r,
												bx_r, by_r, xv_r, yv_r,
												i_start_r, j_start_r, width_i, height_j, B.totalPoints, nx, ny, B.midX, B.midY, domInfo->mid_h);
	kernels::force<<<grid,block>>>(force_x_r, force_y_r, force_pressure_r, force_dudn_r, force_dvdn_r,
												bx_r, by_r,
												B.totalPoints, B.midX, B.midY, nu);
	
	//testForce_p();
	//testForce_dudn();
	
	B.forceX = thrust::reduce(B.force_x.begin(), B.force_x.end());
	B.forceY = thrust::reduce(B.force_y.begin(), B.force_y.end());
}