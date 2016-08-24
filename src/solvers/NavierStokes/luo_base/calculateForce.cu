/***************************************************************************//**
 * \file
 * \author Anush Krishnan (anush@bu.edu),
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke kernels that will calculate the force on the immersed body
 */
#include <solvers/NavierStokes/luo_base.h>

#include <solvers/NavierStokes/FadlunModified/kernels/calculateForce.h>
#include <solvers/NavierStokes/luo_base/kernels/calculateForce.h>

/**
 * \brief Calculates forces acting on an immersed body (on the device).
 *
 * Uses the control volume approach explained by Lai and Peskin (2000).
 * This is a general method that can be used with any immersed boundary method.
 * It uses only the velocity and pressure fields to calculate the forces, and
 * does not involve any body forces on the immersed boundary.
 * Currently works only for one body.
 */
void luo_base::calculateForce()
{
	// Calculating drag
	cusp::array1d<double, cusp::device_memory>
		FxX(B.numCellsY[0]),
		FxY(B.numCellsX[0]+1),
		FxU((B.numCellsX[0]+1)*B.numCellsY[0]);

	double	*FxX_r = thrust::raw_pointer_cast(&FxX[0]),
			*FxY_r = thrust::raw_pointer_cast(&FxY[0]),
			*FxU_r = thrust::raw_pointer_cast(&FxU[0]);

	const int blockSize = 256;

	dim3 dimGrid( int((B.numCellsX[0]+B.numCellsY[0]+1-0.5)/blockSize)+1, 1 );
	dim3 dimBlock(blockSize, 1);
	dim3 dimGridX( int( ( (B.numCellsX[0]+1)*B.numCellsY[0]-0.5 )/blockSize )+1, 1 );
	dim3 dimGridY( int( ( B.numCellsX[0]*(B.numCellsY[0]+1)-0.5 )/blockSize )+1, 1 );

	kernels::dragLeftRight <<<dimGrid, dimBlock>>> (FxX_r, u_r, pressure_r, nu, dx_r, dy_r,
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	kernels::dragBottomTop <<<dimGrid, dimBlock>>> (FxY_r, u_r, nu, dx_r, dy_r,
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	kernels::dragUnsteady <<<dimGridX, dimBlock>>> (FxU_r, u_r, uold_r, ghostTagsUV_r, dx_r, dy_r, dt,
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	fxx = thrust::reduce(FxX.begin(), FxX.end());
	fxy = thrust::reduce(FxY.begin(), FxY.end());
	fxu = thrust::reduce(FxU.begin(), FxU.end());
	B.forceX =  fxx + fxy + fxu;

	// Calculating lift
	cusp::array1d<double, cusp::device_memory>
		FyX(B.numCellsY[0]+1),
		FyY(B.numCellsX[0]),
		FyU((B.numCellsX[0]+1)*B.numCellsY[0]);

	double	*FyX_r = thrust::raw_pointer_cast(&FyX[0]),
			*FyY_r = thrust::raw_pointer_cast(&FyY[0]),
			*FyU_r = thrust::raw_pointer_cast(&FyU[0]);

	kernels::liftLeftRight <<<dimGrid, dimBlock>>> (FyX_r, u_r, nu, dx_r, dy_r,
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	kernels::liftBottomTop <<<dimGrid, dimBlock>>> (FyY_r, u_r, pressure_r, nu, dx_r, dy_r,
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);


	kernels::liftUnsteady <<<dimGridY, dimBlock>>> (FyU_r, u_r, uold_r, ghostTagsUV_r, dx_r, dy_r, dt,
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	B.forceY = thrust::reduce(FyX.begin(), FyX.end()) + thrust::reduce(FyY.begin(), FyY.end()) + thrust::reduce(FyU.begin(), FyU.end());
}

void luo_base::luoForce()
{
	const int blocksize = 256;
	dim3 grid( int( (B.totalPoints-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);

	kernels::force_pressure<<<grid,block>>>(B.force_pressure_r, body_intercept_p_r,
												body_intercept_p_x_r, body_intercept_p_y_r,
												B.x_r, B.y_r, xv_r, yu_r, ghostTagsP_r,
												B.startI_r, B.startI_r, B.numCellsXHost, B.numCellsYHost, B.totalPoints, nx, ny, B.midX, B.midY);
	kernels::force_velocity_x<<<grid,block>>>(B.force_dudn_r, B.uB_r, u_r,
												B.x_r, B.y_r, xu_r, yu_r,
												B.startI_r, B.startJ_r, B.numCellsXHost, B.numCellsYHost, B.totalPoints, nx, ny, B.midX, B.midY, domInfo->mid_h);
	kernels::force_velocity_y<<<grid,block>>>(B.force_dvdn_r, B.vB_r, u_r,
												B.x_r, B.y_r, xv_r, yv_r,
												B.startI_r, B.startJ_r, B.numCellsXHost, B.numCellsYHost, B.totalPoints, nx, ny, B.midX, B.midY, domInfo->mid_h);
	kernels::force<<<grid,block>>>(B.force_x_r, B.force_y_r, B.force_pressure_r, B.force_dudn_r, B.force_dvdn_r,
												B.x_r, B.y_r,
												B.totalPoints, B.midX, B.midY, nu);

	//testForce_p();
	//testForce_dudn();

	B.forceX = thrust::reduce(B.force_x.begin(), B.force_x.end());
	B.forceY = thrust::reduce(B.force_y.begin(), B.force_y.end());
}
