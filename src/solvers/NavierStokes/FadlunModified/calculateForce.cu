/***************************************************************************//**
 * \file
 * \author Anush Krishnan (anush@bu.edu),
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke kernels that will calculate the force on the immersed body
 */
#include <solvers/NavierStokes/fadlunModified.h>

#include <solvers/NavierStokes/FadlunModified/kernels/calculateForce.h>

/**
 * \brief Calculates forces acting on an immersed body (on the device).
 *
 * Uses the control volume approach explained by Lai and Peskin (2000).
 * This is a general method that can be used with any immersed boundary method.
 * It uses only the velocity and pressure fields to calculate the forces, and
 * does not involve any body forces on the immersed boundary.
 * Currently works only for one body.
 */
void fadlunModified::calculateForce()
{
	parameterDB  &db = *NavierStokesSolver::paramDB;

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

	kernels::dragLeftRight <<<dimGrid, dimBlock>>> (FxX_r, u_r, pressure_r, nu, dx_r, dy_r,
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	kernels::dragBottomTop <<<dimGrid, dimBlock>>> (FxY_r, u_r, nu, dx_r, dy_r,
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	dim3 dimGridX( int( ( (B.numCellsX[0]+1)*B.numCellsY[0]-0.5 )/blockSize )+1, 1 );
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

	double *FyX_r = thrust::raw_pointer_cast(&FyX[0]),
	     *FyY_r = thrust::raw_pointer_cast(&FyY[0]),
	     *FyU_r = thrust::raw_pointer_cast(&FyU[0]);

	kernels::liftLeftRight <<<dimGrid, dimBlock>>> (FyX_r, u_r, nu, dx_r, dy_r,
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	kernels::liftBottomTop <<<dimGrid, dimBlock>>> (FyY_r, u_r, pressure_r, nu, dx_r, dy_r,
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	dim3 dimGridY( int( ( B.numCellsX[0]*(B.numCellsY[0]+1)-0.5 )/blockSize )+1, 1 );
	kernels::liftUnsteady <<<dimGridY, dimBlock>>> (FyU_r, u_r, uold_r, ghostTagsUV_r, dx_r, dy_r, dt,
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	B.forceY = thrust::reduce(FyX.begin(), FyX.end()) + thrust::reduce(FyY.begin(), FyY.end()) + thrust::reduce(FyU.begin(), FyU.end());

}
