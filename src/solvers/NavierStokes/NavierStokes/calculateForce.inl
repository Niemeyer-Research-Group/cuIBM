/***************************************************************************//**
 * \file calculateForce.inl
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke kernels that will calculate the force on the immersed body
 */

#include <solvers/NavierStokes/kernels/calculateForce.h>

void NavierStokesSolver::calculateForceFadlun()
{
	double	*uold_r	= thrust::raw_pointer_cast( &(uold[0]) ),
			*Lnew_r	= thrust::raw_pointer_cast( &(Lnew[0]) ),
			*L_r	= thrust::raw_pointer_cast( &(L[0]) ),
			*u_r	= thrust::raw_pointer_cast( &(u[0]) ),
			*force_r= thrust::raw_pointer_cast( &(force[0]) ),
			*N_r	= thrust::raw_pointer_cast( &(N[0]) ),
			*Nold_r	= thrust::raw_pointer_cast( &(Nold[0]) ),
			*dx_r	= thrust::raw_pointer_cast( &(domInfo->dxD[0]) ),
			*dy_r	= thrust::raw_pointer_cast( &(domInfo->dyD[0]) ),
			*ym_r	= thrust::raw_pointer_cast( &(bc[YMINUS][0]) ),
			*yp_r	= thrust::raw_pointer_cast( &(bc[YPLUS][0]) ),
			*xm_r	= thrust::raw_pointer_cast( &(bc[XMINUS][0]) ),
			*xp_r	= thrust::raw_pointer_cast( &(bc[XPLUS][0]) );

	int		*tags_r	= thrust::raw_pointer_cast( &(tagsD[0]) ),
			*tags2_r= thrust::raw_pointer_cast( &(tags2D[0]) );

	int nx = domInfo ->nx,
		ny = domInfo ->ny;

	const int blocksize = 256;

	double	nu = (*paramDB)["flow"]["nu"].get<double>();
	double	dt = (*paramDB)["simulation"]["dt"].get<double>();

	dim3 gridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 blockU(blocksize, 1);
	kernels::Lmidx<<<gridU, blockU>>>(Lnew_r,u_r,dx_r,dy_r,nx,ny,nu);
	kernels::Lbcx<<<gridU, blockU>>>(Lnew_r, u_r, dx_r, dy_r, ym_r, yp_r, xm_r, xp_r, nx, ny, nu);

	dim3 gridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 blockV(blocksize, 1);
	kernels::Lmidy<<<gridV, blockV>>>(Lnew_r,u_r,dx_r,dy_r,nx,ny,nu);
	kernels::Lbcy<<<gridV, blockV>>>(Lnew_r, u_r, dx_r, dy_r, ym_r, yp_r, xm_r, xp_r, nx, ny, nu);

	dim3 gridUV( int( (nx*(ny-1)+(nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 blockUV(blocksize, 1);
	kernels::calcForceFadlun<<<gridUV,blockUV>>>(force_r, L_r, Lnew_r, Nold_r, N_r, u_r, uold_r, tags_r, dt, nx, ny);
}

/**
 * \brief Calculates forces acting on an immersed body (on the device).
 *
 * Uses the control volume approach explained by Lai and Peskin (2000).
 * This is a general method that can be used with any immersed boundary method.
 * It uses only the velocity and pressure fields to calculate the forces, and
 * does not involve any body forces on the immersed boundary.
 * Currently works only for one body.
 */
void NavierStokesSolver::calculateForce()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;

	double	dt = (*paramDB)["simulation"]["dt"].get<double>(),
			nu = (*paramDB)["flow"]["nu"].get<double>();

	double	*u_r		= thrust::raw_pointer_cast(&u[0]),
			*uold_r		= thrust::raw_pointer_cast(&uold[0]),
			*pressure_r	= thrust::raw_pointer_cast(&pressure[0]),
			*dxD		= thrust::raw_pointer_cast(&(domInfo->dxD[0])),
			*dyD		= thrust::raw_pointer_cast(&(domInfo->dyD[0]));

	int		*tagsIn_r	= thrust::raw_pointer_cast( &(tagsInD[0]) );

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

	kernels::dragLeftRight <<<dimGrid, dimBlock>>> (FxX_r, u_r, pressure_r, nu, dxD, dyD, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	kernels::dragBottomTop <<<dimGrid, dimBlock>>> (FxY_r, u_r, nu, dxD, dyD, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	dim3 dimGridX( int( ( (B.numCellsX[0]+1)*B.numCellsY[0]-0.5 )/blockSize )+1, 1 );
	kernels::dragUnsteady <<<dimGridX, dimBlock>>> (FxU_r, u_r, uold_r, tagsIn_r, dxD, dyD, dt, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	B.forceX[0] = thrust::reduce(FxX.begin(), FxX.end()) + thrust::reduce(FxY.begin(), FxY.end()) + thrust::reduce(FxU.begin(), FxU.end());
	//std::cout<<"FxX: "<< thrust::reduce(FxX.begin(), FxX.end())<<"\n";
	//std::cout<<"FxY: "<< thrust::reduce(FxY.begin(), FxY.end())<<"\n";
	//std::cout<<"FxU: "<< thrust::reduce(FxU.begin(), FxU.end())<<"\n";
	///std::cout<<B.forceX[0]<<"\n";
	std::cout<<"Fx: "<< thrust::reduce(FxX.begin(), FxX.end())<<" ";
	std::cout<< thrust::reduce(FxY.begin(), FxY.end())<<" ";
	std::cout<< thrust::reduce(FxU.begin(), FxU.end())<<"\n";

	// Calculating lift
	cusp::array1d<double, cusp::device_memory>
		FyX(B.numCellsY[0]+1),
		FyY(B.numCellsX[0]),
		FyU((B.numCellsX[0]+1)*B.numCellsY[0]);

	double *FyX_r = thrust::raw_pointer_cast(&FyX[0]),
	     *FyY_r = thrust::raw_pointer_cast(&FyY[0]),
	     *FyU_r = thrust::raw_pointer_cast(&FyU[0]);

	kernels::liftLeftRight <<<dimGrid, dimBlock>>> (FyX_r, u_r, nu, dxD, dyD, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	kernels::liftBottomTop <<<dimGrid, dimBlock>>> (FyY_r, u_r, pressure_r, nu, dxD, dyD, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	dim3 dimGridY( int( ( B.numCellsX[0]*(B.numCellsY[0]+1)-0.5 )/blockSize )+1, 1 );
	kernels::liftUnsteady <<<dimGridY, dimBlock>>> (FyU_r, u_r, uold_r, dxD, dyD, dt, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	B.forceY[0] = thrust::reduce(FyX.begin(), FyX.end()) + thrust::reduce(FyY.begin(), FyY.end()) + thrust::reduce(FyU.begin(), FyU.end());
}
