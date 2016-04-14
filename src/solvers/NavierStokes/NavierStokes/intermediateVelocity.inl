/***************************************************************************//**
 * \file intermediateVelocity.inl
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernels that setup the intermediate velocity solve
 */

#include <solvers/NavierStokes/kernels/intermediateVelocity.h>
#include <solvers/NavierStokes/kernels/LHS1.h>
#include <solvers/NavierStokes/kernels/N.h>
#include <solvers/NavierStokes/kernels/L.h>

void NavierStokesSolver::generateRHS1()
{
	double	*u_r    = thrust::raw_pointer_cast( &(u[0]) ),
			*uold_r = thrust::raw_pointer_cast( &(uold[0]) ),
			*N_r    = thrust::raw_pointer_cast( &(N[0]) ),
			*Nold_r = thrust::raw_pointer_cast( &(Nold[0]) ),
			*L_r	= thrust::raw_pointer_cast( &(L[0]) ),
			*rhs_r	= thrust::raw_pointer_cast( &(rhs1[0]) ),
			*bc1_r	= thrust::raw_pointer_cast( &(bc1[0]) ),
			*dxD_r	= thrust::raw_pointer_cast( &(domInfo->dxD[0]) ),
			*dyD_r	= thrust::raw_pointer_cast( &(domInfo->dyD[0]) ),
			*x_r	= thrust::raw_pointer_cast( &(domInfo->x[0]) ),
			*y_r	= thrust::raw_pointer_cast( &(domInfo->y[0]) ),
			*dx_r	= thrust::raw_pointer_cast( &(domInfo->dx[0]) ),
			*dy_r	= thrust::raw_pointer_cast( &(domInfo->dy[0]) ),
			*xp_r	= thrust::raw_pointer_cast( &(bc[XPLUS][0]) ),
			*a_r	= thrust::raw_pointer_cast( &(aD[0]) ),
			*b_r	= thrust::raw_pointer_cast( &(bD[0]) ),
			*uv_r	= thrust::raw_pointer_cast( &(uvD[0]) );
	
	int	*tags_r	= thrust::raw_pointer_cast( &(tagsD[0]) ),
		*tagsIn_r=thrust::raw_pointer_cast( &(tagsInD[0]) );

	double	nu = (*paramDB)["flow"]["nu"].get<double>();
	double	dt = (*paramDB)["simulation"]["dt"].get<double>();

	int nx = domInfo ->nx,
		ny = domInfo ->ny;
	
	//set correct grid and block size
	const int blocksize = 256;
	dim3 dimGridUV( int( ((nx-1)*ny+nx*(ny-1)-0.5)/blocksize ) + 1, 1);
	dim3 dimBlockUV(blocksize, 1);
	dim3 dimGridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 dimBlockU(blocksize, 1);
	dim3 dimGridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 dimBlockV(blocksize, 1);

	updateRobinBoundary();

	Nold = N;
	generateN();

	generateL();

	generateBC1();

	kernels::generateRHS<<<dimGridUV,dimBlockUV>>>(rhs_r, L_r, Nold_r, N_r, u_r, bc1_r, dt, nx, ny);
	if (B.numBodies > 0)
	{
		kernels::updateRHS1forIBX<<<dimGridU,dimBlockU>>>(tags_r, tagsIn_r, rhs_r, a_r, b_r, uv_r, nx, ny);
		kernels::updateRHS1forIBY<<<dimGridV,dimBlockV>>>(tags_r, tagsIn_r, rhs_r, a_r, b_r, uv_r, nx, ny);
	}
}

void NavierStokesSolver::updateRobinBoundary()
{
	logger.startTimer("update Boundary");
	double	*u_r    = thrust::raw_pointer_cast( &(u[0]) ),
			*xp_r	= thrust::raw_pointer_cast( &(bc[XPLUS][0]) ),
			*dxD_r	= thrust::raw_pointer_cast( &(domInfo->dxD[0]) );
			//*dyD_r	= thrust::raw_pointer_cast( &(domInfo->dyD[0]) );
	double 	Uinf = 1, //need a better way to enforce these, ie read from yaml file
			Vinf = 1;
	int 	nx = domInfo ->nx,
			ny = domInfo ->ny;
	double	dt = (*paramDB)["simulation"]["dt"].get<double>();

	const int blocksize = 256;
	dim3 dimGridBCX( int(ny/blocksize) + 1, 1);
	dim3 dimGridBCY( int(nx/blocksize) + 1, 1);
	dim3 dimBlockBC(blocksize, 1);

	if(B.numBodies>0)
	{
		kernels::updateBoundaryX<<<dimGridBCX,dimBlockBC>>>(u_r, xp_r, dxD_r, dt, Uinf, nx, ny);
		kernels::updateBoundaryY<<<dimGridBCY,dimBlockBC>>>(u_r, xp_r, dxD_r, dt, Vinf, nx, ny);
	}
	logger.stopTimer("update Boundary");
}

/**
 * \brief Fills the Convection matrix N
 * \ param nx number of pressure nodes in the x direction
 * \ param ny number of pressure nodes in the y direction
 * \ param dx array containing the X-direction cell widths
 * \ param dy array containing the Y-direction cell heights
 * \ param N  convection matrix (stored as an array)
 * \ param u  velocity matrix (stored as an array)
 */
void NavierStokesSolver::generateN()
{
	logger.startTimer("Advection Terms");
	double	*u_r    = thrust::raw_pointer_cast( &(u[0]) ),
			*N_r    = thrust::raw_pointer_cast( &(N[0]) ),
			*dx_r	= thrust::raw_pointer_cast( &(domInfo->dxD[0]) ),
			*dy_r	= thrust::raw_pointer_cast( &(domInfo->dyD[0]) ),
			*ym_r	= thrust::raw_pointer_cast( &(bc[YMINUS][0]) ),
			*yp_r	= thrust::raw_pointer_cast( &(bc[YPLUS][0]) ),
			*xm_r	= thrust::raw_pointer_cast( &(bc[XMINUS][0]) ),
			*xp_r	= thrust::raw_pointer_cast( &(bc[XPLUS][0]) );

	int nx = domInfo ->nx,
		ny = domInfo ->ny;

	const int blocksize = 256;

	dim3 gridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 blockU(blocksize, 1);
	kernels::Nmidx<<<gridU, blockU>>>(N_r,u_r,dx_r,dy_r,nx,ny);
	kernels::Nbcx<<<gridU, blockU>>>(N_r, u_r, dx_r, dy_r, ym_r, yp_r, xm_r, xp_r, nx, ny);
	
	dim3 gridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 blockV(blocksize, 1);
	kernels::Nmidy<<<gridV, blockV>>>(N_r,u_r,dx_r,dy_r,nx,ny);
	kernels::Nbcy<<<gridV, blockV>>>(N_r, u_r, dx_r, dy_r, ym_r, yp_r, xm_r, xp_r, nx, ny);
	logger.stopTimer("Advection Terms");
}

/**
 * \brief Fills the Laplacian matrix L
 * \ param nx number of pressure nodes in the x direction
 * \ param ny number of pressure nodes in the y direction
 * \ param nu viscosity: effectivly the inverse Reynolds number
 * \ param dx array containing the X-direction cell widths
 * \ param dy array containing the Y-direction cell heights
 * \ param L  laplacian matrix (stored as an array)
 * \ param u  velocity matrix (stored as an array)
 */
void NavierStokesSolver::generateL()
{
	logger.startTimer("L");
	double	*u_r    = thrust::raw_pointer_cast( &(u[0]) ),
			*L_r    = thrust::raw_pointer_cast( &(L[0]) ),
			*dx_r	= thrust::raw_pointer_cast( &(domInfo->dxD[0]) ),
			*dy_r	= thrust::raw_pointer_cast( &(domInfo->dyD[0]) ),
			*ym_r	= thrust::raw_pointer_cast( &(bc[YMINUS][0]) ),
			*yp_r	= thrust::raw_pointer_cast( &(bc[YPLUS][0]) ),
			*xm_r	= thrust::raw_pointer_cast( &(bc[XMINUS][0]) ),
			*xp_r	= thrust::raw_pointer_cast( &(bc[XPLUS][0]) );

	int nx = domInfo ->nx,
		ny = domInfo ->ny;

	const int blocksize = 256;

	double	nu = (*paramDB)["flow"]["nu"].get<double>();

	dim3 gridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 blockU(blocksize, 1);
	kernels::Lmidx<<<gridU, blockU>>>(L_r,u_r,dx_r,dy_r,nx,ny,nu);
	kernels::Lbcx<<<gridU, blockU>>>(L_r, u_r, dx_r, dy_r, ym_r, yp_r, xm_r, xp_r, nx, ny, nu);

	dim3 gridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 blockV(blocksize, 1);
	kernels::Lmidy<<<gridV, blockV>>>(L_r,u_r,dx_r,dy_r,nx,ny,nu);
	kernels::Lbcy<<<gridV, blockV>>>(L_r, u_r, dx_r, dy_r, ym_r, yp_r, xm_r, xp_r, nx, ny, nu);
	logger.stopTimer("L");
}

void NavierStokesSolver::generateLHS1()
{
	logger.startTimer("LHS1");

	double 	*val_r	= thrust::raw_pointer_cast( &(LHS1.values[0])),
			*dx_r	= thrust::raw_pointer_cast( &(domInfo->dxD[0])),
			*dy_r	= thrust::raw_pointer_cast( &(domInfo->dyD[0])),
			*a_r	= thrust::raw_pointer_cast( &(aD[0]) ),
			*b_r	= thrust::raw_pointer_cast( &(bD[0]) );

	int 	*row_r	= thrust::raw_pointer_cast( &(LHS1.row_indices[0])),
			*col_r	= thrust::raw_pointer_cast( &(LHS1.column_indices[0])),
			*tags_r	= thrust::raw_pointer_cast( &(tagsD[0]) ),
			*tags2_r= thrust::raw_pointer_cast( &(tags2D[0]) ),
			*tagsIn_r= thrust::raw_pointer_cast( &(tagsInD[0]) );

	double	nu = (*paramDB)["flow"]["nu"].get<double>();
	double	dt = (*paramDB)["simulation"]["dt"].get<double>();
	int 	nx = domInfo->nx,
			ny = domInfo->ny;

	const int blocksize = 256;
	dim3 gridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 blockU(blocksize, 1);
	dim3 gridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 blockV(blocksize, 1);
	if (B.numBodies > 0)
	{
		kernels::LHS_mid_X<<<gridU,blockU>>>(row_r, col_r, val_r, tags_r, tags2_r, tagsIn_r, a_r, b_r, dx_r, dy_r, dt, nu, nx, ny);
		kernels::LHS_mid_Y<<<gridV,blockV>>>(row_r, col_r, val_r, tags_r, tags2_r, tagsIn_r, a_r, b_r, dx_r, dy_r, dt, nu, nx, ny);
	}
	else
	{
		kernels::LHS_mid_X_nobody<<<gridU,blockU>>>(row_r, col_r, val_r, dx_r, dy_r, dt, nu, nx, ny);
		kernels::LHS_mid_Y_nobody<<<gridV,blockV>>>(row_r, col_r, val_r, dx_r, dy_r, dt, nu, nx, ny);
	}
	kernels::LHS_BC_X<<<gridU,blockU>>>(row_r, col_r, val_r, dx_r, dy_r, dt, nu, nx, ny);
	kernels::LHS_BC_Y<<<gridV,blockV>>>(row_r, col_r, val_r, dx_r, dy_r, dt, nu, nx, ny);

	logger.stopTimer("LHS1");
}

void NavierStokesSolver::generateBC1()
{
	logger.startTimer("BC1");
	double	*u_r    = thrust::raw_pointer_cast( &(u[0]) ),
			*bc1_r	= thrust::raw_pointer_cast( &(bc1[0]) ),
			*dx_r	= thrust::raw_pointer_cast( &(domInfo->dxD[0]) ),
			*dy_r	= thrust::raw_pointer_cast( &(domInfo->dyD[0]) ),
			*ym_r	= thrust::raw_pointer_cast( &(bc[YMINUS][0]) ),
			*yp_r	= thrust::raw_pointer_cast( &(bc[YPLUS][0]) ),
			*xm_r	= thrust::raw_pointer_cast( &(bc[XMINUS][0]) ),
			*xp_r	= thrust::raw_pointer_cast( &(bc[XPLUS][0]) );

	int nx = domInfo ->nx,
		ny = domInfo ->ny;

	double	nu = (*paramDB)["flow"]["nu"].get<double>();
	double	dt = (*paramDB)["simulation"]["dt"].get<double>();

	const int blocksize = 256;

	dim3 gridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 blockU(blocksize, 1);
	kernels::bc1X<<<gridU, blockU>>>(u_r, bc1_r, ym_r, yp_r, xm_r, xp_r, dx_r, dy_r, nu, dt, nx, ny);

	dim3 gridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 blockV(blocksize, 1);
	kernels::bc1Y<<<gridV, blockV>>>(u_r, bc1_r, ym_r, yp_r, xm_r, xp_r, dx_r, dy_r, nu, dt, nx, ny);
	logger.stopTimer("BC1");
}
