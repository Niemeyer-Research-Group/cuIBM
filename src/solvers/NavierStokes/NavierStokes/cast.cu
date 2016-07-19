/***************************************************************************//**
 * \file intermediateVelocity.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief thrust remote cast of all 1darrays
 */

#include <solvers/NavierStokes/NavierStokesSolver.h>

void NavierStokesSolver::cast()
{
	//note: try to resize things before they are cast, it doesn't really like it the other way around
	nu = (*paramDB)["flow"]["nu"].get<double>();
	dt = (*paramDB)["simulation"]["dt"].get<double>();

	nx = domInfo ->nx,
	ny = domInfo ->ny;

	numU = (nx-1)*ny;
	numV = (ny-1)*nx;
	numUV = numU+numV;
	numP = nx*ny;

	//resize arrays
	LHS1.resize(numUV, numUV, (nx-1)*ny*5 - 2*ny-2*(nx-1)       +        (ny-1)*nx*5 - 2*(ny-1) - 2*nx);
	LHS2.resize(nx*ny, nx*ny, 5*nx*ny - 2*ny-2*nx);

	//numUV
	u.resize(numUV);
	uhat.resize(numUV);
	uold.resize(numUV);
	N.resize(numUV);
	Nold.resize(numUV);
	L.resize(numUV);
	Lnew.resize(numUV);
	rhs1.resize(numUV);
	bc1.resize(numUV);

	//numP
	pressure.resize(numP);
	rhs2.resize(numP);
	pressure_old.resize(numP);

	//Boundary conditions
	bc[XMINUS].resize(2*ny-1);
	bc[XPLUS].resize(2*ny-1);
	bc[YMINUS].resize(2*nx-1);
	bc[YPLUS].resize(2*nx-1);

	//cast arrays
	u_r				= thrust::raw_pointer_cast( &(u[0]) ),
	uhat_r			= thrust::raw_pointer_cast( &(uhat[0]) ),
	uold_r			= thrust::raw_pointer_cast( &(uold[0]) ),
	pressure_r		= thrust::raw_pointer_cast( &(pressure[0]) ),
	pressure_old_r	= thrust::raw_pointer_cast( &(pressure_old[0]) ),
	Nold_r			= thrust::raw_pointer_cast( &(Nold[0]) ),
	N_r				= thrust::raw_pointer_cast( &(N[0]) ),
	L_r				= thrust::raw_pointer_cast( &(L[0]) ),
	Lnew_r			= thrust::raw_pointer_cast( &(Lnew[0]) ),
	bc1_r			= thrust::raw_pointer_cast( &(bc1[0]) ),
	rhs1_r			= thrust::raw_pointer_cast( &(rhs1[0]) ),
	rhs2_r			= thrust::raw_pointer_cast( &(rhs2[0]) ),
	xp_r			= thrust::raw_pointer_cast( &(bc[XPLUS][0]) ),
	xm_r			= thrust::raw_pointer_cast( &(bc[XMINUS][0]) ),
	yp_r			= thrust::raw_pointer_cast( &(bc[YPLUS][0]) ),
	ym_r			= thrust::raw_pointer_cast( &(bc[YMINUS][0]) ),
	x_r				= thrust::raw_pointer_cast( &(domInfo->x[0]) ),
	y_r				= thrust::raw_pointer_cast( &(domInfo->y[0]) ),
	dx_r			= thrust::raw_pointer_cast( &(domInfo->dx[0]) ),
	dy_r			= thrust::raw_pointer_cast( &(domInfo->dy[0]) ),
	xu_r			= thrust::raw_pointer_cast( &(domInfo->xu[0]) ),
	xv_r			= thrust::raw_pointer_cast( &(domInfo->xv[0]) ),
	yu_r			= thrust::raw_pointer_cast( &(domInfo->yu[0]) ),
	yv_r			= thrust::raw_pointer_cast( &(domInfo->yv[0]) );
	LHS1_row_r		= thrust::raw_pointer_cast( &(LHS1.row_indices[0]) );
	LHS1_col_r		= thrust::raw_pointer_cast( &(LHS1.column_indices[0]) );
	LHS1_val_r		= thrust::raw_pointer_cast( &(LHS1.values[0]) );
	LHS2_row_r		= thrust::raw_pointer_cast( &(LHS2.row_indices[0]) );
	LHS2_col_r		= thrust::raw_pointer_cast( &(LHS2.column_indices[0]) );
	LHS2_val_r		= thrust::raw_pointer_cast( &(LHS2.values[0]) );
}
