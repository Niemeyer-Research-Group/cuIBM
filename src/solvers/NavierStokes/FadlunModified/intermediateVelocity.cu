/***************************************************************************//**
 * \file
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernels that setup the intermediate velocity solve
 */
#include <solvers/NavierStokes/fadlunModified.h>

#include <solvers/NavierStokes/NavierStokes/kernels/intermediateVelocity.h>
#include <solvers/NavierStokes/FadlunModified/kernels/intermediateVelocity.h>
#include <solvers/NavierStokes/FadlunModified/kernels/LHS1.h>
#include <solvers/NavierStokes/NavierStokes/kernels/LHS1.h>
#include <solvers/NavierStokes/NavierStokes/kernels/N.h>
#include <solvers/NavierStokes/NavierStokes/kernels/L.h>

void fadlunModified::generateRHS1()
{
	const int blocksize = 256;

	dim3 dimGridUV( int( ((nx-1)*ny+nx*(ny-1)-0.5)/blocksize ) + 1, 1);
	dim3 dimBlockUV(blocksize, 1);
	dim3 dimGridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 dimBlockU(blocksize, 1);
	dim3 dimGridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 dimBlockV(blocksize, 1);

	updateRobinBoundary();

	NavierStokesSolver::Nold = NavierStokesSolver::N;
	NavierStokesSolver::generateN();

	NavierStokesSolver::generateL();

	NavierStokesSolver::generateBC1();

	kernels::generateRHS<<<dimGridUV,dimBlockUV>>>(rhs1_r, L_r, Nold_r, N_r, u_r, bc1_r, dt, nx, ny);

	kernels::updateRHS1forIBX<<<dimGridU,dimBlockU>>>(tags_r, tagsIn_r, rhs1_r, distance_from_intersection_to_node_r, distance_between_nodes_at_IB_r, uv_r, nx, ny);
	kernels::updateRHS1forIBY<<<dimGridV,dimBlockV>>>(tags_r, tagsIn_r, rhs1_r, distance_from_intersection_to_node_r, distance_between_nodes_at_IB_r, uv_r, nx, ny);
}

void fadlunModified::updateRobinBoundary()
{
	NavierStokesSolver::logger.startTimer("update Boundary");

	double 	Uinf = 1, //need a better way to enforce these, ie read from yaml file
			Vinf = 1;

	const int blocksize = 256;

	dim3 dimGridBCX( int(ny/blocksize) + 1, 1);
	dim3 dimGridBCY( int(nx/blocksize) + 1, 1);
	dim3 dimBlockBC(blocksize, 1);

	kernels::updateBoundaryX<<<dimGridBCX,dimBlockBC>>>(u_r, xp_r, dx_r, dt, Uinf, nx, ny);
	kernels::updateBoundaryY<<<dimGridBCY,dimBlockBC>>>(u_r, xp_r, dx_r, dt, Vinf, nx, ny);

	NavierStokesSolver::logger.stopTimer("update Boundary");
}

void fadlunModified::generateLHS1()
{
	NavierStokesSolver::logger.startTimer("LHS1");

	const int blocksize = 256;

	dim3 gridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 blockU(blocksize, 1);
	dim3 gridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 blockV(blocksize, 1);

	kernels::LHS_mid_X<<<gridU,blockU>>>(LHS1_row_r, LHS1_col_r, LHS1_val_r, tags_r, tags2_r, tagsIn_r, distance_from_intersection_to_node_r, distance_between_nodes_at_IB_r, dx_r, dy_r, dt, nu, nx, ny);
	kernels::LHS_mid_Y<<<gridV,blockV>>>(LHS1_row_r, LHS1_col_r, LHS1_val_r, tags_r, tags2_r, tagsIn_r, distance_from_intersection_to_node_r, distance_between_nodes_at_IB_r, dx_r, dy_r, dt, nu, nx, ny);

	kernels::LHS_BC_X<<<gridU,blockU>>>(LHS1_row_r, LHS1_col_r, LHS1_val_r, dx_r, dy_r, dt, nu, nx, ny);
	kernels::LHS_BC_Y<<<gridV,blockV>>>(LHS1_row_r, LHS1_col_r, LHS1_val_r, dx_r, dy_r, dt, nu, nx, ny);

	NavierStokesSolver::logger.stopTimer("LHS1");
}
