/***************************************************************************//**
 * \file solveStructure.inl
 * \author Christopher Minar
 * \brief Implementation of the methods of the class \c TCFSISolver
 *        to solve the structure equation on the bodies and solve for the new position
 */

#include <solvers/NavierStokes/kernels/fluidStructureInteraction.h>

template <typename memoryType>
void TCFSISolver<memoryType>::solveStructure()
{
	NavierStokesSolver<memoryType>::logger.startTimer("solveStructure");

	parameterDB  &db = *NavierStokesSolver<memoryType>::paramDB;

	real forcey,
	     Mred,
	     Ured,
	     Cy,
	     dt,
	     alpha_,
	     NumP,
	     tol;
	
	forcey = NSWithBody<memoryType>::B.forceY[0];
	Mred = 2;
	Ured = 3;
	Cy = forcey*2;
	dt  = db["simulation"]["dt"].get<real>();
	alpha_ = 0.1;
	NumP = NSWithBody<memoryType>::B.numPoints[0], //only looks at the first immeresed body
	tol = 0.000001;
	
	real *vBk_r  = thrust::raw_pointer_cast(&NSWithBody<memoryType>::B.vBk[0]),
	     *vB_r   = thrust::raw_pointer_cast(&NSWithBody<memoryType>::B.vB[0]),
	     *y_r    = thrust::raw_pointer_cast(&NSWithBody<memoryType>::B.y[0]),
	     *yk_r   = thrust::raw_pointer_cast(&NSWithBody<memoryType>::B.yk[0]),
	     *ykp1_r = thrust::raw_pointer_cast(&NSWithBody<memoryType>::B.ykp1[0]),
	     *test_r = thrust::raw_pointer_cast(&NSWithBody<memoryType>::B.test[0]);

	bool *con_r  = thrust::raw_pointer_cast(&NSWithBody<memoryType>::B.converged[0]);

	//call fsi kernel
	dim3 dimGrid(NumP,1);
	dim3 dimBlock(1,1);
	kernels::vorticityInducedVibrationsSC<<<dimGrid,dimBlock>>>(vBk_r, vB_r, y_r, ykp1_r, forcey, Mred, Ured, Cy, dt, alpha_);
	
	//check for convergence
	kernels::checkConvergencePosition<<<dimGrid,dimBlock>>>(tol, con_r, yk_r, ykp1_r);

	NavierStokesSolver<memoryType>::logger.stopTimer("solveStructure");
}
