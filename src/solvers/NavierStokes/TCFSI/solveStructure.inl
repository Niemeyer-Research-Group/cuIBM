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
<<<<<<< HEAD

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

=======
	real dt  = db["simulation"]["dt"].get<real>();
	//real dx  = NavierStokesSolver<memoryType>::domInfo->dx[ NSWithBody<memoryType>::B.I[0] ];
	real alpha_ = .25;
	real ratio = 0.125; //rhofluid/rhosolid;

	//this loop could probably be done on the gpu
	real forcex;
	//real forcey;
	//real forcexo;
	forcex = NSWithBody<memoryType>::B.forceX[0];
	//forcexo=NSWithBody<memoryType>::B.forceXold[0];
	//forcey = NSWithBody<memoryType>::B.forceY[0];
	//real Cy = 2*forcey/(1000*db["flow"["uInitial"].get<real>());
	//real Mred = 2;
	//real Ured = 4;

	//find new velocities and positions
	for (int i = 0; i < NSWithBody<memoryType>::B.totalPoints; i++){//this should be done on the gpu
		//new velocity
		//SC
		NSWithBody<memoryType>::B.uBkp1[i] = alpha_*(NSWithBody<memoryType>::B.uB[i] + ratio*dt*(forcex)) + (1-alpha_)*NSWithBody<memoryType>::B.uBk[i];
		//NSWithBody<memoryType>::B.vBkp1[i] = alpha_*(NSWithBody<memoryType>::B.vB[i] + ratio*dt*(forcey)) + (1-alpha_)*NSWithBody<memoryType>::B.vBk[i];
		//NSWithBody<memoryType>::B.vBkp1[i] = alpha_*(NSWithBody<memoryType>::B.vB[i] + dt*(Cy/(2*Mred)-4*3.14159*3.14159*NSWithBody<memoryType>::B.y[0]/(Ured*Ured))) + (1-alpha_)*NSWithBody<memoryType>::B.vBk[i];
		//lC
		//NSWithBody<memoryType>::B.uBk[i] = NSWithBody<memoryType>::B.uB[i] + ratio*dt*forcex;
		//NSWithBody<memoryType>::B.uBkp1[i] = NSWithBody<memoryType>::B.uB[i] + ratio*dt*(forcex);
		//new position
		NSWithBody<memoryType>::B.xk[i] = NSWithBody<memoryType>::B.x[i] + (NSWithBody<memoryType>::B.uB[i]+NSWithBody<memoryType>::B.uBkp1[i])*dt*0.5;
		//NSWithBody<memoryType>::B.yk[i] = NSWithBody<memoryType>::B.y[i] + (NSWithBody<memoryType>::B.vB[i]+NSWithBody<memoryType>::B.vBkp1[i])*dt*0.5;
	}
	//std::cout<<NSWithBody<memoryType>::B.uB[0]<<"\t"<<NSWithBody<memoryType>::B.uBk[0]<<std::endl;
>>>>>>> 62dcb2cbce0388289acb72183e411edebabcbef4
	NavierStokesSolver<memoryType>::logger.stopTimer("solveStructure");
}
