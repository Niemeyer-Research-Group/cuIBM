/***************************************************************************//**
 * \file solvePoisson.inl
 * \author Christopher Minar
 * \brief Implementation of the methods of the class \c TCFSISolver
 *        to solve the Poisson system.
 */

/**
 * \brief Solves the Poisson system for the pressure and body forces
 */


template <typename memoryType>
void TCFSISolver<memoryType>::solvePoisson()
{
	NavierStokesSolver<memoryType>::logger.startTimer("solvePoisson");
	
	parameterDB *DB = NavierStokesSolver<memoryType>::paramDB;

	int  maxIters = (*DB)["PoissonSolve"]["maxIterations"].get<int>();
	real relTol   = (*DB)["PoissonSolve"]["tolerance"].get<real>();
	
	cusp::default_monitor<real> sys2Mon(NavierStokesSolver<memoryType>::rhs2, maxIters, relTol);
	//use lamdak for substep
	cusp::krylov::bicgstab(NavierStokesSolver<memoryType>::C, NavierStokesSolver<memoryType>::lambdak, NavierStokesSolver<memoryType>::rhs2, sys2Mon, *NavierStokesSolver<memoryType>::PC2);
	
		
	NavierStokesSolver<memoryType>::iterationCount2 = sys2Mon.iteration_count();
	if (!sys2Mon.converged())
	{
		std::cout << "ERROR: Solve for Lambda failed at time step " << NavierStokesSolver<memoryType>::timeStep << std::endl;
		std::cout << "Iterations   : " << NavierStokesSolver<memoryType>::iterationCount2 << std::endl;          
		std::cout << "Residual norm: " << sys2Mon.residual_norm() << std::endl;
		std::cout << "Tolerance    : " << sys2Mon.tolerance() << std::endl;
		std::exit(-1);
	}
	
	NavierStokesSolver<memoryType>::logger.stopTimer("solvePoisson");
}

