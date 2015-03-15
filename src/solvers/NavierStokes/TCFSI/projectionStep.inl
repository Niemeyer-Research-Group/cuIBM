/***************************************************************************//**
 * \file projectionStep.inl
 * \author Christopher Minar
 * \brief Implementation of the methods of the class \c TCFSISolver
 *        to calculate q@n+1.
 */


/**
 * \brief Projects the flux onto the divergence-free field.
 * \      calculate qn+1
 */
template <typename memoryType>
void TCFSISolver<memoryType>::projectionStep()
{
	NavierStokesSolver<memoryType>::logger.startTimer("projectionStep");

	cusp::multiply(NavierStokesSolver<memoryType>::Q, NavierStokesSolver<memoryType>::lambdak, NavierStokesSolver<memoryType>::temp1);
	//BN is the inverse of the mass matrix, it doens't need to be remade
	cusp::multiply(NavierStokesSolver<memoryType>::BN, NavierStokesSolver<memoryType>::temp1, NavierStokesSolver<memoryType>::qk);
	cusp::blas::axpby(NavierStokesSolver<memoryType>::qStar, NavierStokesSolver<memoryType>::qk, NavierStokesSolver<memoryType>::qk, 1.0, -1.0);

	NavierStokesSolver<memoryType>::logger.stopTimer("projectionStep");
}
