/***************************************************************************//**
 * \file  TCFSISolver.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class \c TCFSISolver.
 */


#pragma once

#include "NSWithBody.h"


/**
 * \class TCFSISolver
 * \brief Solves the Navier-Stokes equations 
 *        using the Immersed boundary method from Taira and Colonius (2007).
 *
 * <b>The immersed boundary method: a projection approach</b> \n
 * Taira K. and Colonius T. \n
 * Journal of Computational Physics \n
 * Volume 225 Number 2 (2007).
 *
 */

template <typename memoryType>
class TCFSISolver : public NSWithBody<memoryType>
{
private:
	cusp::coo_matrix<int, real, memoryType> 
		E,		///< Interpolation matrix from the Eulerian grid to the Lagrangian points
		ET;		///< Regularization matrix form the Lagrangian points to the Eulerian grid
	bool isSub;
	double alpha;
	
	// generate the transposed gradient matrix and the interpolation matrix
	virtual void generateQT();
	
	// update the interpolation matrix with current location of body points
	void updateQT(bool isSubStep);
	
	// generate the right hand-side of the Poisson system
	virtual void generateBC2();
	
	// update the location of the bodies and re-generate appropriate matrices
	virtual void updateSolverState();
	
	// calculate forces acting on each immersed body
	virtual void calculateForce();
	
	//initialise FSI
	void initialiseFSI(int bodypoints);
	
	//solve poisson eq for pressure
	virtual void solvePoisson();
	
	//project new velocity flux
	virtual void projectionStep();

	//solve structure equation
	void solveStructure();

	//check if substep is converged
	bool checkConvergence();
	/**
	 * \brief Warning: the definition of this method does not exist.
	 */
	void generateE();
	
public:
	// constructor -- copy the database and information about the computational grid
	TCFSISolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	
	// initiliaze the solvers
	virtual void initialise();
	
	// calculate and write forces acting on each immersed body at current time
	virtual void writeData();
	
	//time stepping with fsi fun!
	virtual void stepTime();
	
	/**
	 * \brief Returns the name of the solver.
	 */
	virtual std::string name()
	{
		return "Taira & Colonius /w FSI";
	}
};
