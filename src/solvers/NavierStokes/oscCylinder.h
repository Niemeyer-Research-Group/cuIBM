/***************************************************************************//**
 * \file  oscCylinder.h
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \based on code by Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class oscCylinder.
 */

#pragma once

#include "NavierStokesSolver.h"
#include "luoIBM.h"

class oscCylinder : public luoIBM
{
protected:
	std::ofstream midPositionFile;
<<<<<<< HEAD
	cusp::array1d<double, cusp::device_memory> cfl,
												distance;
	double	cfl_max,
			cfl_I,
			cfl_J,
			cfl_ts;

	//////////////////////////
	//oscCylinder.h
	//////////////////////////
	void updateSolver();
	void moveBody();

	//////////////////////////
	//IntermediateVelocity.inl
	//////////////////////////
	void setVelocityInside();

	//////////////////////////
	//CFL.inl
	//////////////////////////
	void CFL();
	void calcDistance();
=======
>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder
public:
	//constructor -- copy the database and information about the computational grid
	oscCylinder(parameterDB *pDB=NULL, domain *dInfo=NULL);

	//step forward in time
	virtual void stepTime();

	//recalculate LHS matrices
	void updateSolver();

	//write stuff
	virtual void writeData();

	virtual void writeCommon();

	//perform motion calculation
	void moveBody();

	virtual void initialise();

	virtual void shutDown();
};
