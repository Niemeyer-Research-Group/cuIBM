/***************************************************************************//**
 * \file  oscCylinder.h
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \based on code by Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class oscCylinder.
 */

#pragma once

#include "NavierStokesSolver.h"


class oscCylinder : public NavierStokesSolver
{

public:
	//constructor -- copy the database and information about the computational grid
	oscCylinder(parameterDB *pDB=NULL, domain *dInfo=NULL);

	//step forward in time
	virtual void stepTime();

	//recalculate LHS matrices
	void updateSolver();

	//write stuff
	virtual void writeData();

	//perform motion calculation
	void moveBody();

	void initialise();
};
