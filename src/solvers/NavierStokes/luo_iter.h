/***************************************************************************//**
 * \file  luo_iter.h
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \based on code by Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class luo_iter.
 */

#pragma once

#include "NavierStokesSolver.h"
#include "luo_base.h"

class luo_iter: public luo_base
{
protected:
	//////////////////////////
	//luo_iter.h
	//////////////////////////
	void updateSolver();
	void moveBody();

	//////////////////////////
	//IntermediateVelocity
	//////////////////////////
	void setVelocityInside();

public:
	//////////////////////////
	//luo_iter.cu
	//////////////////////////
	luo_iter(parameterDB *pDB=NULL, domain *dInfo=NULL);
	virtual void initialise();
	virtual void stepTime();
	virtual void writeData();
	virtual void writeCommon();
	virtual void shutDown();

	//////////////////////////
	//IntermediateVelocity
	//////////////////////////
	virtual void preRHS1Interpolation();

	virtual void cast();
};
