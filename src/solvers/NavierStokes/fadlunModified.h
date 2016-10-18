/***************************************************************************//**
 * \file  fadlunModified.h
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \based on code by Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class oscCylinder.
 */

#pragma once

#include "NavierStokesSolver.h"
#include "luo_base.h"

class fadlunModified : public luo_base
{
public:
	//constructor -- copy the database and information about the computational grid
	fadlunModified(parameterDB *pDB=NULL, domain *dInfo=NULL);

	//////////////////////////
	//fadlunModified.cu
	//////////////////////////
	virtual void initialise();
	virtual void initialiseLHS();
	virtual void writeData();
	virtual void writeCommon();
	virtual void shutDown();
	virtual void _intermediate_velocity();
	virtual void _pressure();

	//////////////////////////
	//intermediatePressure
	//////////////////////////
	virtual void generateRHS2();
	virtual void generateLHS2();

	//////////////////////////
	//intermediateVelocity
	//////////////////////////
	virtual void generateRHS1();
	virtual void generateLHS1();

	//////////////////////////
	//cast
	//////////////////////////
	virtual void cast();
	void calculateForce();
};
