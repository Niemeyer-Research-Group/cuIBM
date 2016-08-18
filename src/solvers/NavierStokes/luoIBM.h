/***************************************************************************//**
 * \file  luoIBM.h
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \based on code by Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class oscCylinder.
 */

#pragma once

#include "luo_base.h"


class luoIBM : public luo_base
{
protected:
	//////////////////////////
	//intermediateVelocity
	//////////////////////////
	void weightUhat();
	void zeroVelocity();

	//////////////////////////
	//intermediatePressure
	//////////////////////////
	void preRHS2Interpolation();
	void weightPressure();

	//////////////////////////
	//testing
	//////////////////////////

public:
	//constructor -- copy the database and information about the computational grid
	luoIBM(parameterDB *pDB=NULL, domain *dInfo=NULL);

	//////////////////////////
	//luoIBM
	//////////////////////////
	virtual void initialise();
	virtual void initialiseLHS();
	virtual void writeData();
	virtual void writeCommon();
	virtual void _intermediate_velocity();
	virtual void _pressure();
	virtual void shutDown();

	//////////////////////////
	//intermediatePressure
	//////////////////////////
	virtual void generateRHS2();
	virtual void generateLHS2();

	//////////////////////////
	//intermediateVelocity
	//////////////////////////
	virtual void generateRHS1();
	virtual void preRHS1Interpolation();
	virtual void generateLHS1();

	virtual void cast();
};
