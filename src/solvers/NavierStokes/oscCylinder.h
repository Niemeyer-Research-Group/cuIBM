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
	cusp::array1d<double, cusp::device_memory> cfl;
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
public:
	//////////////////////////
	//oscCylinder.cu
	//////////////////////////
	oscCylinder(parameterDB *pDB=NULL, domain *dInfo=NULL);
	virtual void initialise();
	virtual void stepTime();
	virtual void writeData();
	virtual void writeCommon();
	virtual void shutDown();

	//////////////////////////
	//IntermediateVelocity.inl
	//////////////////////////
	virtual void preRHS1Interpolation();
};
