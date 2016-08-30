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

public:
	//////////////////////////
	//luo_iter.cu
	//////////////////////////
	luo_iter(parameterDB *pDB=NULL, domain *dInfo=NULL);
	virtual void initialise();
	virtual void writeData();
	virtual void writeCommon();
	virtual void _intermediate_velocity();
	virtual void _pressure();

	//////////////////////////
	//IntermediateVelocity.cu
	//////////////////////////
	void intermediate_velocity_setup();
		void intermediate_velocity_alpha();
		void intermediate_velocity_interpolation_setup();
		void intermediate_velocity_size_lhs();
		void intermediate_velocity_calculate_lhs();
		void intermediate_velocity_update_rhs();

	//////////////////////////
	//Poisson.cu
	//////////////////////////
	void poisson_setup();
		void poisson_alpha();
		void poisson_interpolation_setup();
		void poisson_size_lhs();
		void poisson_calculate_lhs();
		void poisson_update_rhs();

	virtual void cast();

	//test
	void generateRHS2();
	void weightPressure();
};
