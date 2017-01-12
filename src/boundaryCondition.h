/***************************************************************************//**
 * \file boundaryCondition.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Definition of the class \c boundaryCondition.
 * \This class is outdated from the old cuibm. setting a bc to dirichlet does nothing. All bcs are dirichelt except the ones that are hard coded.
 */


#pragma once

#include <string>
#include <sstream>
#include "bcTypes.h"
#include "parameterDB.h"
#include "boundaryTypes.h"


/**
 * \class boundaryCondition
 * \brief Stores the type of boundary condition and its value.
 */
class boundaryCondition
{
public:
	bcType type; ///< type of boundary condition
	double  value; ///< numerical value associated with the boundary condition
	
	/**
	 * \brief Constructor. Sets Dirichlet boundary condition with value.
	 */
	boundaryCondition() : type(DIRICHLET), value(0) {};

	/**
	 * \brief Constructor overloading. Sets boundary condition to a given type
	 *        with a given value.
	 *
	 * \param _type the type of boundary condition
	 * \param _value the value at the boundary
	 */
	boundaryCondition(bcType _type, double _value) : type(_type), value(_value) {}; 
};
