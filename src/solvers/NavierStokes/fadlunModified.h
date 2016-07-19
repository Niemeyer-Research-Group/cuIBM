/***************************************************************************//**
 * \file  fadlunModified.h
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \based on code by Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class oscCylinder.
 */

#pragma once

#include "NavierStokesSolver.h"


class fadlunModified : public NavierStokesSolver
{
protected:
	cusp::array1d<int, cusp::device_memory>
		tags,		///< indices of velocity nodes near the IB on the device
		tagsOld,
		tagsPOld,
		tags2,		///< indices of 2nd closest velocity node on the device
		tagsIn,	///< indices of velocity nodes inside the IB
		tagsP, 	///< indices of pressure nodes inside the IB
		tagsPOut;	///< indices of pressure nodes outside the IB

	cusp::array1d<double, cusp::device_memory>
		distance_from_intersection_to_node,			///< distance between IB and tagged node on the device
		distance_between_nodes_at_IB,			///< distance between tags and tags2 on the device
		distance_from_u_to_body,
		distance_from_v_to_body,
		uv,		///< velocity at the IB on the device
		test;


	int	*tags_r,
		*tagsOld_r,
		*tagsPOld_r,
		*tags2_r,
		*tagsIn_r,
		*tagsP_r,
		*tagsPOut_r;


	double	*distance_from_intersection_to_node_r,
			*distance_between_nodes_at_IB_r,
			*distance_from_u_to_body_r,
			*distance_from_v_to_body_r,
			*uv_r,
			*test_r;

	bodies 	B;		///< bodies in the flow

	std::ofstream forceFile;

	//////////////////////////
	//calculateForce.inl
	//////////////////////////
	void calculateForce();

	//////////////////////////
	//intermediateVelocity.inl
	//////////////////////////
	void updateRobinBoundary();

	//////////////////////////
	//tagpoints.inl
	//////////////////////////
	void tagPoints();

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
	virtual void stepTime();
	virtual void shutDown();

	//////////////////////////
	//intermediatePressure.inl
	//////////////////////////
	virtual void generateRHS2();
	virtual void generateLHS2();

	//////////////////////////
	//intermediateVelocity.inl
	//////////////////////////
	virtual void generateRHS1();
	virtual void generateLHS1();

	//////////////////////////
	//projectVelocity.inl
	//////////////////////////
	virtual void velocityProjection();

	//////////////////////////
	//cusp.cu
	//////////////////////////
	virtual void cast();
};
