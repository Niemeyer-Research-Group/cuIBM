/***************************************************************************//**
 * \file  luoIBM.h
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \based on code by Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class oscCylinder.
 */

#pragma once

#include "NavierStokesSolver.h"


class luoIBM : public NavierStokesSolver
{
protected:
	cusp::array1d<int, cusp::device_memory> //names are changed to keep consistency with the luo paper, tags the same points as modifiedFadlun
		ghostTagsUV,		///< velocity nodes just inside the boundary  (tagsIn)
		ghostTagsP,			///< pressure nodes just inside the boundary  (tagsP)
		hybridTagsUV,		///< velocity nodes just outside the boundary (tags)
		hybridTagsUV2,		///< velocity nodes 2 outside the boundary    (tags2)
		hybridTagsP;		///< pressure nodes just outside the boundary (tagsPout)


	cusp::array1d<double, cusp::device_memory>
		pressureStar,
		ustar,
		body_intercept_x,
		body_intercept_y,
		image_point_x,
		image_point_y,
		body_intercept_p_x,
		body_intercept_p_y,
		image_point_p_x,
		image_point_p_y,
		distance_from_intersection_to_node,			///< distance between IB and tagged node on the device
		distance_between_nodes_at_IB,			///< distance between tags and tags2 on the device
		distance_from_u_to_body,
		distance_from_v_to_body,
		uv;									///< velocity at the IB on the device

	//testing variables
	cusp::array1d<double, cusp::device_memory>
		x1_ip,
		x2_ip,
		y1_ip,
		y2_ip,
		x1_ip_p,
		x2_ip_p,
		y1_ip_p,
		y2_ip_p,
		ip_u,
		x1,
		x2,
		x3,
		x4,
		y1,
		y2,
		y3,
		y4,
		q1,
		q2,
		q3,
		q4,
		x1_p,
		x2_p,
		x3_p,
		x4_p,
		y1_p,
		y2_p,
		y3_p,
		y4_p,
		q1_p,
		q2_p,
		q3_p,
		q4_p,
		a0,
		a1,
		a2,
		a3;


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
	void weightUhat();
	void preRHS1Interpolation();
	void testInterpX();
	void testInterpY();

	//////////////////////////
	//intermediatePressure.inl
	//////////////////////////
	void preRHS2Interpolation();
	void weightPressure();
	void testInterpP();

	//////////////////////////
	//tagpoints.inl
	//////////////////////////
	void tagPoints();
	void testOutputX();
	void testOutputY();

public:
	//constructor -- copy the database and information about the computational grid
	luoIBM(parameterDB *pDB=NULL, domain *dInfo=NULL);

	//////////////////////////
	//luoIBM.cu
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
};
