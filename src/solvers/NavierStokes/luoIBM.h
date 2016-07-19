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
		body_intercept_p,
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
		image_point_u,
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
		a3,
		dudt,
		ududx,
		vdudy,
		dvdt,
		udvdx,
		vdvdy;

	//interp variables
	cusp::array1d<double, cusp::device_memory>
		detA,
		alpha,
		b11,
		b12,
		b13,
		b14,
		b21,
		b22,
		b23,
		b24,
		b31,
		b32,
		b33,
		b34,
		b41,
		b42,
		b43,
		b44,
		stencilCoef,
		interpCoef;

	cusp::array1d<int, cusp::device_memory>
		countD;

	cusp::array1d<int, cusp::host_memory>
		countH;

	cusp::array1d<bool, cusp::device_memory>
		q1flag,
		q2flag,
		q3flag,
		q4flag;

	cusp::array1d<int, cusp::device_memory>
		index1,
		index2,
		index3,
		index4;


	bodies 	B;		///< bodies in the flow

	std::ofstream forceFile;

	//////////////////////////
	//calculateForce.inl
	//////////////////////////
	void calculateForce();
	void luoForce();

	//////////////////////////
	//intermediateVelocity.inl
	//////////////////////////
	void updateRobinBoundary();
	void weightUhat();
	void zeroVelocity();

	//////////////////////////
	//intermediatePressure.inl
	//////////////////////////
	void preRHS2();
	void interpPGN();
	void sizeLHS2();

	//////////////////////////
	//tagpoints.inl
	//////////////////////////
	void tagPoints();

	//////////////////////////
	//cast.inl
	//////////////////////////
	void castluo();

	//////////////////////////
	//testing.inl
	//////////////////////////
	void divergence();
	void testInterpX(); //x
	void testInterpY(); //y
	void testInterpP(); //for pressure
	void testOutputX(); //for tagpoipnts
	void testOutputY(); //for tagpoints
	void testForce_p();
	void testForce_dudn();


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
	void outputPressure();
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
	virtual void preRHS1Interpolation();
	virtual void generateLHS1();

	//////////////////////////
	//projectVelocity.inl
	//////////////////////////
	virtual void velocityProjection();
};
