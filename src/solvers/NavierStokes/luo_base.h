/***************************************************************************//**
 * \file  luo_base.h
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \based on code by Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class oscCylinder.
 */

#pragma once

#include "NavierStokesSolver.h"


class luo_base : public NavierStokesSolver
{
protected:
	cusp::array1d<int, cusp::device_memory> //names are changed to keep consistency with the luo paper, tags the same points as modifiedFadlun
		ghostTagsUV,		///< velocity nodes just inside the boundary  (tagsIn)
		ghostTagsP,			///< pressure nodes just inside the boundary  (tagsP)
		hybridTagsUV,		///< velocity nodes just outside the boundary (tags)
		hybridTagsUV2,		///< velocity nodes 2 outside the boundary    (tags2)
		hybridTagsP,		///< pressure nodes just outside the boundary (tagsPout)
		index1,
		index2,
		index3,
		index4,
		count;


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
		alpha,
		q1coef,
		q2coef,
		q3coef,
		q4coef,
		ns_rhs,
		interp_rhs,
		dpdn;

	int *ghostTagsUV_r,
		*ghostTagsP_r,
		*hybridTagsUV_r,
		*hybridTagsP_r,
		*hybridTagsUV2_r,
		*index1_r,
		*index2_r,
		*index3_r,
		*index4_r,
		*count_r;

	double	*pressureStar_r,
			*ustar_r,
			*body_intercept_x_r,
			*body_intercept_y_r,
			*image_point_x_r,
			*image_point_y_r,
			*body_intercept_p_x_r,
			*body_intercept_p_y_r,
			*body_intercept_p_r,
			*image_point_p_x_r,
			*image_point_p_y_r,
			*distance_from_intersection_to_node_r,
			*distance_between_nodes_at_IB_r,
			*distance_from_u_to_body_r,
			*distance_from_v_to_body_r,
			*uv_r,
			*alpha_r;

	double	*x1_ip_r,
			*x2_ip_r,
			*y1_ip_r,
			*y2_ip_r,
			*x1_ip_p_r,
			*x2_ip_p_r,
			*y1_ip_p_r,
			*y2_ip_p_r,
			*image_point_u_r,
			*x1_r,
			*x2_r,
			*x3_r,
			*x4_r,
			*y1_r,
			*y2_r,
			*y3_r,
			*y4_r,
			*q1_r,
			*q2_r,
			*q3_r,
			*q4_r,
			*x1_p_r,
			*x2_p_r,
			*x3_p_r,
			*x4_p_r,
			*y1_p_r,
			*y2_p_r,
			*y3_p_r,
			*y4_p_r,
			*q1_p_r,
			*q2_p_r,
			*q3_p_r,
			*q4_p_r,
			*a0_r,
			*a1_r,
			*a2_r,
			*a3_r,
			*q1coef_r,
			*q2coef_r,
			*q3coef_r,
			*q4coef_r,
			*ns_rhs_r,
			*interp_rhs_r,
			*dpdn_r;

	bodies 	B;		///< bodies in the flow

	std::ofstream forceFile;
	std::ofstream midPositionFile;

	//////////////////////////
	//calculateForce
	//////////////////////////
	void calculateForce();
	void luoForce();

	//////////////////////////
	//intermediateVelocity
	//////////////////////////
	void updateRobinBoundary();

	//////////////////////////
	//tagpoints
	//////////////////////////
	void tagPoints();

	//////////////////////////
	//move
	//////////////////////////
	void set_movement();
	void viv_movement_LC();
	void viv_movement_SC();

	//////////////////////////
	//testing
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
	luo_base(parameterDB *pDB=NULL, domain *dInfo=NULL);

	//////////////////////////
	//luo_base
	//////////////////////////
	virtual void initialise();
	virtual void initialiseLHS();
	virtual void writeData();
	virtual void writeCommon();
	virtual void shutDown();
	virtual void stepTime();
	virtual void _pre_step();
	virtual void _intermediate_velocity();
	virtual void _pressure();
	virtual void _update_body();
	virtual void _post_step();
	void moveBody();
	void updateSolver();

	//////////////////////////
	//projectVelocity
	//////////////////////////
	virtual void _project_velocity();

	virtual void cast();
};
