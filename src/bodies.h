/***************************************************************************//**
 * \file bodies.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Declaration of the class \c bodies.
 */
#pragma once

#include "domain.h"
#include "parameterDB.h"
#include "body.h"

/**
 * \class bodies
 * \brief Contains information about bodies in the flow.
 */
class bodies
{
public:
	int  numBodies,   ///< number of bodies
	     totalPoints, ///< total number of boundary points (all bodies)
	     numCellsXHost,	///< number of cells in the x-direction in the bounding box of a body on the host
	     numCellsYHost;

	bool bodiesMove;  ///< tells whether the body is moving or not

	cusp::array1d<int, cusp::device_memory>
		numPoints,    ///< number of points for each body
		offsets;      ///< array index of the first point of each body

	cusp::array1d<int, cusp::device_memory>
		startI,       ///< starting cell index of the bounding box of a body
		startJ,       ///< starting cell index of the bounding box of a body
		numCellsX,    ///< number of cells in the x-direction in the bounding box of a body
		numCellsY,    ///< number of cells in the y-direction in the bounding box of a body
		startI0,		///< starting cell index of the bounding box of a body (original size)
		startJ0,		///< starting cell index of the bounding box of a body (og size)
		numCellsX0,		///< number of cells in the x-direction in the bounding box of a body
		numCellsY0;		///< number of cells in the y-direction in the bounding box of a body

	cusp::array1d<double, cusp::device_memory>
		xmin,  ///< lowest x-coordinate for the bounding box of a body
		xmax,  ///< highest x-coordinate for the bounding box of a body
		ymin,  ///< lowest y-coordinate for the bounding box of a body
		ymax,  ///< highest y-coordinate for the bounding box of a body
		xmin0,  ///< lowest x-coordinate for the bounding box of a body (original size)
		xmax0,  ///< highest x-coordinate for the bounding box of a body (original size)
		ymin0,  ///< lowest y-coordinate for the bounding box of a body (original size)
		ymax0;  ///< highest y-coordinate for the bounding box of a body (original size)

	cusp::array1d<double, cusp::device_memory>
		X,     ///< reference x-coordinates of the boundary points
		Y,     ///< reference y-coordinates of the boundary points
		ds,    ///< vector containing the lengths of all the boundary segments
		ones,  ///< vector of size \link totalPoints \endlink with all elements 1
		x,     ///< actual x-coordinate of the boundary points
		y,     ///< actual y-coordinate of the boundary points
		uB,    ///< x-velocity of the boundary points
		vB;    ///< y-velocity of the boundary points

	cusp::array1d<double, cusp::device_memory>
		uBk,	//x-velocity of boundary points at substep k
		vBk;	//y-velocity of boundary points at substep k

	cusp::array1d<double, cusp::device_memory>
		xleft,	///< min and max values for the body nodes
		xright,
		ybot,
		ytop;

	cusp::array1d<double, cusp::device_memory>
		test,
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
		point_x,
		point_y,
		point2_x,
		point2_y,
		point3_x,
		point3_y;

	cusp::array1d<double, cusp::device_memory>
		force_pressure,
		force_dudn,
		force_dvdn,
		force_x,
		force_y;

	double	centerVelocityU,
			centerVelocityV,
			centerVelocityV0,
			centerVelocityU0,
			midY,
			midY0,
			midX0,
			midX,
			forceX,		///< force acting on a body in the x-direction
			forceY,		///< force acting on a body in the y-direction
			xCoeff,
			uCoeff,
			frequency,
			xPhase,
			uPhase;

	// set initial position and velocity of each body
	void initialise(parameterDB &db, domain &D);

	// store indices of the bounding box of each body
	void calculateBoundingBoxes(parameterDB &db, domain &D);

	// store indices of the non scaled bounding box of each body
	void calculateTightBoundingBoxes(parameterDB &db, domain &D);

	// update position, velocity and neighbors of each body
	void update(parameterDB &db, domain &D, double Time);

	// write body coordinates into a file
	void writeToFile(std::string &caseFolder, int timeStep);

	// write body coordinates into a file called \a bodies
	void writeToFile(double *bx, double *by, std::string &caseFolder, int timeStep);
};

__global__
void boundingBox(double *x, double *y,
				 thrust::device_vector<double>::iterator xmax_in, thrust::device_vector<double>::iterator xmin_in, thrust::device_vector<double>::iterator ymax_in, thrust::device_vector<double>::iterator ymin_in,
				 int *startI, int *startJ, int *numCellsX, int *numCellsY,
				 double *xmax, double *xmin, double *ymax, double *ymin, double scale);
