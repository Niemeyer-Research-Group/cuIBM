/***************************************************************************//**
 * \file domain.h
 * \author Anush Krishnan (anush@bu.edu)
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief Definition of the class \c domain.
 */


#pragma once
#include <cusp/array1d.h>

/**
 * \class domain
 * \brief Stores information about the computational grid.
 */
class domain
{
public:
	int   nx, ///< number of cells in the x-direction
	      ny; ///< number of cells in the y-direction
	
	double mid_h; ///< dx/dy value in the middle

	cusp::array1d<double, cusp::device_memory>
		x,  ///< x-coordinates of the nodes stored on the device
		y,  ///< y-coordinates of the nodes stored on the device
		dx, ///< x- cell widths stored on the device
		dy, ///< y- cell widths stored on the device
		xu,  ///< x-coordinates where the x-components of velocity are evaluated on the device
		yu,  ///< y-coordinates where the x-components of velocity are evaluated on the device
		xv,  ///< x-coordinates where the y-components of velocity are evaluated on the device
		yv;  ///< y-coordinates where the y-components of velocity are evaluated on the device
};
