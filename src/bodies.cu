/***************************************************************************//**
 * \file bodies.cu
 * \author Anush Krishnan (anush@bu.edu)
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief Implementation of the methods of the class \c bodies.
 */


#include "bodies.h"
#include <cusp/blas/blas.h>
#include <iomanip>
#include <fstream>

/**
 * \brief Sets initial position and velocity of each body.
 *
 * \param db database that contains all the simulation parameters
 * \param D information about the computational grid
 */
void bodies::initialise(parameterDB &db, domain &D)
{
	std::cout << "Initialising bodies... ";
	std::vector<body> *B = db["flow"]["bodies"].get<std::vector<body> *>();

	// number of bodies in the flow
	numBodies = B->size();

	// set the sizes of all the arrays
	numPoints.resize(numBodies);
	offsets.resize(numBodies);

	startI.resize(numBodies);
	startJ.resize(numBodies);
	numCellsX.resize(numBodies);
	numCellsY.resize(numBodies);
	startI0.resize(numBodies);
	startJ0.resize(numBodies);
	numCellsX0.resize(numBodies);
	numCellsY0.resize(numBodies);

	xmin.resize(numBodies);
	xmax.resize(numBodies);
	ymin.resize(numBodies);
	ymax.resize(numBodies);
	xmin0.resize(numBodies);
	xmax0.resize(numBodies);
	ymin0.resize(numBodies);
	ymax0.resize(numBodies);
	

	// calculate offsets, number of points in each body and the total number of points
	totalPoints = 0;
	for(int k=0; k<numBodies; k++)
	{
		offsets[k] = totalPoints;
		numPoints[k] = (*B)[k].numPoints;
		totalPoints += numPoints[k];
	}

	// fill up coordinates of body points
	X.resize(totalPoints);
	Y.resize(totalPoints);
	ds.resize(totalPoints);
	ones.resize(totalPoints);
	cusp::blas::fill(ones, 1.0);
	for(int k=0; k<numBodies; k++)
	{
		for(int i=0; i<numPoints[k]; i++)
		{
			X[i+offsets[k]] = (*B)[k].X[i];
			Y[i+offsets[k]] = (*B)[k].Y[i];
		}
	}
	x.resize(totalPoints);
	y.resize(totalPoints);
	uB.resize(totalPoints);
	vB.resize(totalPoints);
	uBk.resize(totalPoints);
	vBk.resize(totalPoints);
	I.resize(totalPoints);
	J.resize(totalPoints);

	force_pressure.resize(totalPoints);
	x1.resize(totalPoints);
	x2.resize(totalPoints);
	x3.resize(totalPoints);
	x4.resize(totalPoints);
	y1.resize(totalPoints);
	y2.resize(totalPoints);
	y3.resize(totalPoints);
	y4.resize(totalPoints);
	q1.resize(totalPoints);
	q2.resize(totalPoints);
	q3.resize(totalPoints);
	q4.resize(totalPoints);
	point_y.resize(totalPoints);
	point_x.resize(totalPoints);
	point2_y.resize(totalPoints);
	point2_x.resize(totalPoints);
	point3_y.resize(totalPoints);
	point3_x.resize(totalPoints);
	centerVelocityU = 0;
	centerVelocityV = 0;

	cusp::blas::fill(vB, 0);
	cusp::blas::fill(uB, 0);
	cusp::blas::fill(vBk, 0);
	cusp::blas::fill(uBk, 0);

	bodiesMove = false;
	for(int k=0; k<numBodies; k++)
	{
		// assume a closed body (closed loop)
		for(int i=offsets[k], j = offsets[k]+numPoints[k]-1; i<offsets[k]+numPoints[k];)
		{
			// calculate the lengths of the boundary segments
			ds[i] = sqrt( (X[i]-X[j])*(X[i]-X[j]) + (Y[i]-Y[j])*(Y[i]-Y[j]) );

			// j takes the value of i, then i is incremented
			j = i++;
		}
		// if the body is moving, set bodiesMove to true
		bodiesMove = bodiesMove || (*B)[k].moving[0] || (*B)[k].moving[1];
	}
	// set initial position of the body
	update(db, D, 0.0);

	if(numBodies)
	{
		calculateCellIndices(D);
		calculateTightBoundingBoxes(db, D);
		calculateBoundingBoxes(db, D);
	}

	midX=0;
	midY=0;
	midX0=0;
	midY0=0;
	for (int i=0;i<totalPoints;i++)
	{
		midX += x[i];
		midY += y[i];
	}
	midX /= totalPoints;
	midY /= totalPoints;
	midX=midX0;
	midY=midY0;
	centerVelocityV = 0;
	centerVelocityU = 0;
	centerVelocityU0= 0;
	centerVelocityV0= 0;
}

/**
 * \brief Stores index of each cell that contains a boundary point.
 *
 * It calculates the index of the x-coordinate and the index of the y-coordinate
 * of the bottom-left node of each cell that contains a boundary point.
 * This information is useful when transferring data between the boundary points
 * and the computational grid.
 *
 * \param D information about the computational grid
 */

void bodies::calculateCellIndices(domain &D)
{
	int	i=0, j=0;

	// find the cell for the zeroth point
	while(D.x[i+1] < x[0])
		i++;
	while(D.y[j+1] < y[0])
		j++;
	I[0] = i;
	J[0] = j;

	for(int k=1; k<totalPoints; k++)
	{
		// if the next boundary point is to the left of the current boundary point
		if(x[k] < x[k-1])
		{
			while(D.x[i] > x[k])
				i--;
		}
		// if the next boundary point is to the right of the current boundary point
		else
		{
			while(D.x[i+1] < x[k])
				i++;
		}
		// if the next boundary point is below the current boundary point
		if(y[k] < y[k-1])
		{
			while(D.y[j] > y[k])
				j--;
		}
		// if the next boundary point is above the current boundary point
		else
		{
			while(D.y[j+1] < y[k])
				j++;
		}
		I[k] = i;
		J[k] = j;
	}
}

/**
 * \brief Calculates indices of the bounding box of each body in the flow.
 *
 * First the bounding box is scaled by a coefficient stored in the database.
 * Then, indices of the x-coordinate and y-coordinate of the bottom left cell
 * of the bounding box are stored. Finally, the number of cells in the x- and y-
 * directions are calculated.
 *
 * \param db database that contains all the simulation parameters
 * \param D information about the computational grid
 */
void bodies::calculateBoundingBoxes(parameterDB &db, domain &D)
{
	double scale = db["simulation"]["scaleCV"].get<double>(),
	     dx, dy;
	int  i, j;
	for(int k=0; k<numBodies; k++)
	{
		xmin[k] = x[offsets[k]];
		xmax[k] = xmin[k];
		ymin[k] = y[offsets[k]];
		ymax[k] = ymin[k];
		for(int l=offsets[k]+1; l<offsets[k]+numPoints[k]; l++)
		{
			if(x[l] < xmin[k]) xmin[k] = x[l];
			if(x[l] > xmax[k]) xmax[k] = x[l];
			if(y[l] < ymin[k]) ymin[k] = y[l];
			if(y[l] > ymax[k]) ymax[k] = y[l];
		}
		dx = xmax[k]-xmin[k];
		dy = ymax[k]-ymin[k];
		xmax[k] += 0.5*dx*(scale-1.0);
		xmin[k] -= 0.5*dx*(scale-1.0);
		ymax[k] += 0.5*dy*(scale-1.0);
		ymin[k] -= 0.5*dy*(scale-1.0);
		
		i=0; j=0;
		while(D.x[i+1] < xmin[k])
			i++;
		while(D.y[j+1] < ymin[k])
			j++;
		startI[k] = i;
		startJ[k] = j;
		
		while(D.x[i] < xmax[k])
			i++;
		while(D.y[j] < ymax[k])
			j++;
		numCellsX[k] = i - startI[k];
		numCellsY[k] = j - startJ[k];
	}
}

void bodies::calculateTightBoundingBoxes(parameterDB &db, domain &D)
{
	double scale = db["simulation"]["scaleCV"].get<double>();
	int  i, j;
	for(int k=0; k<numBodies; k++)
	{
		xmin0[k] = x[offsets[k]];
		xmax0[k] = xmin[k];
		ymin0[k] = y[offsets[k]];
		ymax0[k] = ymin[k];
		for(int l=offsets[k]+1; l<offsets[k]+numPoints[k]; l++)
		{
			if(x[l] < xmin0[k]) xmin0[k] = x[l];
			if(x[l] > xmax0[k]) xmax0[k] = x[l];
			if(y[l] < ymin0[k]) ymin0[k] = y[l];
			if(y[l] > ymax0[k]) ymax0[k] = y[l];
		}

		i=0; j=0;
		while(D.x[i+1] < xmin0[k])
			i++;
		while(D.y[j+1] < ymin0[k])
			j++;
		startI0[k] = i;
		startJ0[k] = j;

		while(D.x[i] < xmax[k])
			i++;
		while(D.y[j] < ymax[k])
			j++;
		numCellsX0[k] = i - startI0[k];
		numCellsY0[k] = j - startJ0[k];
	}
}

/**
 * \brief Updates position, velocity and neighbors of each body.
 *
 * This is done using the formulae:
 *
 * \f$ x_{i,m} = X^c_m + (X_{i,m} - X^0_m) \cos\theta - (Y_{i,m} - Y^0_m) \sin\theta \f$
 *
 * and
 *
 * \f$ y_{i,m} = Y^c_m + (X_{i,m} - X^0_m) \sin\theta + (Y_{i,m} - Y^0_m) \cos\theta \f$
 *
 * \param db database that contains all the simulation parameters
 * \param D information about the computational grid
 * \param Time the time
 */
void bodies::update(parameterDB &db, domain &D, double Time)
{
	typedef typename cusp::array1d<double, cusp::device_memory> Array;
	typedef typename Array::iterator                 Iterator;
	typedef cusp::array1d_view<Iterator>             View;

	// views of the vectors that store the coordinates and velocities of all the body points
	View    XView, YView, xView, yView, onesView, uBView, vBView;

	// body data
	std::vector<body> *B = db["flow"]["bodies"].get<std::vector<body> *>();

	for(int l=0; l<numBodies; l++)
	{
		// update the location and velocity of the body
		(*B)[l].update(Time);

		// create the views for the current body
		if(l < numBodies-1)
		{
			XView    = View(X.begin()+offsets[l], X.begin()+offsets[l+1]);
			YView    = View(Y.begin()+offsets[l], Y.begin()+offsets[l+1]);
			onesView = View(ones.begin()+offsets[l], ones.begin()+offsets[l+1]);
			uBView   = View(uB.begin()+offsets[l], uB.begin()+offsets[l+1]);
			vBView   = View(vB.begin()+offsets[l], vB.begin()+offsets[l+1]);
			xView    = View(x.begin()+offsets[l], x.begin()+offsets[l+1]);
			yView    = View(y.begin()+offsets[l], y.begin()+offsets[l+1]);
		}
		else
		{
			XView    = View(X.begin()+offsets[l], X.end());
			YView    = View(Y.begin()+offsets[l], Y.end());
			onesView = View(ones.begin()+offsets[l], ones.end());
			xView    = View(x.begin()+offsets[l], x.end());
			yView    = View(y.begin()+offsets[l], y.end());
			uBView   = View(uB.begin()+offsets[l], uB.end());
			vBView   = View(vB.begin()+offsets[l], vB.end());

		}

		// update postitions
		// x-coordinates
		cusp::blas::axpbypcz( onesView, XView, onesView, xView, (*B)[l].Xc[0],  cos((*B)[l].Theta), -(*B)[l].X0[0]*cos((*B)[l].Theta) );
		cusp::blas::axpbypcz( xView,    YView, onesView, xView,           1.0, -sin((*B)[l].Theta),  (*B)[l].X0[1]*sin((*B)[l].Theta) );
		// y-coordinates
		cusp::blas::axpbypcz( onesView, XView, onesView, yView, (*B)[l].Xc[1],  sin((*B)[l].Theta), -(*B)[l].X0[0]*sin((*B)[l].Theta) );
		cusp::blas::axpbypcz( yView,    YView, onesView, yView,           1.0,  cos((*B)[l].Theta), -(*B)[l].X0[1]*cos((*B)[l].Theta) );

		// update velocities
		// x-velocities
		cusp::blas::axpbypcz(onesView, yView, onesView, uBView, (*B)[l].vel[0], -(*B)[l].angVel,  (*B)[l].angVel*(*B)[l].Xc[1]);
		// y-velocities
		cusp::blas::axpbypcz(onesView, xView, onesView, vBView, (*B)[l].vel[1],  (*B)[l].angVel, -(*B)[l].angVel*(*B)[l].Xc[0]);
	}

	if(numBodies)
		calculateCellIndices(D);
}


/**
 * \brief Writes body coordinates into a file (using data from the device).
 *
 * \param caseFolder directory of the simulation
 * \param timeStep time-step of the simulation
 */
void bodies::writeToFile(std::string &caseFolder, int timeStep)
{
	cusp::array1d<double, cusp::host_memory>
		xHost = x,
		yHost = y;
	double *bx = thrust::raw_pointer_cast(&(xHost[0])),
	     *by = thrust::raw_pointer_cast(&(yHost[0]));
	writeToFile(bx, by, caseFolder, timeStep);
}

/**
 * \brief Writes body coordinates into a file called \a bodies.
 *
 * \param bx x-coordinate of all points of all bodies
 * \param by y-coordinate of all points of all bodies
 * \param caseFolder directory of the simulation
 * \param timeStep time-step of the simulation
 */
void bodies::writeToFile(double *bx, double *by, std::string &caseFolder, int timeStep)
{
	std::string       path;
	std::stringstream out;
	out << caseFolder << '/' << std::setfill('0') << std::setw(7) << timeStep << "/bodies";
	std::ofstream file(out.str().c_str());;
	file << '#' << std::setw(19) << "x-coordinate" << std::setw(20) << "y-coordinate" << std::endl;
	for (int l=0; l < totalPoints; l++)
	{
		file << bx[l] << '\t' << by[l] << '\n';
	}
	file.close();
}
