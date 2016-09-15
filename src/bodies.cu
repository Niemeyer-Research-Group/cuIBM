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
	std::vector<body> *B = db["flow"]["bodies"].get<std::vector<body> *>();
	// number of bodies in the flow
	numBodies = B->size();

	//oscylinder
	xCoeff = (*B)[0].xCoefficient;
	yCoeff = (*B)[0].yCoefficient;
	uCoeff = (*B)[0].uCoefficient;
	vCoeff = (*B)[0].vCoefficient;
	xfrequency = (*B)[0].xfrequency;
	yfrequency = (*B)[0].yfrequency;
	xPhase = (*B)[0].xPhase;
	yPhase = (*B)[0].yPhase;
	uPhase = (*B)[0].uPhase;
	vPhase = (*B)[0].vPhase;

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
	
	xleft.resize(numBodies);
	xright.resize(numBodies);
	ytop.resize(numBodies);
	ybot.resize(numBodies);


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
	xk.resize(totalPoints);
	yk.resize(totalPoints);
	uB.resize(totalPoints);
	vB.resize(totalPoints);
	uBk.resize(totalPoints);
	vBk.resize(totalPoints);

	force_dudn.resize(totalPoints);
	force_dvdn.resize(totalPoints);
	force_pressure.resize(totalPoints);
	force_x.resize(totalPoints);
	force_y.resize(totalPoints);
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
	cast();
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
		calculateTightBoundingBoxes(db, D);
		calculateBoundingBoxes(db, D);
		numCellsXHost = numCellsX[0];
		numCellsYHost = numCellsY[0];
	}

	midX=0;
	midY=0;
	midX0=0;
	midY0=0;
	midXk=0;
	midYk=0;
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
	centerVelocityVk = 0;
	centerVelocityUk = 0;
	centerVelocityU0= 0;
	centerVelocityV0= 0;
}

void bodies::cast()
{
	numPoints_r		= thrust::raw_pointer_cast ( &(numPoints[0]));
	offsets_r		= thrust::raw_pointer_cast ( &(offsets[0]));
	startI_r		= thrust::raw_pointer_cast ( &(startI[0]));
	startJ_r		= thrust::raw_pointer_cast ( &(startJ[0]));
	numCellsX_r		= thrust::raw_pointer_cast ( &(numCellsX[0]));
	numCellsY_r		= thrust::raw_pointer_cast ( &(numCellsY[0]));
	startI0_r		= thrust::raw_pointer_cast ( &(startI0[0]));
	startJ0_r		= thrust::raw_pointer_cast ( &(startJ0[0]));
	numCellsX0_r	= thrust::raw_pointer_cast ( &(numCellsX0[0]));
	numCellsY0_r	= thrust::raw_pointer_cast ( &(numCellsY0[0]));

	xmin_r			= thrust::raw_pointer_cast ( &(xmin[0]));
	xmax_r			= thrust::raw_pointer_cast ( &(xmax[0]));
	ymin_r			= thrust::raw_pointer_cast ( &(ymin[0]));
	ymax_r			= thrust::raw_pointer_cast ( &(ymax[0]));
	xmin0_r			= thrust::raw_pointer_cast ( &(xmin0[0]));
	xmax0_r			= thrust::raw_pointer_cast ( &(xmax0[0]));
	ymin0_r			= thrust::raw_pointer_cast ( &(ymin0[0]));
	ymax0_r			= thrust::raw_pointer_cast ( &(xmax0[0]));
	X_r				= thrust::raw_pointer_cast ( &(X[0]));
	Y_r				= thrust::raw_pointer_cast ( &(Y[0]));
	ds_r			= thrust::raw_pointer_cast ( &(ds[0]));
	ones_r			= thrust::raw_pointer_cast ( &(ones[0]));
	x_r				= thrust::raw_pointer_cast ( &(x[0]));
	y_r				= thrust::raw_pointer_cast ( &(y[0]));
	//xk_r			= thrust::raw_pointer_cast ( &(xk[0]));
	//yk_r			= thrust::raw_pointer_cast ( &(yk[0]));
	uB_r			= thrust::raw_pointer_cast ( &(uB[0]));
	vB_r			= thrust::raw_pointer_cast ( &(vB[0]));
	uBk_r			= thrust::raw_pointer_cast ( &(uBk[0]));
	vBk_r			= thrust::raw_pointer_cast ( &(vBk[0]));
	xleft_r			= thrust::raw_pointer_cast ( &(xleft[0]));
	xright_r		= thrust::raw_pointer_cast ( &(xright[0]));
	ybot_r			= thrust::raw_pointer_cast ( &(ybot[0]));
	ytop_r			= thrust::raw_pointer_cast ( &(ytop[0]));
	test_r			= thrust::raw_pointer_cast ( &(test[0]));
	x1_r			= thrust::raw_pointer_cast ( &(x1[0]));
	x2_r			= thrust::raw_pointer_cast ( &(x2[0]));
	x3_r			= thrust::raw_pointer_cast ( &(x3[0]));
	x4_r			= thrust::raw_pointer_cast ( &(x4[0]));
	y1_r			= thrust::raw_pointer_cast ( &(y1[0]));
	y2_r			= thrust::raw_pointer_cast ( &(y2[0]));
	y3_r			= thrust::raw_pointer_cast ( &(y3[0]));
	y4_r			= thrust::raw_pointer_cast ( &(y4[0]));
	q1_r			= thrust::raw_pointer_cast ( &(q1[0]));
	q2_r			= thrust::raw_pointer_cast ( &(q2[0]));
	q3_r			= thrust::raw_pointer_cast ( &(q3[0]));
	q4_r			= thrust::raw_pointer_cast ( &(q4[0]));
	point_x_r		= thrust::raw_pointer_cast ( &(point_x[0]));
	point_y_r		= thrust::raw_pointer_cast ( &(point_y[0]));
	point2_x_r		= thrust::raw_pointer_cast ( &(point2_x[0]));
	point2_y_r		= thrust::raw_pointer_cast ( &(point2_y[0]));
	point3_x_r		= thrust::raw_pointer_cast ( &(point3_x[0]));
	point3_y_r		= thrust::raw_pointer_cast ( &(point3_y[0]));
	force_pressure_r= thrust::raw_pointer_cast ( &(force_pressure[0]));
	force_dudn_r	= thrust::raw_pointer_cast ( &(force_dudn[0]));
	force_dvdn_r	= thrust::raw_pointer_cast ( &(force_dvdn[0]));
	force_x_r		= thrust::raw_pointer_cast ( &(force_x[0]));
	force_y_r		= thrust::raw_pointer_cast ( &(force_y[0]));
}

//flag this isn't setup to work with multiple bodies
//flag this kernel is setup to be called recursivly to handle body sizes larger than than
//flag this kernel isn't working for smaller body node spacing (0.005 works, 0.004 does not), disabling for now
//the maximum number of points allowable in a block
__global__
void boundingBox(double *x, double *y,
				 thrust::device_vector<double>::iterator xmax_in, thrust::device_vector<double>::iterator xmin_in, thrust::device_vector<double>::iterator ymax_in, thrust::device_vector<double>::iterator ymin_in,
				 int *startI, int *startJ, int *numCellsX, int *numCellsY,
				 double *xmax, double *xmin, double *ymax, double *ymin, double scale)
{
	if (threadIdx.x > 0)
		return;

	xmax[0] = *xmax_in;
	xmin[0] = *xmin_in;
	ymax[0] = *ymax_in;
	ymin[0] = *ymin_in;
	double	dx = xmax[0]-xmin[0],
			dy = ymax[0]-ymin[0];
	xmax[0] += 0.5*dx*(scale-1.0);
	xmin[0] -= 0.5*dx*(scale-1.0);
	ymax[0] += 0.5*dy*(scale-1.0);
	ymin[0] -= 0.5*dy*(scale-1.0);
	int i=0,
		j=0;
	while(x[i+1] < xmin[0])
		i++;
	while(y[j+1] < ymin[0])
		j++;
	startI[0] = i;
	startJ[0] = j;

	while(x[i] < xmax[0])
		i++;
	while(y[j] < ymax[0])
		j++;
	numCellsX[0] = i - startI[0];
	numCellsY[0] = j - startJ[0];
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
	double scale = db["simulation"]["scaleCV"].get<double>();
	double	*x_r = thrust::raw_pointer_cast( &(D.x[0]) ),
			*y_r = thrust::raw_pointer_cast( &(D.y[0]) );

	thrust::device_vector<double>::iterator iter_xmax,
											iter_xmin,
											iter_ymax,
											iter_ymin;

	iter_xmax = thrust::max_element(x.begin(),x.end());
	iter_xmin = thrust::min_element(x.begin(),x.end());
	iter_ymax = thrust::max_element(y.begin(),y.end());
	iter_ymin = thrust::min_element(y.begin(),y.end());


	const int blocksize = 1;
	dim3 grid(1, 1);
	dim3 block(blocksize, 1);
	boundingBox<<<grid,block>>>(x_r,y_r,
								iter_xmax, iter_xmin, iter_ymax, iter_ymin,
								startI_r, startJ_r, numCellsX_r, numCellsY_r,
								xmax_r,xmin_r,ymax_r,ymin_r, scale);
}

void bodies::calculateTightBoundingBoxes(parameterDB &db, domain &D) //flag this should be merged into the normal calculate bounding box function
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
