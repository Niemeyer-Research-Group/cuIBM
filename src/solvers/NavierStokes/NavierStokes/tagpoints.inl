/***************************************************************************//**
 * \file tagPoints.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c DirectForcingSolver to tag
 *        points near the immersed boundary using a ray-tracing algorithm.
 */

/**
 * \brief Tags the forcing nodes among the velocity nodes, i.e. the nodes at 
 *        which the velocity interpolation is performed.
 */

void NavierStokesSolver::tagPoints()
{
	logger.startTimer("tagPoints");
	// transferring boundary point coordinates to the host
	cusp::array1d<double, cusp::host_memory> bxH(B.totalPoints), byH(B.totalPoints), uBH(B.totalPoints), vBH(B.totalPoints);
	bxH = B.x;
	byH = B.y;
	uBH = B.uB;
	vBH = B.vB;

	// creating raw pointers
	double	*bx = thrust::raw_pointer_cast(&(bxH[0])),
			*by = thrust::raw_pointer_cast(&(byH[0])),
			*uB = thrust::raw_pointer_cast(&(uBH[0])),
			*vB = thrust::raw_pointer_cast(&(vBH[0]));

	tagPoints(bx, by, uB, vB);

	// transferring tag and a data to the device
	tagsD    = tags;
	tags2D   = tags2;
	tagsInD	 = tagsIn;
	tagsPD	 = tagsP;
	tagsPOutD= tagsPOut;
	aD       = a;
	bD       = b;
	distance_from_u_to_bodyD = distance_from_u_to_body;
	distance_from_v_to_bodyD = distance_from_v_to_body;
	uvD      = uv;

	logger.stopTimer("tagPoints");
}

// Bilinear Fadlun1c-type interpolation outside the body, for a moving body.
/**
 * \brief Tags all the forcing nodes required for the type of linear
 *        interpolation explained in the paper by Fadlun et al. (2000).
 *
 * It uses a raytracing algorithm to detect points that are near the boundary,
 * and just outside it. For more information about the algorithm, read the 
 * section on ray-crossings in the Search and Intersection chapter of the 
 * book Computational Geometry in C by Joseph O'Rourke.
 *
 * \param bx host array of the x-coordinates of the boundary points
 * \param by host array of the y-coordinates of the boundary points
 * \param uB host array of the x-components of the boundary velocities
 * \param vB host array of the y-components of the boundary velocities
 */

void NavierStokesSolver::tagPoints(double *bx, double *by, double *uB, double *vB)
{
	int  nx = NavierStokesSolver::domInfo->nx,
	     ny = NavierStokesSolver::domInfo->ny;
	
	double	*xu = thrust::raw_pointer_cast(&(NavierStokesSolver::domInfo->xu[0])),
			*yu = thrust::raw_pointer_cast(&(NavierStokesSolver::domInfo->yu[0]));

	parameterDB &db = *NavierStokesSolver::paramDB;
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();

	int  I, ip;
	int  bottom, top, left, right;
	double eps = 1.e-10;
	bool outsideX, outsideY;
	int  bdryFlagX, bdryFlagY, bdryFlag2X, bdryFlag2Y;
	double Xa, Xb, Ya, Yb;
	double uvX, uvY;
	bool flag;
	double x, y;
	double midX, midY;
	int  k, l;
	int  totalPoints = B.totalPoints;
	cusp::blas::fill(tagsIn,-1);
	//cusp::blas::fill(distance_from_u_to_body, 1);
	//cusp::blas::fill(distance_from_v_to_body, 1);

	// tags pressure nodes inside the body
	for (int i=0;i<totalPoints;i++)
	{
		midX += bx[i];
		midY += by[i];
	}
	midX /= totalPoints;
	midY /= totalPoints;

	// tag points at which u is evaluated
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			// index of the current point on the u-grid
			I = j*(nx-1)+i;
			ip = j*nx + i;
			
			// tags and coefficients
			tags[I] = -1;
			tags2[I] = -1;
			a[I] = 1.0;
			b[I] = 1.0;
			
			// initial indices of the points on the body that define the segment under consideration
			k = totalPoints-1;
			l = 0;
			
			outsideX = true;
			outsideY = true;
			bdryFlagX = -1;  // stores if a point is near the boundary
			bdryFlagY = -1;
			bdryFlag2X = -1;
			bdryFlag2Y = -1;
			uvX = 0.0;
			uvY = 0.0;
			Xa = 1.0;
			Ya = 1.0;
			Xb = 1.0;
			Yb = 1.0;
			flag = false;
			// cycle through all the segments on the body surface
			while(l<totalPoints && !flag)
			{
				// figure out which of the two end points of the segment are at the bottom and the left
				if (by[k] > by[l])
				{
					bottom = l;
					top = k;
				}
				else
				{
					bottom = k;
					top = l;
				}
				if (bx[k] > bx[l])
				{
					left = l;
					right = k;
				}
				else
				{
					left = k;
					right = l;
				}

				// consider rays along the x-direction
				// if the ray intersects the boundary segment (top endpoint must be strictly above the ray; bottom can be on or below the ray)
				if (by[bottom]-eps < yu[j] && by[top]-eps > yu[j])
				{
					// if the segment is not parallel to the x-direction 
					if (fabs(by[l]-by[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						x = bx[k] + (bx[l]-bx[k]) * (yu[j]-by[k])/(by[l]-by[k]);
						
						// calculate the body velocity at the point of intersection
						uvX = uB[k] + (uB[l]-uB[k]) * (yu[j]-by[k])/(by[l]-by[k]);

						// if the point of intersection coincides with the grid point
						if (fabs(x-xu[i]) < eps)
						{
							outsideX  = true;
							bdryFlagX = I;
							tagsIn[I] = I;
							Xa        = 0.0;
							Xb        = 1.0;
							flag      = true; // flag is true when the point of intersection coincides with the grid point
							std::cout<<"tagpoints warning\n\n\n\n\n";
						}
						// if the point of intersection lies to the right of the grid point (right-facing ray intersects the boundary)
				 		else if (x > xu[i]+eps)
							outsideX = !outsideX;
						
						// if the point of intersection is in the cell to the immediate left of the grid point
						if (x>xu[i-1]+eps && x<xu[i]-eps)
						{
							bdryFlagX  = I;
							bdryFlag2X = I+1;
							//if (outsideX)
							if (tags[I-1]==-1)
								tagsIn[I-1]	= I-1;
							Xa = xu[i]-x;
							Xb = xu[i+1]-xu[i];
							if (x > midX)
								distance_from_u_to_body[ip] = Xa;
							//case 1
						}
						// if the point of intersection is in the cell to the immediate right of the grid point
						else if (x>xu[i]+eps && x<xu[i+1]-eps)
						{
							bdryFlagX  = I;
							bdryFlag2X = I-1;
							//if (outsideX)
							if(tags[I+1] == -1)
								tagsIn[I+1]	= I+1;
							Xa = x-xu[i];
							Xb = xu[i]-xu[i-1];
							if (x < midX)
								distance_from_u_to_body[ip+1] = Xa;
							//case 2
						}
					}
				}
				// consider rays along the y-direction
				// if the ray intersects the boundary segment (right endpoint must be strictly to the right of ray; left can be on or to the left of the ray)
				if ( (bx[left]-eps < xu[i]) && (bx[right]-eps > xu[i]) && ( !flag ) ) // no need to do this part if flag is false
				{
					// if the segment is not parallel to the y-direction
					if (fabs(bx[l]-bx[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						y = by[k] + (by[l]-by[k]) * (xu[i]-bx[k]) / (bx[l]-bx[k]);
						
						// calculate the body velocity at the point of intersection
						uvY = uB[k] + (uB[l]-uB[k]) * (xu[i]-bx[k])/(bx[l]-bx[k]);
						
						// if the point of intersection coincides with the grid point
						if (fabs(y-yu[j]) < eps)
						{
							outsideY  = true; // then the point is considered to be outside the grid
							bdryFlagY = I;    // the point is considered to be a forcing point, with index I
							bdryFlag2Y= I;
							tagsIn[I]=I;
							Ya        = 0.0;  // the coefficient for the linear interpolation during forcing
							Yb        = 1.0;
							flag      = true; // flag is true when the point of intersection coincides with the grid point
							std::cout<<"tagpoitns warning\n\n\n\n\n";
						}
						// if the point of intersection lies to the top of the grid point
				 		else if (y > yu[j]+eps)
							outsideY = !outsideY; // then flip if inside or outside (start with true, i.e. outside) //this seems wrong too

						// if point of intersection is just below the concerned grid point
						if (y>yu[j-1]+eps && y<yu[j]-eps)
						{
							bdryFlagY = I;
							bdryFlag2Y= I+(nx-1);
							//if (outsideY)
							if(tags[I-nx+1]==-1)
								tagsIn[I-(nx-1)]=I-(nx-1);
							Ya = yu[j]-y;
							Yb = yu[j+1]-yu[j];
						}
						// if point of intersection is just above the concerned grid point
						else if (y>yu[j]+eps && y<yu[j+1]-eps)
						{
							bdryFlagY = I;
							bdryFlag2Y= I-(nx-1);
							//if (outsideY)
							if (tags[I+nx-1]==-1)
								tagsIn[I+(nx-1)]=I+(nx-1);
							Ya = y-yu[j];
							Yb = yu[j]-yu[j-1];
						}
					}
				}
				k = l;
				l = l+1;
			}

			if (outsideX && bdryFlagX>=0)
			{
				tagsIn[I] = -1;
				tags[I]		= bdryFlagX;
				tags2[I]	= bdryFlag2X;
				//distance_from_u_to_body[I] = Xa;
				a[I]		= Xa;
				b[I]		= Xb;
				uv[I]		= uvX;
			}
			else if (outsideY && bdryFlagY>=0)
			{
				tagsIn[I] = -1;
				tags[I]		= bdryFlagY;
				tags2[I]	= bdryFlag2Y;
				//distance_from_u_to_body[I] = Xa;
				a[I]		= Ya;
				b[I]		= Yb;
				uv[I]		= uvY;
			}
		}
	}
	double	*xv = thrust::raw_pointer_cast(&(domInfo->xv[0])),
			*yv = thrust::raw_pointer_cast(&(domInfo->yv[0]));
	
	// tag points at which v is evaluated
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			// index of the current point on the u-grid
			I = j*nx+i + (nx-1)*ny;
			ip = j*nx + i;

			// tags and coefficients
			tags[I] = -1;
			a[I] = 1.0;

			// initial indices of the points on the body that define the segment under consideration
			k = totalPoints-1;
			l = 0;
			
			outsideX = true;
			outsideY = true;
			bdryFlagX = -1;
			bdryFlagY = -1;
			uvX = 0.0;
			uvY = 0.0;
			Xa = 1.0;
			Ya = 1.0;
			Xb = 1.0;
			Yb = 1.0;
			flag = false;
			
			while(l<totalPoints)
			{
				if (by[k] > by[l])
				{
					bottom = l;
					top = k;
				}
				else
				{
					bottom = k;
					top = l;
				}
				if (bx[k] > bx[l])
				{
					left = l;
					right = k;
				}
				else
				{
					left = k;
					right = l;
				}
				// consider rays along the x-direction
				// if the ray intersects the boundary segment top endpoint must be strictly above the ray bottom can be on or below the ray
				if (by[bottom]-eps < yv[j] && by[top]-eps > yv[j] && !flag)
				{
					// if the segment is not parallel to the ray
					if (fabs(by[l]-by[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						x = bx[k] + (bx[l]-bx[k]) * (yv[j]-by[k])/(by[l]-by[k]);
						// calculate the body velocity at the point of intersection
						uvX = vB[k] + (vB[l]-vB[k]) * (yv[j]-by[k])/(by[l]-by[k]);

						// if the point of intersection coincides with the grid point
						if (fabs(x-xv[i]) < eps)
						{
							outsideX   = true;							
							bdryFlagX = I;
							bdryFlag2X= I;
							tagsIn[I] = I;
							Xa        = 0.0;
							Xb        = 1.0;
							flag      = true;
							std::cout<<"tagpoints warning\n\n\n\n\n";
						}
						// if the point of intersection lies to the right of the grid point
				 		else if (x > xv[i]+eps)
							outsideX = !outsideX;

						// if the point of intersection is in the cell to the immediate right of the grid point
						if (x>xv[i-1]+eps && x<xv[i]-eps)
						{
							bdryFlagX  = I;
							bdryFlag2X = I+1;
							if (tags[I-1]==-1)
								tagsIn[I-1] = I-1;
							Xa = xv[i]-x;
							Xb = xv[i+1]-xv[i];
						}
						// if the point of intersection is in the cell to the immediate left of the grid point
						else if (x>xv[i]+eps && x<xv[i+1]-eps)
						{
							bdryFlagX  = I;
							bdryFlag2X = I-1;
							if (tags[I+1] == -1)
								tagsIn[I+1] = I+1;
							Xa = x-xv[i];
							Xb = xv[i]-xv[i-1];
						}
					}
				}
				// consider rays along the y-direction
				if (bx[left]-eps < xv[i] && bx[right]-eps > xv[i] && !flag)
				{
					// if the segment is not parallel to the ray
					if (fabs(bx[l]-bx[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						y = by[k] + (by[l]-by[k]) * (xv[i]-bx[k])/(bx[l]-bx[k]);
						// calculate the body velocity at the point of intersectioin
						uvY = vB[k] + (vB[l]-vB[k]) * (xv[i]-bx[k])/(bx[l]-bx[k]);

						// if the point of intersection coincides with the grid point
						if (fabs(y-yv[j]) < eps)
						{
							outsideY   = true;
							bdryFlagY = I;
							bdryFlag2Y= I;
							tagsIn[I] = I;
							Ya        = 0.0;
							Yb        = 1.0;
							flag      = true;
							std::cout<<"tagpoints warning\n\n\n\n\n";
						}
						// if the point of intersection lies to the top of the grid point
				 		else if (y > yv[j]+eps)
							outsideY = !outsideY;

						if (y>yv[j-1]+eps && y<yv[j]-eps)
						{
							bdryFlagY  = I;
							bdryFlag2Y = I+nx;
							if (tags[I-nx] == -1)
								tagsIn[I-nx] = I-nx;
							Ya = yv[j]-y;
							Yb = yv[j+1]-yv[j];
							//case 3
							if (outsideY)
								distance_from_v_to_body[ip] = Ya;
						}
						else if (y>yv[j]+eps && y<yv[j+1]-eps)
						{
							bdryFlagY  = I;
							bdryFlag2Y = I-nx;
							if(tags[I+nx] == -1)
								tagsIn[I+nx] = I+nx;
							Ya = y-yv[j];
							Yb = yv[j]-yv[j-1];
							//case 4
							if (outsideY)
								distance_from_v_to_body[ip+nx] = Ya;
						}
					}
				}
				k = l;
				l = l+1;
			}
			if (outsideY && bdryFlagY>=0)
			{
				tagsIn[I] = -1;
				tags[I]    = bdryFlagY;
				tags2[I]   = bdryFlag2Y;
				a[I]  = Ya;
				b[I] = Yb;
				//distance_from_v_to_body[I] = Ya;
				uv[I]      = uvY;
			}
			else if (outsideX && bdryFlagX>=0)
			{
				tagsIn[I] = -1;
				tags[I]    = bdryFlagX;
				tags2[I]   = bdryFlag2X;
				a[I]  = Xa;
				b[I] = Xb;
				//distance_from_v_to_body[I] = Ya;
				uv[I]      = uvX;
			}
		}
	}

	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx; i++)
		{
			// index of the current point on the u-grid
			I = j*nx+i;

			// tags
			tagsP[I] = -1;

			// initial indices of the points on the body that define the segment under consideration
			k = totalPoints-1;
			l = 0;

			while(l<totalPoints)
			{
				if (by[k] > by[l])
				{
					bottom = l;
					top = k;
				}
				else
				{
					bottom = k;
					top = l;
				}
				if (bx[k] > bx[l])
				{
					left = l;
					right = k;
				}
				else
				{
					left = k;
					right = l;
				}
				// consider rays along the x-direction
				// if the ray intersects the boundary segment top endpoint must be strictly above the ray bottom can be on or below the ray
				if (by[bottom]-eps < yu[j] && by[top]-eps > yu[j] && !flag)
				{
					// if the segment is not parallel to the ray
					if (fabs(by[l]-by[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						x = bx[k] + (bx[l]-bx[k]) * (yu[j]-by[k])/(by[l]-by[k]);

						// if the point of intersection coincides with the grid point
						if (fabs(x-xv[i]) < eps)
						{
							std::cout<<"tagpoints warning\n\n\n\n\n";
						}
						// just inside, right of mid
						if (x > midX + eps && x > xv[i] + eps && x < xv[i+1] + eps)
						{
							tagsP[I] = I;
						}
						// just inside, left of mid
						else if (x < midX + eps && x < xv[i] +eps && x > xv[i-1])
						{
							tagsP[I] = I;
						}
					}
				}
				// consider rays along the y-direction
				if (bx[left]-eps < xv[i] && bx[right]-eps > xv[i] && !flag)
				{
					// if the segment is not parallel to the ray
					if (fabs(bx[l]-bx[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						y = by[k] + (by[l]-by[k]) * (xv[i]-bx[k])/(bx[l]-bx[k]);

						// if the point of intersection coincides with the grid point
						if (fabs(y-yu[j]) < eps)
						{
							std::cout<<"tagpoints warning\n\n\n\n\n";
						}

						// just inside, north of mid
						if (y > midY + eps && y > yu[j] +eps && y < yu[j+1])
						{
							tagsP[I] = I;
						}
						//just inside, south of mid
						else if (y < midY + eps && y < yu[j] +eps && y > yu[j-1] +eps)
						{
							tagsP[I] = I;
						}
					}
				}
				k = l;
				l = l+1;
			}// end while
		}//end i
	}//end j

	// tags pressure nodes outside the body
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx; i++)
		{
			// index of the current point on the u-grid
			I = j*nx+i;

			// tags
			tagsPOut[I] = -1;

			// initial indices of the points on the body that define the segment under consideration
			k = totalPoints-1;
			l = 0;

			while(l<totalPoints)
			{
				if (by[k] > by[l])
				{
					bottom = l;
					top = k;
				}
				else
				{
					bottom = k;
					top = l;
				}
				if (bx[k] > bx[l])
				{
					left = l;
					right = k;
				}
				else
				{
					left = k;
					right = l;
				}
				// consider rays along the x-direction
				// if the ray intersects the boundary segment top endpoint must be strictly above the ray bottom can be on or below the ray
				if (by[bottom]-eps < yu[j] && by[top]-eps > yu[j] && !flag)
				{
					// if the segment is not parallel to the ray
					if (fabs(by[l]-by[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						x = bx[k] + (bx[l]-bx[k]) * (yu[j]-by[k])/(by[l]-by[k]);

						// if the point of intersection coincides with the grid point
						if (fabs(x-xv[i]) < eps)
						{
							std::cout<<"tagpoints warning\n\n\n\n\n";
						}
						// just inside, right of mid
						if (x > midX + eps && x < xv[i] + eps && x > xv[i-1] + eps)
						{
							tagsPOut[I] = I;
						}
						// just inside, left of mid
						else if (x < midX + eps && x > xv[i] +eps && x < xv[i+1] + eps)
						{
							tagsPOut[I] = I;
						}
					}
				}
				// consider rays along the y-direction
				if (bx[left]-eps < xv[i] && bx[right]-eps > xv[i] && !flag)
				{
					// if the segment is not parallel to the ray
					if (fabs(bx[l]-bx[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						y = by[k] + (by[l]-by[k]) * (xv[i]-bx[k])/(bx[l]-bx[k]);

						// if the point of intersection coincides with the grid point
						if (fabs(y-yu[j]) < eps)
						{
							std::cout<<"tagpoints warning\n\n\n\n\n";
						}

						// just inside, north of mid
						if (y > midY + eps && y < yu[j] +eps && y > yu[j-1])
						{
							tagsPOut[I] = I;
						}
						//just inside, south of mid
						else if (y < midY + eps && y > yu[j] +eps && y < yu[j+1] +eps)
						{
							tagsPOut[I] = I;
						}
					}
				}
				k = l;
				l = l+1;
			}// end while
		}//end i
	}//end j


	//set inside of pressure to 0
	for (int j=0; j<ny-1; j++)
	{
		bool rowIsntDone = true;
		bool flag = false;
		for (int i=1; i<nx-1; i++)
		{
			I = j*nx+i;

			if (tagsP[I-1] >= 0)
			{
				if (tagsP[I]!=-1 && tagsP[I-1]==0)
				{
					rowIsntDone = false;
				}
				if(tagsP[I]==-1 && rowIsntDone)
				{
					tagsP[I]=0;
				}
				if(tagsP[I]==0 && tagsP[I-1] !=0)
				{
					int k = i;
					while (k < nx)
					{
						k++;
						if (tagsP[j*nx+k] != -1 )
						{
							flag = true;
						}
					}
					if (!flag)
					{
						tagsP[I]=-1;
						rowIsntDone = false;
					}
				}
			}
		}
	}

	//set tagsin x inside to 0
	for (int j=0; j<ny-1; j++)
		{
			bool rowIsntDone = true;
			bool flag = false;
			for (int i=1; i<nx-1; i++)
			{
				I = j*(nx-1)+i;

				if (tagsIn[I-1] >= 0)
				{
					if (tagsIn[I]!=-1 && tagsIn[I-1]==0)
					{
						rowIsntDone = false;
					}
					if(tagsIn[I]==-1 && rowIsntDone)
					{
						tagsIn[I]=0;
					}
					if(tagsIn[I]==0 && tagsIn[I-1] !=0)
					{
						int k = i;
						while (k < nx-1)
						{
							k++;
							if (tagsIn[j*(nx-1)+k] != -1 )
							{
								flag = true;
							}
						}
						if (!flag)
						{
							tagsIn[I]=-1;
							rowIsntDone = false;
						}
					}
				}
			}
		}

	//set tagsin y inside to 0
	for (int j=0; j<ny-1; j++)
	{
		bool rowIsntDone = true;
		bool flag = false;
		for (int i=1; i<nx-1; i++)
		{
			I = j*(nx)+i + (nx-1)*ny;
			if (tagsIn[I-1] >= 0)
			{
				if (tagsIn[I]!=-1 && tagsIn[I-1]==0)
				{
					rowIsntDone = false;
				}
				if(tagsIn[I]==-1 && rowIsntDone)
				{
					tagsIn[I]=0;
				}
				if(tagsIn[I]==0 && tagsIn[I-1] !=0)
				{
					int k = i;
					while (k < nx-1)
					{
						k++;
						if (tagsIn[j*(nx)+k + (nx-1)*ny] != -1 )
						{
							flag = true;
						}
					}
					if (!flag)
					{
						tagsIn[I]=-1;
						rowIsntDone = false;
					}
				}
			}
		}
	}
}// end tagpoints
