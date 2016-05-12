/***************************************************************************//**
 * \file checkTags.inl
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief checks to see if tags are coincident with body edge
 */

void NavierStokesSolver::checkPoints()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;	

	int  I;
	int  bottom, top, left, right;
	double x, y;
	bool flag;
	int  k, l;
	int  totalPoints = B.totalPoints;
	double eps = 1.e-10;
	
	// tag points at which u is evaluated
	for(int j=B.startJ0[0]-2; j<B.startJ0[0]+B.numCellsY0[0]+2; j++)
	{
		std::cout<<j<<std::endl;
		for(int i=B.startI0[0]-2; i<B.startI0[0]+B.numCellsX0[0]+2; i++)
		{
			// index of the current point on the u-grid
			//std::cout<< j*(nx-1)+i <<std::endl;
			
			// initial indices of the points on the body that define the segment under consideration
			k = totalPoints-1;
			l = 0;
			flag = false;
			
			// cycle through all the segments on the body surface
			while(l<totalPoints && !flag)
			{
				// figure out which of the two end points of the segment are at the bottom and the left
				if (B.y[k] > B.y[l])
				{
					bottom = l;
					top = k;
				}
				else
				{
					bottom = k;
					top = l;
				}
				if (B.x[k] > B.x[l])
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
				if (B.y[bottom]-eps < domInfo->yu[j] && B.y[top]-eps > domInfo->yu[j])
				{
					// if the segment is not parallel to the x-direction
					if (fabs(B.y[l]-B.y[k]) > eps)
					{
						// calculate the point of intersection of the doB.uBle ray with the boundary
						x = B.x[k] + (B.x[l]-B.x[k]) * (domInfo->yu[j]-B.y[k])/(B.y[l]-B.y[k]);;
					
						// if the point of intersection coincides with the grid point
						if (fabs(x-domInfo->xu[i]) < eps)
							std::cout<<"tagPoints warning\n";
							flag = true;
					}
				}
				// consider rays along the y-direction
				// if the ray intersects the boundary segment (right endpoint must be strictly to the right of ray; left can be on or to the left of the ray)
				if ( (B.x[left]-eps < domInfo->xu[i]) && (B.x[right]-eps > domInfo->xu[i]) && ( !flag ) ) // no need to do this part if flag is false
				{
					// if the segment is not parallel to the y-direction
					if (fabs(B.x[l]-B.x[k]) > eps)
					{
						// calculate the point of intersection of the doB.uBle ray with the boundary
						y = B.y[k] + (B.y[l]-B.y[k]) * (domInfo->xu[i]-B.x[k]) / (B.x[l]-B.x[k]);
						
						// if the point of intersection coincides with the grid point
						if (fabs(y-domInfo->yu[j]) < eps)
						{
							std::cout<<"tagPoints warning\n";
							flag = true;
						}
					}
				}
				k = l;
				l = l+1;
			}
		}
	}
std::cout<<"xdonzeo\n";
	// tag points at which v is evaluated
	for(int j=B.startJ0[0]-2; j<B.startJ0[0]+B.numCellsY0[0]+2; j++)
	{
		std::cout<<j<<std::endl;
		for(int i=B.startI0[0]-2; i<B.startI0[0]+B.numCellsX0[0]+2; i++)
		{
			// index of the current point on the u-grid
			I = j*nx+i + (nx-1)*ny;

			// initial indices of the points on the body that define the segment under consideration
			k = totalPoints-1;
			l = 0;
			
			while(l<totalPoints)
			{
				if (B.y[k] > B.y[l])
				{
					bottom = l;
					top = k;
				}
				else
				{
					bottom = k;
					top = l;
				}
				if (B.x[k] > B.x[l])
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
				/**
				* if the ray intersects the boundary segment
				* top endpoint must be strictly above the ray
				* bottom can be on or below the ray
				*/
				if (B.y[bottom]-eps < domInfo->yv[j] && B.y[top]-eps > domInfo->yv[j] && !flag)
				{
					// if the segment is not parallel to the ray
					if (fabs(B.y[l]-B.y[k]) > eps)
					{
						// calculate the point of intersection of the doB.uBle ray with the boundary
						x = B.x[k] + (B.x[l]-B.x[k]) * (domInfo->yv[j]-B.y[k])/(B.y[l]-B.y[k]);

						// if the point of intersection coincides with the grid point
						if (fabs(x-domInfo->xv[i]) < eps)
						{
							std::cout<<"tagPoints warning\n";
							flag = true;
						}
					}
				}
				// consider rays along the y-direction
				if (B.x[left]-eps < domInfo->xv[i] && B.x[right]-eps > domInfo->xv[i] && !flag)
				{
					// if the segment is not parallel to the ray
					if (fabs(B.x[l]-B.x[k]) > eps)
					{
						// calculate the point of intersection of the doB.uBle ray with the boundary
						y = B.y[k] + (B.y[l]-B.y[k]) * (domInfo->xv[i]-B.x[k])/(B.x[l]-B.x[k]);

						// if the point of intersection coincides with the grid point
						if (fabs(y-domInfo->yv[j]) < eps)
						{
							std::cout<<"tagPoints warning\n";
							flag = true;
						}
					}
				}
				k = l;
				l = l+1;
			}
		}
	}
	std::cout<<"donzo\n";
}