/***************************************************************************//**
 * \file projectVelocity.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \CPU Author, Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c DirectForcingSolver to tag
 *        points near the immersed boundary using a ray-tracing algorithm.
 */

#include "tagPoints.h"

namespace kernels
{
__global__
void tag_u(int *tags, int *tagsIn, int *tags2, double *bx, double *by, double *uB, double *vB, double *yu, double *xu,
			   double *distance_from_intersection_to_node, double *distance_between_nodes_at_IB, double *distance_from_u_to_body, double *distance_from_v_to_body, double *uv,
			   int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints, double midX, double midY)
{
	// calculate indicies indices
	int idx	= threadIdx.x + blockDim.x * blockIdx.x,
		i	= idx % (i_end-i_start),
		j	= idx / (i_end-i_start),
		I	= i_start + i,
		J	= j_start + j,
		iu = J*(nx-1) + I,
		ip = J*nx + I;

	// return if out of bounds of the array
	if (iu >= (nx-1)*ny)
			return;

	// initial indices of the points on the body that define the segment under consideration
	int 	k = totalPoints-1,
			l = 0;

	// logic for the segment
	bool	outsideX = true,
			outsideY = true,
			flag = false;
	int		bdryFlagX = -1,  // stores if a point is near the boundary
			bdryFlagY = -1,
			bdryFlag2X = -1,
			bdryFlag2Y = -1,
			bottom,
			top,
			left,
			right;
	double	uvX = 0.0,
			uvY = 0.0,
			Xa = 1.0,
			Ya = 1.0,
			Xb = 1.0,
			Yb = 1.0,
			eps = 1.e-10,
			x,
			y;
	double testx,
			testy;//flag
	distance_from_intersection_to_node[iu]		= xu[I];
	distance_from_u_to_body[ip] = 0;
	distance_from_u_to_body[ip] = 0;
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
		if (by[bottom]-eps < yu[J] && by[top]-eps > yu[J])
		{
			// if the segment is not parallel to the x-direction
			if (fabs(by[l]-by[k]) > eps)
			{
				// calculate the point of intersection of the double ray with the boundary
				x = bx[k] + (bx[l]-bx[k]) * (yu[J]-by[k])/(by[l]-by[k]);

				// calculate the body velocity at the point of intersection
				uvX = uB[k] + (uB[l]-uB[k]) * (yu[J]-by[k])/(by[l]-by[k]);

				// if the point of intersection coincides with the grid point
				if (fabs(x-xu[I]) < eps)
				{
					outsideX  = true;
					bdryFlagX = iu;
					tagsIn[iu] = iu;
					Xa        = 0.0;
					Xb        = 1.0;
					flag      = true; // flag is true when the point of intersection coincides with the grid point
					//printf("tagpoints warning\n\n\n\n\n");
				}
				// if the point of intersection lies to the right of the grid point (right-facing ray intersects the boundary)
				else if (x > xu[I]+eps)
					outsideX = !outsideX;

				// if the point of intersection is in the cell to the immediate left of the grid point
				if (x>xu[I-1]+eps && x<xu[I]-eps)
				{
					bdryFlagX  = iu;
					testx=x;//flag
					bdryFlag2X = iu+1;
					if (tags[iu-1]==-1)
						tagsIn[iu-1]	= iu-1;
					Xa = xu[I]-x;
					Xb = xu[I+1]-xu[I];
					if (x > midX)
						distance_from_u_to_body[ip] = Xa;
					//case 1
				}
				// if the point of intersection is in the cell to the immediate right of the grid point
				else if (x>xu[I]+eps && x<xu[I+1]-eps)
				{
					bdryFlagX  = iu;
					testx =x;//flag
					bdryFlag2X = iu-1;
					if(tags[iu+1] == -1)
						tagsIn[iu+1]	= iu+1;
					Xa = x-xu[I];
					Xb = xu[I]-xu[I-1];
					if (x < midX)
						distance_from_u_to_body[ip+1] = Xa;
					//case 2
				}
			}
		}
		// consider rays along the y-direction
		// if the ray intersects the boundary segment (right endpoint must be strictly to the right of ray; left can be on or to the left of the ray)
		if ( (bx[left]-eps < xu[I]) && (bx[right]-eps > xu[I]) && ( !flag ) ) // no need to do this part if flag is false
		{
			// if the segment is not parallel to the y-direction
			if (fabs(bx[l]-bx[k]) > eps)
			{
				// calculate the point of intersection of the double ray with the boundary
				y = by[k] + (by[l]-by[k]) * (xu[I]-bx[k]) / (bx[l]-bx[k]);

				// calculate the body velocity at the point of intersection
				uvY = uB[k] + (uB[l]-uB[k]) * (xu[I]-bx[k])/(bx[l]-bx[k]);

				// if the point of intersection coincides with the grid point
				if (fabs(y-yu[J]) < eps)
				{
					outsideY  = true; // then the point is considered to be outside the grid
					bdryFlagY = iu;    // the point is considered to be a forcing point, with index iu
					bdryFlag2Y= iu;
					tagsIn[iu]=iu;
					Ya        = 0.0;  // the coefficient for the linear interpolation during forcing
					Yb        = 1.0;
					flag      = true; // flag is true when the point of intersection coincides with the grid point
					//printf("tagpoitns warning\n\n\n\n\n");
				}
				// if the point of intersection lies to the top of the grid point
				else if (y > yu[J]+eps)
					outsideY = !outsideY; // then flip if inside or outside (start with true, i.e. outside) //this seems wrong too

				// if point of intersection is just below the concerned grid point
				if (y>yu[J-1]+eps && y<yu[J]-eps)
				{
					bdryFlagY = iu;
					testy = yu[J];//flag
					bdryFlag2Y= iu+(nx-1);
					//if (outsideY)
					if(tags[iu-nx+1]==-1)
						tagsIn[iu-(nx-1)]=iu-(nx-1);
					Ya = yu[J]-y;
					Yb = yu[J+1]-yu[J];
				}
				// if point of intersection is just above the concerned grid point
				else if (y>yu[J]+eps && y<yu[J+1]-eps)
				{
					bdryFlagY = iu;
					testy=1;//flag
					bdryFlag2Y= iu-(nx-1);
					//if (outsideY)
					if (tags[iu+nx-1]==-1)
						tagsIn[iu+(nx-1)]=iu+(nx-1);
					Ya = y-yu[J];
					Yb = yu[J]-yu[J-1];
				}
			}
		}
		k = l;
		l = l+1;
	}

	if (outsideX && bdryFlagX>=0)
	{
		tagsIn[iu]	= -1;
		tags[iu]	= bdryFlagX;
		tags2[iu]	= bdryFlag2X;
		distance_from_intersection_to_node[iu]		= Xa;
		//distance_from_intersection_to_node[iu]		= testx;
		distance_between_nodes_at_IB[iu]		= Xb;
		uv[iu]		= uvX;
	}
	else if (outsideY && bdryFlagY>=0)
	{
		tagsIn[iu]	= -1;
		tags[iu]	= bdryFlagY;
		tags2[iu]	= bdryFlag2Y;
		distance_from_intersection_to_node[iu]		= Ya;
		//distance_from_intersection_to_node[iu] = testy;
		distance_between_nodes_at_IB[iu]		= Yb;
		uv[iu]		= uvY;
	}
}

__global__
void tag_v(int *tags, int *tagsIn, int *tags2, double *bx, double *by, double *uB, double *vB, double *yv, double *xv,
		   double *distance_from_intersection_to_node, double *distance_between_nodes_at_IB, double *distance_from_u_to_body, double *distance_from_v_to_body, double *uv,
		   int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints, double midX, double midY)
{
	// calculate indicies indices
	int idx	= threadIdx.x + blockDim.x * blockIdx.x,
		i	= idx % (i_end-i_start),
		j	= idx / (i_end-i_start),
		I	= i_start + i,
		J	= j_start + j,
		iv	= J*nx + I + (nx-1)*ny,
		ip	= J*nx + I;

	// return if out of bounds of the array
	if (iv >= (nx-1)*ny + nx*(ny-1))
			return;

	// initial indices of the points on the body that define the segment under consideration
	int 	k = totalPoints-1,
			l = 0;

	// logic for the segment
	bool	outsideX = true,
			outsideY = true,
			flag = false;
	int		bdryFlagX = -1,  // stores if a point is near the boundary
			bdryFlagY = -1,
			bdryFlag2X = -1,
			bdryFlag2Y = -1,
			bottom,
			top,
			left,
			right;
	double	uvX = 0.0,
			uvY = 0.0,
			Xa = 1.0,
			Ya = 1.0,
			Xb = 1.0,
			Yb = 1.0,
			eps = 1.e-10,
			x,
			y;

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
		if (by[bottom]-eps < yv[J] && by[top]-eps > yv[J] && !flag)
		{
			// if the segment is not parallel to the ray
			if (fabs(by[l]-by[k]) > eps)
			{
				// calculate the point of intersection of the double ray with the boundary
				x = bx[k] + (bx[l]-bx[k]) * (yv[J]-by[k])/(by[l]-by[k]);
				// calculate the body velocity at the point of intersection
				uvX = vB[k] + (vB[l]-vB[k]) * (yv[J]-by[k])/(by[l]-by[k]);

				// if the point of intersection coincides with the grid point
				if (fabs(x-xv[I]) < eps)
				{
					outsideX   = true;
					bdryFlagX = iv;
					bdryFlag2X= iv;
					tagsIn[iv] = iv;
					Xa        = 0.0;
					Xb        = 1.0;
					flag      = true;
				}
				// if the point of intersection lies to the right of the grid point
				else if (x > xv[I]+eps)
					outsideX = !outsideX;

				// if the point of intersection is in the cell to the immediate right of the grid point
				if (x>xv[I-1]+eps && x<xv[I]-eps)
				{
					bdryFlagX  = iv;
					bdryFlag2X = iv+1;
					if (tags[iv-1]==-1)
						tagsIn[iv-1] = iv-1;
					Xa = xv[I]-x;
					Xb = xv[I+1]-xv[I];
				}
				// if the point of intersection is in the cell to the immediate left of the grid point
				else if (x>xv[I]+eps && x<xv[I+1]-eps)
				{
					bdryFlagX  = iv;
					bdryFlag2X = iv-1;
					if (tags[iv+1] == -1)
						tagsIn[iv+1] = iv+1;
					Xa = x-xv[I];
					Xb = xv[I]-xv[I-1];
				}
			}
		}
		// consider rays along the y-direction
		if (bx[left]-eps < xv[I] && bx[right]-eps > xv[I] && !flag)
		{
			// if the segment is not parallel to the ray
			if (fabs(bx[l]-bx[k]) > eps)
			{
				// calculate the point of intersection of the double ray with the boundary
				y = by[k] + (by[l]-by[k]) * (xv[I]-bx[k])/(bx[l]-bx[k]);
				// calculate the body velocity at the point of intersectioin
				uvY = vB[k] + (vB[l]-vB[k]) * (xv[I]-bx[k])/(bx[l]-bx[k]);

				// if the point of intersection coincides with the grid point
				if (fabs(y-yv[J]) < eps)
				{
					outsideY   = true;
					bdryFlagY = iv;
					bdryFlag2Y= iv;
					tagsIn[iv] = iv;
					Ya        = 0.0;
					Yb        = 1.0;
					flag      = true;
				}
				// if the point of intersection lies to the top of the grid point
				else if (y > yv[J]+eps)
					outsideY = !outsideY;

				if (y>yv[J-1]+eps && y<yv[J]-eps)
				{
					bdryFlagY  = iv;
					bdryFlag2Y = iv+nx;
					if (tags[iv-nx] == -1)
						tagsIn[iv-nx] = iv-nx;
					Ya = yv[J]-y;
					Yb = yv[J+1]-yv[J];
					//case 3
					if (outsideY)
						distance_from_v_to_body[ip] = Ya;
				}
				else if (y>yv[J]+eps && y<yv[J+1]-eps)
				{
					bdryFlagY  = iv;
					bdryFlag2Y = iv-nx;
					if(tags[iv+nx] == -1)
						tagsIn[iv+nx] = iv+nx;
					Ya = y-yv[J];
					Yb = yv[J]-yv[J-1];
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
		tagsIn[iv] = -1;
		tags[iv]    = bdryFlagY;
		tags2[iv]   = bdryFlag2Y;
		distance_from_intersection_to_node[iv]  = Ya;
		distance_between_nodes_at_IB[iv] = Yb;
		uv[iv]      = uvY;
	}
	else if (outsideX && bdryFlagX>=0)
	{
		tagsIn[iv] = -1;
		tags[iv]    = bdryFlagX;
		tags2[iv]   = bdryFlag2X;
		distance_from_intersection_to_node[iv]  = Xa;
		distance_between_nodes_at_IB[iv] = Xb;
		uv[iv]      = uvX;
	}
}

__global__
void tag_p(int *tagsP, int *tagsPOut, double *bx, double *by, double *yu, double *xv,
		   int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints, double midX, double midY)
{
	// calculate indicies indices
	int idx	= threadIdx.x + blockDim.x * blockIdx.x,
		i	= idx % (i_end-i_start),
		j	= idx / (i_end-i_start),
		I	= i_start + i,
		J	= j_start + j,
		ip	= J*nx + I;

	// return if out of bounds of the array
	if (ip >= nx*ny)
			return;

	// initial indices of the points on the body that define the segment under consideration
	int 	k = totalPoints-1,
			l = 0;

	// logic for the segment
	bool	flag = false;
	int		bottom,
			top,
			left,
			right;
	double	eps = 1.e-10,
			x,
			y;

	//tagsP[ip] = 1000;
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
		if (by[bottom]-eps < yu[J] && by[top]-eps > yu[J] && !flag)
		{
			// if the segment is not parallel to the ray
			if (fabs(by[l]-by[k]) > eps)
			{
				// calculate the point of intersection of the double ray with the boundary
				x = bx[k] + (bx[l]-bx[k]) * (yu[J]-by[k])/(by[l]-by[k]);

				// if the point of intersection coincides with the grid point
				/*if (fabs(x-xv[I]) < eps)
				{
					std::cout<<"tagpoints warning\n\n\n\n\n";
				}*/
				// just inside, right of mid
				if (x > midX + eps && x > xv[I] + eps && x < xv[I+1] + eps)
				{
					tagsP[ip] = ip;
				}
				// just inside, left of mid
				else if (x < midX + eps && x < xv[I] +eps && x > xv[I-1])
				{
					tagsP[ip] = ip;
				}

				// just inside, right of mid
				if (x > midX + eps && x < xv[I] + eps && x > xv[I-1] + eps)
				{
					tagsPOut[ip] = ip;
				}
				// just inside, left of mid
				else if (x < midX + eps && x > xv[I] +eps && x < xv[I+1] + eps)
				{
					tagsPOut[ip] = ip;
				}
			}
		}
		// consider rays along the y-direction
		if (bx[left]-eps < xv[I] && bx[right]-eps > xv[I] && !flag)
		{
			// if the segment is not parallel to the ray
			if (fabs(bx[l]-bx[k]) > eps)
			{
				// calculate the point of intersection of the double ray with the boundary
				y = by[k] + (by[l]-by[k]) * (xv[I]-bx[k])/(bx[l]-bx[k]);

				// if the point of intersection coincides with the grid point
				/*if (fabs(y-yu[J]) < eps)
				{
					std::cout<<"tagpoints warning\n\n\n\n\n";
				}*/
				// just inside, north of mid
				if (y > midY + eps && y > yu[J] +eps && y < yu[J+1])
				{
					tagsP[ip] = ip;
				}
				//just inside, south of mid
				else if (y < midY + eps && y < yu[J] +eps && y > yu[J-1] +eps)
				{
					tagsP[ip] = ip;
				}

				// just inside, north of mid
				if (y > midY + eps && y < yu[J] +eps && y > yu[J-1])
				{
					tagsPOut[ip] = ip;
				}
				//just inside, south of mid
				else if (y < midY + eps && y > yu[J] +eps && y < yu[J+1] +eps)
				{
					tagsPOut[ip] = ip;
				}
			}
		}
		k = l;
		l = l+1;
	}// end while
}

__global__
void zero_pressure(int *tagsP,  int i_start, int j_start, int i_end, int j_end, int nx, int ny)
{
	// calculate indicies indices
	int j	= threadIdx.x + blockDim.x * blockIdx.x,
		J	= j_start + j,
		I;

	if (J > j_end)
		return;

	bool 	rowIsntDone = true,
			flag 		= false;

	for (int i=i_start; i<i_end; i++)
	{
		I = J*nx+i;
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
					if (tagsP[J*nx+k] != -1 )
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

__global__
void zero_x(int *tagsIn,  int i_start, int j_start, int i_end, int j_end, int nx, int ny)
{
	// calculate indicies indices
	int j	= threadIdx.x + blockDim.x * blockIdx.x,
		J	= j_start + j,
		I;

	if (J > j_end)
		return;

	bool 	rowIsntDone = true,
			flag 		= false;
	for (int i=i_start; i<i_end; i++)
	{
		I = J*(nx-1)+i;

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
					if (tagsIn[J*(nx-1)+k] != -1 )
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

__global__
void zero_y(int *tagsIn,  int i_start, int j_start, int i_end, int j_end, int nx, int ny)
{
	// calculate indicies indices
	int j	= threadIdx.x + blockDim.x * blockIdx.x,
		J	= j_start + j,
		I;

	if (J > j_end)
		return;

	bool 	rowIsntDone = true,
			flag 		= false;

	for (int i=i_start; i<i_end; i++)
	{
		I = J*(nx)+i + (nx-1)*ny;
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
					if (tagsIn[J*(nx)+k + (nx-1)*ny] != -1 )
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




}
