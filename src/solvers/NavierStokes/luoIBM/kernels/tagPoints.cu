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
void tag_u_luo(int *hybridTagsUV, int *ghostTagsUV, int *hybridTagsUV2, double *bx, double *by, double *uB, double *vB, double *yu, double *xu,
			double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
			double *x1, double *y1, double *x2, double *y2, //testing
			double *distance_from_intersection_to_node, double *distance_between_nodes_at_IB, double *distance_from_u_to_body, double *distance_from_v_to_body, double *uv,
			int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints, double midX, double midY)
{
	// calculate indices
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
			y,
			p,o,b,a,
			theta_321;

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
					ghostTagsUV[iu] = iu;
					Xa        = 0.0;
					Xb        = 1.0;
					flag      = true; // flag is true when the point of intersection coincides with the grid point
					//printf("tagpoints warning\n\n\n\n\n");
				}
				// if the point of intersection lies to the right of the grid point (right-facing ray intersects the boundary)
				else if (x > xu[I]+eps)
					outsideX = !outsideX;

				// if the point of intersection is in the cell to the immediate left of the grid point
				//right of body
				if (x>xu[I-1]+eps && x<xu[I]-eps)
				{
					bdryFlagX  = iu;
					bdryFlag2X = iu+1;
					if (x>midX+eps)
					{
						ghostTagsUV[iu-1]	= iu-1;
						//calculate image point and body point for ghost node

						x1[iu-1] = bx[l];
						y1[iu-1] = by[l];
						x2[iu-1] = bx[k];
						y2[iu-1] = by[k];

						p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
						o = sqrt(pow((bx[l]-xu[I-1]),2) + pow((by[l]-yu[J]),2)); //distance from 3 to 1
						b = sqrt(pow((xu[I-1]-bx[k]),2) + pow((yu[J]-by[k]),2)); //distance from 2 to 3

						theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

						a = sin(theta_321)*b;//distance from 3 to 4

						body_intercept_y[iu-1] = by[k] + (by[l]-by[k])/p * a/tan(theta_321);
						body_intercept_x[iu-1] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
						image_point_y[iu-1] = body_intercept_y[iu-1]+(body_intercept_y[iu-1]-yu[J]);//flag change so it doesn't access memory twice
						image_point_x[iu-1] = body_intercept_x[iu-1]+(body_intercept_x[iu-1]-xu[I-1]);
						//calculate ip and bp for hybrid node
						x1[iu] = bx[l];
						y1[iu] = by[l];
						x2[iu] = bx[k];
						y2[iu] = by[k];

						p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
						o = sqrt(pow((bx[l]-xu[I]),2) + pow((by[l]-yu[J]),2)); //distance from 3 to 1
						b = sqrt(pow((xu[I]-bx[k]),2) + pow((yu[J]-by[k]),2)); //distance from 2 to 3

						theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

						a = sin(theta_321)*b;//distance from 3 to 4

						body_intercept_y[iu] = by[k] + (by[l]-by[k])/p * a/tan(theta_321);
						body_intercept_x[iu] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
						image_point_y[iu] = 2*yu[J] - body_intercept_y[iu];
						image_point_x[iu] = 2*xu[I] - body_intercept_x[iu];
					}
					Xa = xu[I]-x;
					Xb = xu[I+1]-xu[I];
					if (x > midX)
						distance_from_u_to_body[ip] = Xa;
					//case 1
				}
				// if the point of intersection is in the cell to the immediate right of the grid point
				//left of body
				else if (x>xu[I]+eps && x<xu[I+1]-eps)
				{
					bdryFlagX  = iu;
					bdryFlag2X = iu-1;
					if(x<midX-eps)
					{
						ghostTagsUV[iu+1]	= iu+1;
						x1[iu+1] = bx[l];
						y1[iu+1] = by[l];
						x2[iu+1] = bx[k];
						y2[iu+1] = by[k];

						p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
						o = sqrt(pow((bx[l]-xu[I+1]),2) + pow((by[l]-yu[J]),2)); //distance from 3 to 1
						b = sqrt(pow((xu[I+1]-bx[k]),2) + pow((yu[J]-by[k]),2)); //distance from 2 to 3

						theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

						a = sin(theta_321)*b;//distance from 3 to 4

						body_intercept_y[iu+1] = by[k] + (by[l]-by[k])/p * a/tan(theta_321);
						body_intercept_x[iu+1] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
						image_point_y[iu+1] = body_intercept_y[iu+1]+(body_intercept_y[iu+1]-yu[J]);
						image_point_x[iu+1] = body_intercept_x[iu+1]+(body_intercept_x[iu+1]-xu[I+1]);
						//calculate ip and bp for hybrid node
						x1[iu] = bx[l];
						y1[iu] = by[l];
						x2[iu] = bx[k];
						y2[iu] = by[k];

						p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
						o = sqrt(pow((bx[l]-xu[I]),2) + pow((by[l]-yu[J]),2)); //distance from 3 to 1
						b = sqrt(pow((xu[I]-bx[k]),2) + pow((yu[J]-by[k]),2)); //distance from 2 to 3

						theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

						a = sin(theta_321)*b;//distance from 3 to 4

						body_intercept_y[iu] = by[k] + (by[l]-by[k])/p * a/tan(theta_321);
						body_intercept_x[iu] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
						image_point_y[iu] = 2*yu[J] - body_intercept_y[iu];
						image_point_x[iu] = 2*xu[I] - body_intercept_x[iu];

					}
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
					ghostTagsUV[iu]=iu;
					Ya        = 0.0;  // the coefficient for the linear interpolation during forcing
					Yb        = 1.0;
					flag      = true; // flag is true when the point of intersection coincides with the grid point
					//printf("tagpoitns warning\n\n\n\n\n");
				}
				// if the point of intersection lies to the top of the grid point
				else if (y > yu[J]+eps)
					outsideY = !outsideY; // then flip if inside or outside (start with true, i.e. outside)

				// if point of intersection is just below the concerned grid point
				//above body
				if (y>yu[J-1]+eps && y<yu[J]-eps)
				{
					bdryFlagY = iu;
					bdryFlag2Y= iu+(nx-1);
					//if (outsideY)
					if(y>midY+eps)
					{
						ghostTagsUV[iu-(nx-1)]=iu-(nx-1);
						x1[iu-(nx-1)] = bx[l];
						y1[iu-(nx-1)] = by[l];
						x2[iu-(nx-1)] = bx[k];
						y2[iu-(nx-1)] = by[k];

						p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
						o = sqrt(pow((bx[l]-xu[I]),2) + pow((by[l]-yu[J-1]),2)); //distance from 3 to 1
						b = sqrt(pow((xu[I]-bx[k]),2) + pow((yu[J-1]-by[k]),2)); //distance from 2 to 3

						theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

						a = sin(theta_321)*b;//distance from 3 to 4

						body_intercept_y[iu-(nx-1)] = by[k] + (by[l]-by[k])/p * a/tan(theta_321);
						body_intercept_x[iu-(nx-1)] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
						image_point_y[iu-(nx-1)] = body_intercept_y[iu-(nx-1)]+(body_intercept_y[iu-(nx-1)]-yu[J-1]);
						image_point_x[iu-(nx-1)] = body_intercept_x[iu-(nx-1)]+(body_intercept_x[iu-(nx-1)]-xu[I]);
						//calculate ip and bp for hybrid node
						x1[iu] = bx[l];
						y1[iu] = by[l];
						x2[iu] = bx[k];
						y2[iu] = by[k];

						p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
						o = sqrt(pow((bx[l]-xu[I]),2) + pow((by[l]-yu[J]),2)); //distance from 3 to 1
						b = sqrt(pow((xu[I]-bx[k]),2) + pow((yu[J]-by[k]),2)); //distance from 2 to 3

						theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

						a = sin(theta_321)*b;//distance from 3 to 4

						body_intercept_y[iu] = by[k] + (by[l]-by[k])/p * a/tan(theta_321);
						body_intercept_x[iu] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
						image_point_y[iu] = 2*yu[J] - body_intercept_y[iu];
						image_point_x[iu] = 2*xu[I] - body_intercept_x[iu];
					}
					Ya = yu[J]-y;
					Yb = yu[J+1]-yu[J];
				}
				// if point of intersection is just above the concerned grid point
				//below body
				else if (y>yu[J]+eps && y<yu[J+1]-eps)
				{
					bdryFlagY = iu;
					bdryFlag2Y= iu-(nx-1);
					//if (outsideY)
					if (y<midY-eps)
					{
						ghostTagsUV[iu+(nx-1)]=iu+(nx-1);
						x1[iu+(nx-1)] = bx[l];
						y1[iu+(nx-1)] = by[l];
						x2[iu+(nx-1)] = bx[k];
						y2[iu+(nx-1)] = by[k];

						p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
						o = sqrt(pow((bx[l]-xu[I]),2) + pow((by[l]-yu[J+1]),2)); //distance from 3 to 1
						b = sqrt(pow((xu[I]-bx[k]),2) + pow((yu[J+1]-by[k]),2)); //distance from 2 to 3

						theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

						a = sin(theta_321)*b;//distance from 3 to 4

						body_intercept_y[iu+(nx-1)] = by[k] + (by[l]-by[k])/p * a/tan(theta_321);
						body_intercept_x[iu+(nx-1)] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
						image_point_y[iu+(nx-1)] = body_intercept_y[iu+(nx-1)]+(body_intercept_y[iu+(nx-1)]-yu[J+1]);
						image_point_x[iu+(nx-1)] = body_intercept_x[iu+(nx-1)]+(body_intercept_x[iu+(nx-1)]-xu[I]);
						//calculate ip and bp for hybrid node
						x1[iu] = bx[l];
						y1[iu] = by[l];
						x2[iu] = bx[k];
						y2[iu] = by[k];

						p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
						o = sqrt(pow((bx[l]-xu[I]),2) + pow((by[l]-yu[J]),2)); //distance from 3 to 1
						b = sqrt(pow((xu[I]-bx[k]),2) + pow((yu[J]-by[k]),2)); //distance from 2 to 3

						theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

						a = sin(theta_321)*b;//distance from 3 to 4

						body_intercept_y[iu] = by[k] + (by[l]-by[k])/p * a/tan(theta_321);
						body_intercept_x[iu] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
						image_point_y[iu] = 2*yu[J] - body_intercept_y[iu];
						image_point_x[iu] = 2*xu[I] - body_intercept_x[iu];
					}
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
		ghostTagsUV[iu]	= -1;
		hybridTagsUV[iu]	= bdryFlagX;
		hybridTagsUV2[iu]	= bdryFlag2X;
		distance_from_intersection_to_node[iu]		= Xa;
		distance_between_nodes_at_IB[iu]		= Xb;
		uv[iu]		= uvX;
	}
	else if (outsideY && bdryFlagY>=0)
	{
		ghostTagsUV[iu]	= -1;
		hybridTagsUV[iu]	= bdryFlagY;
		hybridTagsUV2[iu]	= bdryFlag2Y;
		distance_from_intersection_to_node[iu]		= Ya;
		distance_between_nodes_at_IB[iu]		= Yb;
		uv[iu]		= uvY;
	}
}

__global__
void tag_v_luo(int *hybridTagsUV, int *ghostTagsUV, int *hybridTagsUV2, double *bx, double *by, double *uB, double *vB, double *yv, double *xv,
				double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y, double *x1, double *y1, double *x2, double *y2,
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
			y,
			p,o,b,a,
			theta_321;
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
					ghostTagsUV[iv] = iv;
					Xa        = 0.0;
					Xb        = 1.0;
					flag      = true;
				}
				// if the point of intersection lies to the right of the grid point
				else if (x > xv[I]+eps)
					outsideX = !outsideX;

				//right of body
				if (x>xv[I-1]+eps && x<xv[I]-eps)
				{
					bdryFlagX  = iv;
					bdryFlag2X = iv+1;
					if (x>midX+eps)
					{
						ghostTagsUV[iv-1] = iv-1;
						//ghost node interpolation
						x1[iv-1] = bx[l];
						y1[iv-1] = by[l];
						x2[iv-1] = bx[k];
						y2[iv-1] = by[k];

						p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
						o = sqrt(pow((bx[l]-xv[I-1]),2) + pow((by[l]-yv[J]),2)); //distance from 3 to 1
						b = sqrt(pow((xv[I-1]-bx[k]),2) + pow((yv[J]-by[k]),2)); //distance from 2 to 3

						theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

						a = sin(theta_321)*b;//distance from 3 to 4

						body_intercept_y[iv-1] = by[k] + (by[l]-by[k])/p * a/tan(theta_321); //body intercept point, forms a right angle with lines 1-2 and 4-3
						body_intercept_x[iv-1] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
						image_point_y[iv-1] = body_intercept_y[iv-1]+(body_intercept_y[iv-1]-yv[J]); //image point, mirror the ghost node accross line 1-2
						image_point_x[iv-1] = body_intercept_x[iv-1]+(body_intercept_x[iv-1]-xv[I-1]);
						//hybrid node interpolation
						x1[iv] = bx[l];
						y1[iv] = by[l];
						x2[iv] = bx[k];
						y2[iv] = by[k];

						p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
						o = sqrt(pow((bx[l]-xv[I]),2) + pow((by[l]-yv[J]),2)); //distance from 3 to 1
						b = sqrt(pow((xv[I]-bx[k]),2) + pow((yv[J]-by[k]),2)); //distance from 2 to 3

						theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

						a = sin(theta_321)*b;//distance from 3 to 4

						body_intercept_y[iv] = by[k] + (by[l]-by[k])/p * a/tan(theta_321); //body intercept point, forms a right angle with lines 1-2 and 4-3
						body_intercept_x[iv] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
						image_point_y[iv] = 2*yv[J] - body_intercept_y[iv];
						image_point_x[iv] = 2*xv[I] - body_intercept_x[iv];
					}
					Xa = xv[I]-x;
					Xb = xv[I+1]-xv[I];
				}
				// left of body
				else if (x>xv[I]+eps && x<xv[I+1]-eps)
				{
					bdryFlagX  = iv;
					bdryFlag2X = iv-1;
					if (x<midX-eps)
					{
						ghostTagsUV[iv+1] = iv+1;
						x1[iv+1] = bx[l];
						y1[iv+1] = by[l];
						x2[iv+1] = bx[k];
						y2[iv+1] = by[k];

						p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
						o = sqrt(pow((bx[l]-xv[I+1]),2) + pow((by[l]-yv[J]),2)); //distance from 3 to 1
						b = sqrt(pow((xv[I+1]-bx[k]),2) + pow((yv[J]-by[k]),2)); //distance from 2 to 3

						theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

						a = sin(theta_321)*b;//distance from 3 to 4

						body_intercept_y[iv+1] = by[k] + (by[l]-by[k])/p * a/tan(theta_321);
						body_intercept_x[iv+1] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
						image_point_y[iv+1] = body_intercept_y[iv+1]+(body_intercept_y[iv+1]-yv[J]);
						image_point_x[iv+1] = body_intercept_x[iv+1]+(body_intercept_x[iv+1]-xv[I+1]);
						//hybrid node interpolation
						x1[iv] = bx[l];
						y1[iv] = by[l];
						x2[iv] = bx[k];
						y2[iv] = by[k];

						p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
						o = sqrt(pow((bx[l]-xv[I]),2) + pow((by[l]-yv[J]),2)); //distance from 3 to 1
						b = sqrt(pow((xv[I]-bx[k]),2) + pow((yv[J]-by[k]),2)); //distance from 2 to 3

						theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

						a = sin(theta_321)*b;//distance from 3 to 4

						body_intercept_y[iv] = by[k] + (by[l]-by[k])/p * a/tan(theta_321); //body intercept point, forms a right angle with lines 1-2 and 4-3
						body_intercept_x[iv] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
						image_point_y[iv] = 2*yv[J] - body_intercept_y[iv];
						image_point_x[iv] = 2*xv[I] - body_intercept_x[iv];
					}
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
					ghostTagsUV[iv] = iv;
					Ya        = 0.0;
					Yb        = 1.0;
					flag      = true;
				}
				// if the point of intersection lies to the top of the grid point
				else if (y > yv[J]+eps)
					outsideY = !outsideY;

				//above body
				if (y>yv[J-1]+eps && y<yv[J]-eps)
				{
					bdryFlagY  = iv;
					bdryFlag2Y = iv+nx;
					if (y>midY+eps)
					{
						ghostTagsUV[iv-nx] = iv-nx;
						x1[iv-nx] = bx[l];
						y1[iv-nx] = by[l];
						x2[iv-nx] = bx[k];
						y2[iv-nx] = by[k];

						p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
						o = sqrt(pow((bx[l]-xv[I]),2) + pow((by[l]-yv[J-1]),2)); //distance from 3 to 1
						b = sqrt(pow((xv[I]-bx[k]),2) + pow((yv[J-1]-by[k]),2)); //distance from 2 to 3

						theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

						a = sin(theta_321)*b;//distance from 3 to 4

						body_intercept_y[iv-nx] = by[k] + (by[l]-by[k])/p * a/tan(theta_321);
						body_intercept_x[iv-nx] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
						image_point_y[iv-nx] = body_intercept_y[iv-nx]+(body_intercept_y[iv-nx]-yv[J-1]);
						image_point_x[iv-nx] = body_intercept_x[iv-nx]+(body_intercept_x[iv-nx]-xv[I]);
						//hybrid node interpolation
						x1[iv] = bx[l];
						y1[iv] = by[l];
						x2[iv] = bx[k];
						y2[iv] = by[k];

						p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
						o = sqrt(pow((bx[l]-xv[I]),2) + pow((by[l]-yv[J]),2)); //distance from 3 to 1
						b = sqrt(pow((xv[I]-bx[k]),2) + pow((yv[J]-by[k]),2)); //distance from 2 to 3

						theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

						a = sin(theta_321)*b;//distance from 3 to 4

						body_intercept_y[iv] = by[k] + (by[l]-by[k])/p * a/tan(theta_321); //body intercept point, forms a right angle with lines 1-2 and 4-3
						body_intercept_x[iv] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
						image_point_y[iv] = 2*yv[J] - body_intercept_y[iv];
						image_point_x[iv] = 2*xv[I] - body_intercept_x[iv];
					}
					Ya = yv[J]-y;
					Yb = yv[J+1]-yv[J];
					//case 3
					if (outsideY)
						distance_from_v_to_body[ip] = Ya;
				}
				// if point of intersection is just below the concerned grid point
				//below body
				else if (y>yv[J]+eps && y<yv[J+1]-eps)
				{
					bdryFlagY  = iv;
					bdryFlag2Y = iv-nx;
					if(y<midY-eps)
					{
						ghostTagsUV[iv+nx] = iv+nx;
						x1[iv+nx] = bx[l];
						y1[iv+nx] = by[l];
						x2[iv+nx] = bx[k];
						y2[iv+nx] = by[k];

						p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
						o = sqrt(pow((bx[l]-xv[I]),2) + pow((by[l]-yv[J+1]),2)); //distance from 3 to 1
						b = sqrt(pow((xv[I]-bx[k]),2) + pow((yv[J+1]-by[k]),2)); //distance from 2 to 3

						theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

						a = sin(theta_321)*b;//distance from 3 to 4

						body_intercept_y[iv+nx] = by[k] + (by[l]-by[k])/p * a/tan(theta_321);
						body_intercept_x[iv+nx] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
						image_point_y[iv+nx] = body_intercept_y[iv+nx]+(body_intercept_y[iv+nx]-yv[J+1]);
						image_point_x[iv+nx] = body_intercept_x[iv+nx]+(body_intercept_x[iv+nx]-xv[I]);
						//hybrid node interpolation
						x1[iv] = bx[l];
						y1[iv] = by[l];
						x2[iv] = bx[k];
						y2[iv] = by[k];

						p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
						o = sqrt(pow((bx[l]-xv[I]),2) + pow((by[l]-yv[J]),2)); //distance from 3 to 1
						b = sqrt(pow((xv[I]-bx[k]),2) + pow((yv[J]-by[k]),2)); //distance from 2 to 3

						theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

						a = sin(theta_321)*b;//distance from 3 to 4

						body_intercept_y[iv] = by[k] + (by[l]-by[k])/p * a/tan(theta_321); //body intercept point, forms a right angle with lines 1-2 and 4-3
						body_intercept_x[iv] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
						image_point_y[iv] = 2*yv[J] - body_intercept_y[iv];
						image_point_x[iv] = 2*xv[I] - body_intercept_x[iv];
					}
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
		ghostTagsUV[iv] = -1;
		hybridTagsUV[iv]    = bdryFlagY;
		hybridTagsUV2[iv]   = bdryFlag2Y;
		distance_from_intersection_to_node[iv]  = Ya;
		distance_between_nodes_at_IB[iv] = Yb;
		uv[iv]      = uvY;
	}
	else if (outsideX && bdryFlagX>=0)
	{
		ghostTagsUV[iv] = -1;
		hybridTagsUV[iv]    = bdryFlagX;
		hybridTagsUV2[iv]   = bdryFlag2X;
		distance_from_intersection_to_node[iv]  = Xa;
		distance_between_nodes_at_IB[iv] = Xb;
		uv[iv]      = uvX;
	}
}

__global__
void tag_p_luo(int *ghostTagsP, int *hybridTagsP, double *bx, double *by, double *yu, double *xv,
				double *body_intercept_p_x, double *body_intercept_p_y, double *image_point_p_x, double *image_point_p_y, double *x1_p, double *y1_p, double *x2_p, double *y2_p,
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
			y,
			p,o,b,a,
			theta_321;

	//ghostTagsP[ip] = 1000;
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

				// just inside, right of mid
				if (x > midX + eps && x > xv[I] + eps && x < xv[I+1] + eps)
				{
					ghostTagsP[ip] = ip;
					x1_p[ip] = bx[l];
					y1_p[ip] = by[l];
					x2_p[ip] = bx[k];
					y2_p[ip] = by[k];

					p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
					o = sqrt(pow((bx[l]-xv[I]),2) + pow((by[l]-yu[J]),2)); //distance from 3 to 1
					b = sqrt(pow((xv[I]-bx[k]),2) + pow((yu[J]-by[k]),2)); //distance from 2 to 3

					theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

					a = sin(theta_321)*b;//distance from 3 to 4

					body_intercept_p_y[ip] = by[k] + (by[l]-by[k])/p * a/tan(theta_321); //body intercept point, forms a right angle with lines 1-2 and 4-3
					body_intercept_p_x[ip] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
					image_point_p_y[ip] = body_intercept_p_y[ip]+(body_intercept_p_y[ip]-yu[J]); //image point, mirror the ghost node accross line 1-2
					image_point_p_x[ip] = body_intercept_p_x[ip]+(body_intercept_p_x[ip]-xv[I]);
				}
				// just inside, left of mid
				else if (x < midX + eps && x < xv[I] +eps && x > xv[I-1])
				{
					ghostTagsP[ip] = ip;
					x1_p[ip] = bx[l];
					y1_p[ip] = by[l];
					x2_p[ip] = bx[k];
					y2_p[ip] = by[k];

					p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
					o = sqrt(pow((bx[l]-xv[I]),2) + pow((by[l]-yu[J]),2)); //distance from 3 to 1
					b = sqrt(pow((xv[I]-bx[k]),2) + pow((yu[J]-by[k]),2)); //distance from 2 to 3

					theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

					a = sin(theta_321)*b;//distance from 3 to 4

					body_intercept_p_y[ip] = by[k] + (by[l]-by[k])/p * a/tan(theta_321); //body intercept point, forms a right angle with lines 1-2 and 4-3
					body_intercept_p_x[ip] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
					image_point_p_y[ip] = body_intercept_p_y[ip]+(body_intercept_p_y[ip]-yu[J]); //image point, mirror the ghost node accross line 1-2
					image_point_p_x[ip] = body_intercept_p_x[ip]+(body_intercept_p_x[ip]-xv[I]);
				}

				// just inside, right of mid
				if (x > midX + eps && x < xv[I] + eps && x > xv[I-1] + eps)
				{
					hybridTagsP[ip] = ip;
					//hybrid node interpolation
					x1_p[ip] = bx[l];
					y1_p[ip] = by[l];
					x2_p[ip] = bx[k];
					y2_p[ip] = by[k];

					p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
					o = sqrt(pow((bx[l]-xv[I]),2) + pow((by[l]-yu[J]),2)); //distance from 3 to 1
					b = sqrt(pow((xv[I]-bx[k]),2) + pow((yu[J]-by[k]),2)); //distance from 2 to 3

					theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

					a = sin(theta_321)*b;//distance from 3 to 4

					body_intercept_p_y[ip] = by[k] + (by[l]-by[k])/p * a/tan(theta_321); //body intercept point, forms a right angle with lines 1-2 and 4-3
					body_intercept_p_x[ip] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
					image_point_p_y[ip] = 2*yu[J] - body_intercept_p_y[ip];
					image_point_p_x[ip] = 2*xv[I] - body_intercept_p_x[ip];
				}
				// just inside, left of mid
				else if (x < midX + eps && x > xv[I] +eps && x < xv[I+1] + eps)
				{
					hybridTagsP[ip] = ip;
					//hybrid node interpolation
					x1_p[ip] = bx[l];
					y1_p[ip] = by[l];
					x2_p[ip] = bx[k];
					y2_p[ip] = by[k];

					p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
					o = sqrt(pow((bx[l]-xv[I]),2) + pow((by[l]-yu[J]),2)); //distance from 3 to 1
					b = sqrt(pow((xv[I]-bx[k]),2) + pow((yu[J]-by[k]),2)); //distance from 2 to 3

					theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

					a = sin(theta_321)*b;//distance from 3 to 4

					body_intercept_p_y[ip] = by[k] + (by[l]-by[k])/p * a/tan(theta_321); //body intercept point, forms a right angle with lines 1-2 and 4-3
					body_intercept_p_x[ip] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
					image_point_p_y[ip] = 2*yu[J] - body_intercept_p_y[ip];
					image_point_p_x[ip] = 2*xv[I] - body_intercept_p_x[ip];
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
					ghostTagsP[ip] = ip;
					x1_p[ip] = bx[l];
					y1_p[ip] = by[l];
					x2_p[ip] = bx[k];
					y2_p[ip] = by[k];

					p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
					o = sqrt(pow((bx[l]-xv[I]),2) + pow((by[l]-yu[J]),2)); //distance from 3 to 1
					b = sqrt(pow((xv[I]-bx[k]),2) + pow((yu[J]-by[k]),2)); //distance from 2 to 3

					theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

					a = sin(theta_321)*b;//distance from 3 to 4

					body_intercept_p_y[ip] = by[k] + (by[l]-by[k])/p * a/tan(theta_321); //body intercept point, forms a right angle with lines 1-2 and 4-3
					body_intercept_p_x[ip] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
					image_point_p_y[ip] = body_intercept_p_y[ip]+(body_intercept_p_y[ip]-yu[J]); //image point, mirror the ghost node accross line 1-2
					image_point_p_x[ip] = body_intercept_p_x[ip]+(body_intercept_p_x[ip]-xv[I]);
				}
				//just inside, south of mid
				else if (y < midY + eps && y < yu[J] +eps && y > yu[J-1] +eps)
				{
					ghostTagsP[ip] = ip;
					x1_p[ip] = bx[l];
					y1_p[ip] = by[l];
					x2_p[ip] = bx[k];
					y2_p[ip] = by[k];

					p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
					o = sqrt(pow((bx[l]-xv[I]),2) + pow((by[l]-yu[J]),2)); //distance from 3 to 1
					b = sqrt(pow((xv[I]-bx[k]),2) + pow((yu[J]-by[k]),2)); //distance from 2 to 3

					theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

					a = sin(theta_321)*b;//distance from 3 to 4

					body_intercept_p_y[ip] = by[k] + (by[l]-by[k])/p * a/tan(theta_321); //body intercept point, forms a right angle with lines 1-2 and 4-3
					body_intercept_p_x[ip] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
					image_point_p_y[ip] = body_intercept_p_y[ip]+(body_intercept_p_y[ip]-yu[J]); //image point, mirror the ghost node accross line 1-2
					image_point_p_x[ip] = body_intercept_p_x[ip]+(body_intercept_p_x[ip]-xv[I]);
				}

				// just inside, north of mid
				if (y > midY + eps && y < yu[J] +eps && y > yu[J-1])
				{
					hybridTagsP[ip] = ip;
					//hybrid node interpolation
					x1_p[ip] = bx[l];
					y1_p[ip] = by[l];
					x2_p[ip] = bx[k];
					y2_p[ip] = by[k];

					p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
					o = sqrt(pow((bx[l]-xv[I]),2) + pow((by[l]-yu[J]),2)); //distance from 3 to 1
					b = sqrt(pow((xv[I]-bx[k]),2) + pow((yu[J]-by[k]),2)); //distance from 2 to 3

					theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

					a = sin(theta_321)*b;//distance from 3 to 4

					body_intercept_p_y[ip] = by[k] + (by[l]-by[k])/p * a/tan(theta_321); //body intercept point, forms a right angle with lines 1-2 and 4-3
					body_intercept_p_x[ip] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
					image_point_p_y[ip] = 2*yu[J] - body_intercept_p_y[ip];
					image_point_p_x[ip] = 2*xv[I] - body_intercept_p_x[ip];
				}
				//just inside, south of mid
				else if (y < midY + eps && y > yu[J] +eps && y < yu[J+1] +eps)
				{
					hybridTagsP[ip] = ip;
					//hybrid node interpolation
					x1_p[ip] = bx[l];
					y1_p[ip] = by[l];
					x2_p[ip] = bx[k];
					y2_p[ip] = by[k];

					p = sqrt(pow((bx[l]-bx[k]),2) + pow((by[l]-by[k]),2)); //distance from 1 to 2
					o = sqrt(pow((bx[l]-xv[I]),2) + pow((by[l]-yu[J]),2)); //distance from 3 to 1
					b = sqrt(pow((xv[I]-bx[k]),2) + pow((yu[J]-by[k]),2)); //distance from 2 to 3

					theta_321 = acos((p*p+b*b-o*o)/(2*p*b)); //angle format: start point,vertex,endpoint

					a = sin(theta_321)*b;//distance from 3 to 4

					body_intercept_p_y[ip] = by[k] + (by[l]-by[k])/p * a/tan(theta_321); //body intercept point, forms a right angle with lines 1-2 and 4-3
					body_intercept_p_x[ip] = bx[k] + (bx[l]-bx[k])/p * a/tan(theta_321);
					image_point_p_y[ip] = 2*yu[J] - body_intercept_p_y[ip];
					image_point_p_x[ip] = 2*xv[I] - body_intercept_p_x[ip];
				}
			}
		}
		k = l;
		l = l+1;
	}// end while
}

__global__
void zero_pressure_luo(int *ghostTagsP,  int i_start, int j_start, int i_end, int j_end, int nx, int ny)
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
		if (ghostTagsP[I-1] >= 0)
		{
			if (ghostTagsP[I]!=-1 && ghostTagsP[I-1]==0)
			{
				rowIsntDone = false;
			}
			if(ghostTagsP[I]==-1 && rowIsntDone)
			{
				ghostTagsP[I]=0;
			}
			if(ghostTagsP[I]==0 && ghostTagsP[I-1] !=0)
			{
				int k = i;
				while (k < nx)
				{
					k++;
					if (ghostTagsP[J*nx+k] != -1 )
					{
						flag = true;
					}
				}
				if (!flag)
				{
					ghostTagsP[I]=-1;
					rowIsntDone = false;
				}
			}
		}
	}
}

__global__
void zero_x_luo(int *ghostTagsUV,  int i_start, int j_start, int i_end, int j_end, int nx, int ny)
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

		if (ghostTagsUV[I-1] >= 0)
		{
			if (ghostTagsUV[I]!=-1 && ghostTagsUV[I-1]==0)
			{
				rowIsntDone = false;
			}
			if(ghostTagsUV[I]==-1 && rowIsntDone)
			{
				ghostTagsUV[I]=0;
			}
			if(ghostTagsUV[I]==0 && ghostTagsUV[I-1] !=0)
			{
				int k = i;
				while (k < nx-1)
				{
					k++;
					if (ghostTagsUV[J*(nx-1)+k] != -1 )
					{
						flag = true;
					}
				}
				if (!flag)
				{
					ghostTagsUV[I]=-1;
					rowIsntDone = false;
				}
			}
		}
	}
}

__global__
void zero_y_luo(int *ghostTagsUV,  int i_start, int j_start, int i_end, int j_end, int nx, int ny)
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
		if (ghostTagsUV[I-1] >= 0)
		{
			if (ghostTagsUV[I]!=-1 && ghostTagsUV[I-1]==0)
			{
				rowIsntDone = false;
			}
			if(ghostTagsUV[I]==-1 && rowIsntDone)
			{
				ghostTagsUV[I]=0;
			}
			if(ghostTagsUV[I]==0 && ghostTagsUV[I-1] !=0)
			{
				int k = i;
				while (k < nx-1)
				{
					k++;
					if (ghostTagsUV[J*(nx)+k + (nx-1)*ny] != -1 )
					{
						flag = true;
					}
				}
				if (!flag)
				{
					ghostTagsUV[I]=-1;
					rowIsntDone = false;
				}
			}
		}
	}
}




}
