#include <solvers/NavierStokes/oscCylinder.h>

void oscCylinder::cast()
{
	cfl.resize(nx*ny);
	distance.resize((nx-1)*ny + (ny-1)*nx);

	cfl_r		= thrust::raw_pointer_cast( &(cfl[0]) );
	distance_r	= thrust::raw_pointer_cast( &(distance[0]) );
}
