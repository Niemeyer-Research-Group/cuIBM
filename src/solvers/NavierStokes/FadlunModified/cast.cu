#include <solvers/NavierStokes/fadlunModified.h>

void fadlunModified::cast()
{
	//resize numUV
	tags.resize(numUV);
	tagsOld.resize(numUV);
	tagsPOld.resize(numP);
	tags2.resize(numUV);
	tagsIn.resize(numUV);
	distance_from_intersection_to_node.resize(numUV);
	distance_between_nodes_at_IB.resize(numUV);
	uv.resize(numUV);

	//resize numP
	tagsP.resize(numP);
	tagsPOut.resize(numP);
	distance_from_u_to_body.resize(numP);
	distance_from_v_to_body.resize(numP);
	test.resize(numP);

	//cast cusp arrays
	tags_r		= thrust::raw_pointer_cast( &(tags[0]) );
	tagsOld_r	= thrust::raw_pointer_cast( &(tagsOld[0]) );
	tagsPOld_r	= thrust::raw_pointer_cast( &(tagsPOld[0]) );
	tags2_r		= thrust::raw_pointer_cast( &(tags2[0]) );
	tagsIn_r	= thrust::raw_pointer_cast( &(tagsIn[0]) );
	tagsP_r		= thrust::raw_pointer_cast( &(tagsP[0]) );
	tagsPOut_r	= thrust::raw_pointer_cast( &(tagsPOut[0]) );

	distance_from_intersection_to_node_r	= thrust::raw_pointer_cast( &(distance_from_intersection_to_node[0]) );
	distance_between_nodes_at_IB_r			= thrust::raw_pointer_cast( &(distance_between_nodes_at_IB[0]) );
	distance_from_u_to_body_r				= thrust::raw_pointer_cast( &(distance_from_u_to_body[0]) );
	distance_from_v_to_body_r				= thrust::raw_pointer_cast( &(distance_from_v_to_body[0]) );
	uv_r									= thrust::raw_pointer_cast( &(uv[0]) );
	test_r									= thrust::raw_pointer_cast( &(test[0]) );
}
