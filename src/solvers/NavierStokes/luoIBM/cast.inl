void castluo()
{
	/*
	double pressureStar_r = thrust::raw_pointer_cast ( &(pressureStar[0]) );
	double ustar_r = thrust::raw_pointer_cast ( &(ustar[0]) );
	double body_intercept_x_r = thrust::raw_pointer_cast ( &(body_intercept_x[0]) );
	double body_intercept_y_r = thrust::raw_pointer_cast ( &(body_intercept_y[0]) );
	double image_point_x_r = thrust::raw_pointer_cast ( &(image_point_x[0]) );
	double image_point_y_r = thrust::raw_pointer_cast ( &(image_point_y[0]) );
	double body_intercept_p_x_r = thrust::raw_pointer_cast ( &(body_intercept_p_x[0]) );
	double body_intercept_p_y_r = thrust::raw_pointer_cast ( &(body_intercept_p_y[0]) );
	double body_intercept_p_r = thrust::raw_pointer_cast ( &(body_intercept_p[0]) );
	double image_point_p_x_r = thrust::raw_pointer_cast ( &(image_point_p_x[0]) );
	double image_point_p_y_r = thrust::raw_pointer_cast ( &(image_point_p_y[0]) );
	double distance_from_intersection_to_node_r = thrust::raw_pointer_cast ( &(distance_from_intersection_to_node[0]) );
	double distance_between_nodes_at_IB_r = thrust::raw_pointer_cast ( &(distance_between_nodes_at_IB[0]) );
	double uv_r = thrust::raw_pointer_cast ( &(uv[0]) );
	double x1_ip_r = thrust::raw_pointer_cast ( &(x1_ip[0]) );
	double x2_ip_r = thrust::raw_pointer_cast ( &(x2_ip[0]) );
	double x3_ip_r = thrust::raw_pointer_cast ( &(x3_ip[0]) );
	double x4_ip_r = thrust::raw_pointer_cast ( &(x4_ip[0]) );
	double y1_ip_r = thrust::raw_pointer_cast ( &(y1_ip[0]) );
	double y2_ip_r = thrust::raw_pointer_cast ( &(y2_ip[0]) );
	double y3_ip_r = thrust::raw_pointer_cast ( &(y3_ip[0]) );
	double y4_ip_r = thrust::raw_pointer_cast ( &(y4_ip[0]) );
	double image_point_u_r = thrust::raw_pointer_cast ( &(image_point_u[0]) );
	double x1_r = thrust::raw_pointer_cast ( &(x1[0]) );
	double x2_r = thrust::raw_pointer_cast ( &(x2[0]) );
	double x3_r = thrust::raw_pointer_cast ( &(x3[0]) );
	double x4_r = thrust::raw_pointer_cast ( &(x4[0]) );
	double y1_r = thrust::raw_pointer_cast ( &(y1[0]) );
	double y2_r = thrust::raw_pointer_cast ( &(y2[0]) );
	double y3_r = thrust::raw_pointer_cast ( &(y3[0]) );
	double y4_r = thrust::raw_pointer_cast ( &(y4[0]) );
	double q1_r = thrust::raw_pointer_cast ( &([0]) );
	double q2_r = thrust::raw_pointer_cast ( &([0]) );
	double q3_r = thrust::raw_pointer_cast ( &([0]) );
	double q4_r = thrust::raw_pointer_cast ( &([0]) );
	double x1_p_r = thrust::raw_pointer_cast ( &(x1_p[0]) );
	double x2_p_r = thrust::raw_pointer_cast ( &(x2_p[0]) );
	double x3_p_r = thrust::raw_pointer_cast ( &(x3_p[0]) );
	double x4_p_r = thrust::raw_pointer_cast ( &(x4_p[0]) );
	double y1_p_r = thrust::raw_pointer_cast ( &(y1_p[0]) );
	double y2_p_r = thrust::raw_pointer_cast ( &(y2_p[0]) );
	double y3_p_r = thrust::raw_pointer_cast ( &(y3_p[0]) );
	double y4_p_r = thrust::raw_pointer_cast ( &(y4_p[0]) );
	double q1_p_r = thrust::raw_pointer_cast ( &(q1_p[0]) );
	double q2_p_r = thrust::raw_pointer_cast ( &(q2_p[0]) );
	double q3_p_r = thrust::raw_pointer_cast ( &(q3_p[0]) );
	double q4_p_r = thrust::raw_pointer_cast ( &(q4_p[0]) );
	double a0_r = thrust::raw_pointer_cast ( &(a0[0]) );
	double a1_r = thrust::raw_pointer_cast ( &(a1[0]) );
	double a2_r = thrust::raw_pointer_cast ( &(a2[0]) );
	double a3_r = thrust::raw_pointer_cast ( &(a3[0]) );
	double dudt_r = thrust::raw_pointer_cast ( &(dudt[0]) );
	double ududx_r = thrust::raw_pointer_cast ( &(ududx[0]) );
	double vdudy_r = thrust::raw_pointer_cast ( &(vdudy[0]) );
	double dvdt_r = thrust::raw_pointer_cast ( &(dvdt[0]) );
	double udvdx_r = thrust::raw_pointer_cast ( &(udvdx[0]) );
	double vdvdy_r = thrust::raw_pointer_cast ( &(vdvdy[0]) );

	double detA_r = thrust::raw_pointer_cast ( &(detA[0]) );
	double alpha_r = thrust::raw_pointer_cast ( &(alpha[0]) );
	
	double b11_r = thrust::raw_pointer_cast ( &(b11[0]) );
	double b12_r = thrust::raw_pointer_cast ( &(b12[0]) );
	double b13_r = thrust::raw_pointer_cast ( &(b13[0]) );
	double b14_r = thrust::raw_pointer_cast ( &(b14[0]) );
	double b21_r = thrust::raw_pointer_cast ( &(b21[0]) );
	double b22_r = thrust::raw_pointer_cast ( &(b22[0]) );
	double b23_r = thrust::raw_pointer_cast ( &(b23[0]) );
	double b24_r = thrust::raw_pointer_cast ( &(b24[0]) );
	double b31_r = thrust::raw_pointer_cast ( &(b31[0]) );
	double b32_r = thrust::raw_pointer_cast ( &(b32[0]) );
	double b33_r = thrust::raw_pointer_cast ( &(b33[0]) );
	double b34_r = thrust::raw_pointer_cast ( &(b34[0]) );
	double b41_r = thrust::raw_pointer_cast ( &(b41[0]) );
	double b42_r = thrust::raw_pointer_cast ( &(b42[0]) );
	double b43_r = thrust::raw_pointer_cast ( &(b43[0]) );
	double b44_r = thrust::raw_pointer_cast ( &(b44[0]) );
	

	int index1_r = thrust::raw_pointer_cast ( &(index1[0]) );
	int index2_r = thrust::raw_pointer_cast ( &(index2[0]) );
	int index3_r = thrust::raw_pointer_cast ( &(index3[0]) );
	int index4_r = thrust::raw_pointer_cast ( &(index4[0]) );
	
	int ghostTagsUV_r = thrust::raw_pointer_cast ( &(ghostTagsUV[0]) );
	int ghostTagsP_r = thrust::raw_pointer_cast ( &(ghostTagsP[0]) );
	int hybridTagsUV_r = thrust::raw_pointer_cast ( &(hybridTagsUV[0]) );
	int hybridTagsUV2_r = thrust::raw_pointer_cast ( &(hybridTagsUV2[0]) );
	int hybridTagsP_r = thrust::raw_pointer_cast ( &(hybridTagsP[0]) );
	int Tags_r = thrust::raw_pointer_cast ( &(Tags[0]) );
	
	bool q1flag_r = thrust::raw_pointer_cast ( &(q1flag[0]) );
	bool q2flag_r = thrust::raw_pointer_cast ( &(q2flag[0]) );
	bool q3flag_r = thrust::raw_pointer_cast ( &(q3flag[0]) );
	bool q4flag_r = thrust::raw_pointer_cast ( &(q4flag[0]) );
	
	double distance_from_u_to_body_r = thrust::raw_pointer_cast ( &(distance_from_u_to_body[0]) );
	double distance_from_v_to_body_r = thrust::raw_pointer_cast ( &(distance_from_v_to_body[0]) );*/
}