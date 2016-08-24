#include <solvers/NavierStokes/luo_base.h>


void luo_base::divergence()
{
	std::cout<<"Outputing divergence\n";
	int ip, iv, iu;
	int i_start = B.startI[0],
		j_start = B.startJ[0],
		width_i = B.numCellsX[0],
		height_j = B.numCellsY[0],
		i_end = i_start + width_i,
		j_end = j_start + height_j;
	std::ofstream test_output;
	std::string folder = (*paramDB)["inputs"]["caseFolder"].get<std::string>();
	std::stringstream out;
	out << folder << "/Divergence.csv";
	test_output<<	"v_n\t"
					"dy_n\t"
					"u_e\t"
					"dx_e\t"
					"v_s\t"
					"dy_s\t"
					"u_w\t"
					"dx_w\t"
					"x\t"
					"y\n";

	test_output.open(out.str().c_str());
	for (int J=j_start;  J<j_end;  J++)
	{
		for (int I=i_start;  I<i_end;  I++)
		{
			ip = J*nx + I;
			iu = J*(nx-1)+I;
			iv = J*nx + I + (nx-1)*ny;
			if (hybridTagsP[ip]>0)
			{
				test_output <<-u[iv] << "\t";
				test_output <<domInfo->dy[J] << "\t";
				test_output <<-u[iu] << "\t";
				test_output <<domInfo->dx[I] << "\t";
				test_output <<u[iv-nx] << "\t";
				test_output <<domInfo->dy[J] << "\t";
				test_output <<u[iu-1] << "\t";
				test_output <<domInfo->dx[I] << "\t";
				test_output <<domInfo->xv[I] << "\t";
				test_output <<domInfo->yu[J] << "\n";
			}
		}
	}
	test_output.close();
}

void luo_base::testInterpX()//flag split this into ghost node and hybrid node function
{
	std::cout<<"Outputing for interpolation of the u values\n";
	int iu;
	int i_start = B.startI[0],
		j_start = B.startJ[0],
		width_i = B.numCellsX[0],
		height_j = B.numCellsY[0],
		i_end = i_start + width_i,
		j_end = j_start + height_j;
	std::ofstream body_nodes;
	parameterDB  &db = *NavierStokesSolver::paramDB;
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();
	std::stringstream out;
	out << folder << "/interp_testX.csv";
	body_nodes.open(out.str().c_str());
	body_nodes << "BN_X1\t"
					"BN_Y1\t"
					"BN_X2\t"
					"BN_Y2\t"
					"GN_X\t"
					"GN_Y\t"
					"BI_X\t"
					"BI_Y\t"
					"IP_X\t"
					"IP_Y\t"
					"x1\t"
					"x2\t"
					"x3\t"
					"x4\t"
					"y1\t"
					"y2\t"
					"y3\t"
					"y4\t"
					"q1\t"
					"q2\t"
					"q3\t"
					"q4\t"
					"GN_U\t"
					"image_point_u\n"
			;
	for (int J=j_start;  J<j_end;  J++)
	{
		for (int I=i_start;  I<i_end;  I++)
		{
			iu = J*(nx-1) + I;
			if (ghostTagsUV[iu] >0)//for inside
			//if (hybridTagsUV[iu] >0)//for outside
			{
				body_nodes << x1_ip[iu]<<"\t";
				body_nodes << y1_ip[iu]<<"\t";
				body_nodes << x2_ip[iu]<<"\t";
				body_nodes << y2_ip[iu]<<"\t";
				body_nodes << domInfo->xu[I]<<"\t";
				body_nodes << domInfo->yu[J]<<"\t";
				body_nodes << body_intercept_x[iu] <<"\t";
				body_nodes << body_intercept_y[iu] <<"\t";
				body_nodes << image_point_x[iu] <<"\t";
				body_nodes << image_point_y[iu] <<"\t";
				body_nodes << x1[iu] <<"\t";
				body_nodes << x2[iu] <<"\t";
				body_nodes << x3[iu] <<"\t";
				body_nodes << x4[iu] <<"\t";
				body_nodes << y1[iu] <<"\t";
				body_nodes << y2[iu] <<"\t";
				body_nodes << y3[iu] <<"\t";
				body_nodes << y4[iu] <<"\t";
				body_nodes << q1[iu] <<"\t";
				body_nodes << q2[iu] <<"\t";
				body_nodes << q3[iu] <<"\t";
				body_nodes << q4[iu] <<"\t";
				body_nodes << u[iu] <<"\t";//inside
				//body_nodes << ustar[iu] <<"\t";//outside
				body_nodes << image_point_u[iu]<<"\n";
			}
		}
	}
	body_nodes.close();
}

void luo_base::testInterpY()
{
	std::cout<<"Outputing for interpolation of the v values\n";
	int iv;
	int i_start = B.startI[0],
		j_start = B.startJ[0],
		width_i = B.numCellsX[0],
		height_j = B.numCellsY[0],
		i_end = i_start + width_i,
		j_end = j_start + height_j;
	std::ofstream body_nodes;
	parameterDB  &db = *NavierStokesSolver::paramDB;
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();
	std::stringstream out;
	out << folder << "/interp_testY.csv";
	body_nodes.open(out.str().c_str());
	body_nodes <<	"BN_X1\t"
					"BN_Y1\t"
					"BN_X2\t"
					"BN_Y2\t"
					"GN_X\t"
					"GN_Y\t"
					"BI_X\t"
					"BI_Y\t"
					"IP_X\t"
					"IP_Y\t"
					"x1\t"
					"x2\t"
					"x3\t"
					"x4\t"
					"y1\t"
					"y2\t"
					"y3\t"
					"y4\t"
					"q1\t"
					"q2\t"
					"q3\t"
					"q4\t"
					"GN_U\t"
					"image_point_u\n";
	for (int J=j_start;  J<j_end;  J++)
	{
		for (int I=i_start;  I<i_end;  I++)
		{
			iv = J*nx + I  +  ny*(nx-1);
			if (ghostTagsUV[iv] >0)//for inside
			//if (hybridTagsUV[iv] >0)//for outside
			{
				//std::cout<<I<<"\t"<<J<<"\t"<<iv<<"\n";
				body_nodes << x1_ip[iv]<<"\t";
				body_nodes << y1_ip[iv]<<"\t";
				body_nodes << x2_ip[iv]<<"\t";
				body_nodes << y2_ip[iv]<<"\t";
				body_nodes << domInfo->xv[I]<<"\t";
				body_nodes << domInfo->yv[J]<<"\t";
				body_nodes << body_intercept_x[iv] <<"\t";
				body_nodes << body_intercept_y[iv] <<"\t";
				body_nodes << image_point_x[iv] <<"\t";
				body_nodes << image_point_y[iv] <<"\t";
				body_nodes << x1[iv] <<"\t";
				body_nodes << x2[iv] <<"\t";
				body_nodes << x3[iv] <<"\t";
				body_nodes << x4[iv] <<"\t";
				body_nodes << y1[iv] <<"\t";
				body_nodes << y2[iv] <<"\t";
				body_nodes << y3[iv] <<"\t";
				body_nodes << y4[iv] <<"\t";
				body_nodes << q1[iv] <<"\t";
				body_nodes << q2[iv] <<"\t";
				body_nodes << q3[iv] <<"\t";
				body_nodes << q4[iv] <<"\t";
				body_nodes << u[iv]  <<"\t";//inside
				//body_nodes << ustar[iv] <<"\t";//outside
				body_nodes << image_point_u[iv]<<"\n";
			}
		}
	}
	body_nodes.close();
}

void luo_base::testInterpP()
{
	std::cout<<"Outputing for interpolation of the p values\n";
	int ip;
	int i_start = B.startI[0],
		j_start = B.startJ[0],
		width_i = B.numCellsX[0],
		height_j = B.numCellsY[0],
		i_end = i_start + width_i,
		j_end = j_start + height_j;
	std::ofstream body_nodes;
	parameterDB  &db = *NavierStokesSolver::paramDB;
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();
	std::stringstream out;
	out << folder << "/interp_testP.csv";
	body_nodes.open(out.str().c_str());
	body_nodes <<	"BN_X1\t"
					"BN_Y1\t"
					"BN_X2\t"
					"BN_Y2\t"
					"GN_X\t"
					"GN_Y\t"
					"BI_X\t"
					"BI_Y\t"
					"IP_X\t"
					"IP_Y\t"
					"x1\t"
					"x2\t"
					"x3\t"
					"x4\t"
					"y1\t"
					"y2\t"
					"y3\t"
					"y4\t"
					"q1\t"
					"q2\t"
					"q3\t"
					"q4\t"
					"p*\t"
					"a0\t"
					"a1\t"
					"a2\t"
					"a3\t"
					"BI_p\t"
					"dudt\t"
					"ududx\t"
					"vdudy\t"
					"dvdt\t"
					"udvdx\t"
					"vdvdy\n"
					;
	for (int J=j_start;  J<j_end;  J++)
	{
		for (int I=i_start;  I<i_end;  I++)
		{
			ip = J*nx + I;
			if (ghostTagsP[ip] >0)//for inside
			//if (hybridTagsP[ip] >0)//for outside
			{
				//std::cout<<I<<"\t"<<J<<"\t"<<iv<<"\n";
				body_nodes << x1_ip_p[ip]<<"\t";
				body_nodes << y1_ip_p[ip]<<"\t";
				body_nodes << x2_ip_p[ip]<<"\t";
				body_nodes << y2_ip_p[ip]<<"\t";
				body_nodes << domInfo->xv[I]<<"\t";
				body_nodes << domInfo->yu[J]<<"\t";
				body_nodes << body_intercept_p_x[ip] <<"\t";
				body_nodes << body_intercept_p_y[ip] <<"\t";
				body_nodes << image_point_p_x[ip] <<"\t";
				body_nodes << image_point_p_y[ip] <<"\t";
				body_nodes << x1_p[ip] <<"\t";
				body_nodes << x2_p[ip] <<"\t";
				body_nodes << x3_p[ip] <<"\t";
				body_nodes << x4_p[ip] <<"\t";
				body_nodes << y1_p[ip] <<"\t";
				body_nodes << y2_p[ip] <<"\t";
				body_nodes << y3_p[ip] <<"\t";
				body_nodes << y4_p[ip] <<"\t";
				body_nodes << q1_p[ip] <<"\t";
				body_nodes << q2_p[ip] <<"\t";
				body_nodes << q3_p[ip] <<"\t";
				body_nodes << q4_p[ip] <<"\t";
				//body_nodes << pressureStar[ip] <<"\t";//outside
				body_nodes << pressure[ip] <<"\t";//inside
				body_nodes << a0[ip] <<"\t";
				body_nodes << a1[ip] <<"\t";
				body_nodes << a2[ip] <<"\t";
				body_nodes << a3[ip] <<"\t";
				body_nodes << body_intercept_p[ip] << "\n";
			}
		}
	}
	body_nodes.close();
}

void luo_base::testOutputX()
{
	int iu;
	int i_start = B.startI[0],
		j_start = B.startJ[0],
		width_i = B.numCellsX[0],
		height_j = B.numCellsY[0],
		i_end = i_start + width_i,
		j_end = j_start + height_j;
	std::ofstream body_nodes;
	parameterDB  &db = *NavierStokesSolver::paramDB;
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();
	std::stringstream out;
	out << folder << "/body_nodesX.csv";
	body_nodes.open(out.str().c_str());
	body_nodes << "BN_X1\tBN_Y1\tBN_X2\tBN_Y2\tGN_X\tGN_Y\tBI_X\tBI_Y\tIP_X\tIP_Y\n";
	for (int J=j_start;  J<j_end;  J++)
	{
		for (int I=i_start;  I<i_end;  I++)
		{
			iu = J*(nx-1) + I;
			if (hybridTagsUV[iu] >0) //for testing outside interpolation
			//if (ghostTagsUV[iu] >0) //for testing inside interpolation
			{
				std::cout<<iu<<std::endl;
				body_nodes << x1_ip[iu]<<"\t";
				body_nodes << y1_ip[iu]<<"\t";
				body_nodes << x2_ip[iu]<<"\t";
				body_nodes << y2_ip[iu]<<"\t";
				body_nodes << domInfo->xu[I]<<"\t";
				body_nodes << domInfo->yu[J]<<"\t";
				body_nodes << body_intercept_x[iu] <<"\t";
				body_nodes << body_intercept_y[iu] <<"\t";
				body_nodes << image_point_x[iu] <<"\t";
				body_nodes << image_point_y[iu] <<"\n";
			}
		}
	}
	body_nodes.close();
}

void luo_base::testOutputY()
{
	//test bi, ip
		int iv;
		int	i_start = B.startI[0],
			j_start = B.startJ[0],
			width_i = B.numCellsX[0],
			height_j = B.numCellsY[0],
			i_end = i_start + width_i,
			j_end = j_start + height_j;
		std::ofstream body_nodes;
		parameterDB  &db = *NavierStokesSolver::paramDB;
		std::string folder = db["inputs"]["caseFolder"].get<std::string>();
		std::stringstream out;
		out << folder << "/body_nodesY.csv";
		body_nodes.open(out.str().c_str());
		body_nodes << "x1\ty1\tx2\ty2\tg_x\tg_y\tbi_x\tbi_y\tip_x\tip_y\n";
		arrayprint(ghostTagsUV,"gnuvY","y",-1);
		arrayprint(ghostTagsUV,"gnuvX","x",-1);
		for (int J=j_start;  J<j_end;  J++)
		{
			for (int I=i_start;  I<i_end;  I++)
			{
				iv = J*(nx) + I + ny*(nx-1);
				if (hybridTagsUV[iv] >0) //for testing outside interpolation
				//if (ghostTagsUV[iv] >0) //for testing inside interpolation
				{
					std::cout<<iv<<std::endl;
					body_nodes << x1_ip[iv]<<"\t";
					body_nodes << y1_ip[iv]<<"\t";
					body_nodes << x2_ip[iv]<<"\t";
					body_nodes << y2_ip[iv]<<"\t";
					body_nodes << domInfo->xv[I]<<"\t";
					body_nodes << domInfo->yv[J]<<"\t";
					body_nodes << body_intercept_x[iv] <<"\t";
					body_nodes << body_intercept_y[iv] <<"\t";
					body_nodes << image_point_x[iv] <<"\t";
					body_nodes << image_point_y[iv] <<"\n";
				}
			}
		}
		body_nodes.close();
}

void luo_base::testForce_p()
{
	std::cout<<"Outputing for interpolation of the pressure force values\n";
	std::ofstream body_nodes;
	parameterDB  &db = *NavierStokesSolver::paramDB;
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();
	std::stringstream out;
	out << folder << "/interp_test_force.csv";
	body_nodes.open(out.str().c_str());
	body_nodes <<	"bx\t"
					"by\t"
					"body node p\t"
					"close\t"
					"close2\n"
					;
	for (int i=0;  i<B.totalPoints;  i++)
	{
		body_nodes << B.x[i] <<"\t";
		body_nodes << B.y[i] <<"\t";
		body_nodes << B.force_pressure[i] <<"\n";
	}
	body_nodes.close();
}

void luo_base::testForce_dudn()
{
	std::cout<<"Outputing for interpolation of the dudn force values\n";
	std::ofstream body_nodes;
	parameterDB  &db = *NavierStokesSolver::paramDB;
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();
	std::stringstream out;
	out << folder << "/interp_test_force_dudn.csv";
	body_nodes.open(out.str().c_str());
	body_nodes <<	"bx\t"
					"by\t"
					"body node dudn\t"
					"pressure component\t"
					"dudn component\n"
					;
	for (int i=0;  i<B.totalPoints;  i++)
	{
		body_nodes << B.x[i] <<"\t";
		body_nodes << B.y[i] <<"\t";
		body_nodes << B.force_dudn[i] <<"\t";
		body_nodes << B.force_x[i] <<"\t";
		body_nodes << B.force_y[i] <<"\n";
	}
	body_nodes.close();
}
