#pragma once

/**
 * \enum  bcType
 * \brief Specifies the type of boundary condition.
 * \not used, outdated from cuibm
 */
enum bcType
{
	DIRICHLET,  ///< Dirichlet boundary condition
	NEUMANN,    ///< Neumann boundary condition		//flag
	CONVECTIVE, ///< convective boundary condition
	PERIODIC,   ///< periodic boundary condition		//flag
	SPECIAL							//flag
};
