#pragma once

#include <cstdlib>
#include <cmath>

#include "Array.h"
#include "Eddy.h"

struct XSEM_eddy : public eddy {
	double radius;


	XSEM_eddy(const size_t& max_nodes)
		: eddy(max_nodes)
	{

	}


};
