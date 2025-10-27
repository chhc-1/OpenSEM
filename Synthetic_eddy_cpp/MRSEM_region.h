#pragma once

#include "MRSEM_eddy.h"
#include "region.h"
#include "Array.h"
#include "interpolate.h"

#include <cstdlib>

class MRSEM_region : public region {
public:
	Array<MRSEM_eddy> eddies;


	MRSEM_region()
		: region()
	{

	}

	void increment_eddies() {

	}

	void increment_eddy() {

	}

private:
	
public:


};
