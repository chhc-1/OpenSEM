#pragma once

#include "Array.h"
#include "Eddy.h"

#include <climits>
#include <malloc.h>
#include <initializer_list>
#include <cassert>
#include <cmath>
#include <exception>
#include <random>
#include <iostream>
#include <map>
#include <fstream>
#include <string>

struct ISEM1_eddy : public eddy {
public:
	double radius;

	ISEM1_eddy()
		: eddy()
	{

	}

	ISEM1_eddy(const size_t& max_nodes)
		: eddy(max_nodes)
	{

	}

	void reset(const double& new_x, const double& new_y, const double& new_z, const double& new_r, const double& volume_sqrt, const double& u, const double& dt, const Array<double>& y_plane, const Array<double>& z_plane, const Array<int>& _eps) {
		// reset eddy position
		position(0) = new_x;
		position(1) = new_y;
		position(2) = new_z;
		radius = new_r;
		shape_scaling_factor = volume_sqrt / pow(new_r, 1.5);
		increment = u * dt;

		// insert nodes that need to be checked now
		// nodes.set(0);
		num_nodes = 0;
		double temp_dist{ 0 };
		for (size_t j{ 0 }; j < y_plane.shape[0]; j++) {
			for (size_t k{ 0 }; k < z_plane.shape[1]; k++) {
				temp_dist = pow(pow(y_plane(j, k) - new_y, 2) + pow(z_plane(j, k) - new_z, 2), 0.5);
				if (temp_dist < radius) {
					nodes(num_nodes, 0) = j;
					nodes(num_nodes, 1) = k;
					nodes_pos(num_nodes, 0) = y_plane(j, k);
					nodes_pos(num_nodes, 1) = z_plane(j, k);
					num_nodes += 1;
				}
			}
		}

		for (size_t j{ 0 }; j < 3; j++) {
			epsilon(j) = _eps(j);
		}
	}

	double shape_fn(const double& ref_pos, const double& eddy_pos) {
		theta = (ref_pos - eddy_pos) / radius;

		if (std::abs(theta) < 1) {
			return pow(1.5, 0.5) * (1 - std::abs(theta));
		}
		else {
			return 0;
		}

	}

private:

public:

};