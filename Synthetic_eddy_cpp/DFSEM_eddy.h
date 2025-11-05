#pragma once

#include "Array.h"
#include "Eddy.h"

#include <numbers>
#include <cmath>
#include <algorithm>
//#include <math.h>
#include <cstdlib>

#define M_PI 3.14159265358979323746

struct DFSEM_eddy : public eddy {
public:
	double radius[3];
	//double r_norm;
	double rk[3];
	double r_magn2; // norm of r vector squared
	double K;

	// eigenvalues to evaluate range of random distribution function
	double eigvals[3];
	double alpha[3]; // can remove or replace epsilon...
	//double eigvect1[3];
	//double eigvect2[3];
	//double eigvect3[3];

	DFSEM_eddy()
		:eddy()
	{

	}

	DFSEM_eddy(const size_t& max_nodes)
		: eddy(max_nodes)
	{
		//rk.resize({ 3 });
		//radius.resize({ 3 });
		K = 1.0;
	}

	void reset(const double& new_x, const double& new_y, const double& new_z, const double new_r[3], const double& volume, const double& u, const double& dt, const Array<double>& y_plane, const Array<double>& z_plane, const double _alpha[3]) {
		position(0) = new_x;
		position(1) = new_y;
		position(2) = new_z;

		//radius = new_r;
		for (size_t m{ 0 }; m < 3; m++) {
			radius[m] = new_r[m];
		}

		//shape_scaling_factor = K * pow(16 * volume / (15 * M_PI * pow(new_r, 3)), 0.5);//volume_sqrt / pow(new_r, 3 / 2);
		shape_scaling_factor = 1.0;
		increment = u * dt;

		// insert nodes that need to be checked now
		// nodes.set(0);
		
		num_nodes = 0;
		double temp_dist{ 0 };
		double rad_max = std::max(radius[0], std::max(radius[1], radius[2]));
		for (size_t j{ 0 }; j < y_plane.shape[0]; j++) {
			for (size_t k{ 0 }; k < z_plane.shape[1]; k++) {
				temp_dist = pow(pow(y_plane(j, k) - new_y, 2) + pow(z_plane(j, k) - new_z, 2), 0.5);
				if (temp_dist < rad_max) {
				//if (temp_dist < radius) {
					nodes(num_nodes, 0) = j;
					nodes(num_nodes, 1) = k;
					nodes_pos(num_nodes, 0) = y_plane(j, k);
					nodes_pos(num_nodes, 1) = z_plane(j, k);
					num_nodes += 1;
				}
			}
		}
		
		for (size_t d{ 0 }; d < 3; d++) {
			alpha[d] = _alpha[d];
		}
	}

	void calc_r(const double& x, const double& y, const double& z) {
		rk[0] = (x - position(0)) / radius[0];
		rk[1] = (y - position(1)) / radius[1];
		rk[2] = (z - position(2)) / radius[2];

		theta_magn = sqrt(pow(rk[0], 2) + pow(rk[1], 2) + pow(rk[2], 2));

		r_magn2 = std::min(1.0, pow(theta_magn, 2));
	}
	/*
	const double shape_fn(const double& _radius) {
		if (theta_magn < 1) {
			return _radius * (1 - r_magn2);
		}
		else {
			return 0;
		}
	}*/

private:

public:

};

